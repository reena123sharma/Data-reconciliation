%%-------------Hydrogen generation Electrolyzer(PHOEBUS) SPECS-----------
Rt = 0.167; %CW^-1
Ct = 625; %KJC^-1
Ta = 20; %ambient temp
Tcw = 14.5; %cooling water inlet temp.
T_in = 51.7; %initial temp of electrolyzer 
eta_c = 21; %no. of cells in series
Ccw = 0.004186; %KJ.C^-1.g^-1
F = 0.6; %Nm^3/h 
tao_t = 29; %h
U_tn = 1.482; %thermal Voltage

rUI = [ 1.1211e-005  -3.0071e-004];
sUI = [-0.0024638   0.2772353];
tUI = [ -3.9527e-004  5.6183e-002  -1.4421e+000];
UrevUI = [-8.182e-004, 1.249];

%%--------------DATA-----------%%
load('Ts_1_2.mat');

%%----GENERATING MEASURED DATA-------------%%
%%assuming no error in electrolyzer's voltage
MT = T + 2*randn(m_T,1) ;%% normal noise of standard deviation 0.2
MI_electro = I_electro + 3*randn(m,1); %% normal noise of standard deviation 0.1
MP_electro = V_electro(1:m).*MI_electro/1000; %%(kW)

%%----------PARAMETERS FOR STATE ESTIMATION-------%%

                   %--ELECTROLYZER

Ts =1/3600;%%minorsampling time(in hrs)
Ts_T = 2/3600;%%minorsampling time
mRR = Ts_T/Ts;
Np = 20
H = [1,0;0,1];
H_f = [0,1];
P = [0.5 , 0 ; 0, 0.1];
S = chol(P);
Q = 0.05;
R = [2,0;0,3];
R_f = 2;
EI_elec = zeros(m,1);
ET = zeros(m,1);
wp = zeros(Np,1);
iT = 1;

tic();
%%---------STATE ESTIMATION--------------%%

for i=1:m

             %%-----Hydrogen generation Electrolyzer(PHOEBUS)-----%%
%% assume that we are getting correct values of voltage
 if(i==1)
  ET(i) = MT(i);
  UI = @(x)polyval(UrevUI,ET(i)) + polyval(rUI,ET(i))*x + polyval(sUI,ET(i))*log(polyval(tUI,ET(i))*x + 1) - V_electro(i);
  EI_elec(i) = (fzero(UI, 45))*2.5;
  x_elec = [ET(1);EI_elec(1)];
% Initial Particles
  parfor j=1:Np
   yk = mvnrnd([0;0],eye(2));
   x_elecp(j,:) = x_elec + S*yk';
  end
   ET_p = x_elecp(:,1);
   EI_elecp = x_elecp(:,2);
   
 else 
  if(mod((i-1),mRR)~=0)          %minor sampling time
  i
  v1 = ET_p;
  v2 = EI_elecp;
%Prediction step%  
  parfor j=1:Np
   M = diag([1,0]);
   y0 = zeros(2,1);
   y0(1) = v1(j);
   y0(2) = v2(j);
   options=odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-2);
   tspan = [t(i-1) t(i)];
   [t_,y]=ode23t(@(t,y) ElecDae2(t,y,V_electro(i-1)),tspan,y0,options);
   
   ET_p(j) = y(end,1) + Q*randn();
   T_copy = ET_p(j);
   UI = @(x)polyval(UrevUI,T_copy) + polyval(rUI,T_copy)*x + polyval(sUI,T_copy)*log(polyval(tUI,T_copy)*x + 1) - V_electro(i-1);
   EI_elecp(j) = (fzero(UI, 45))*2.5;
   x_elecp(j,:) = [ET_p(j),EI_elecp(j)];
  
% Importance weight calculation
   wp1(j) = (1/(det(R_f))^0.5)*exp( (-([MI_electro(i)] - H_f*x_elecp(j,:)')'*inv(R_f)*([MI_electro(i)] - H_f*x_elecp(j,:)'))/2 );
  end
  
  elseif(mod((i-1),mRR)==0)               %%major sampling time
   
   iT = ((i-1)/mRR) + 1;
   i
   v1 = ET_p;
   v2 = EI_elecp;
%Prediction step%  
  parfor j=1:Np
   M = diag([1,0]);
   y0 = zeros(2,1);
   y0(1) = v1(j);
   y0(2) = v2(j);
   options=odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-2);
   tspan = [t(i-1) t(i)];
   [t_,y]=ode23t(@(t,y) ElecDae2(t,y,V_electro(i-1)),tspan,y0,options);
   
   ET_p(j) = y(end,1) + Q*randn();
   T_copy = ET_p(j);
   UI = @(x)polyval(UrevUI,T_copy) + polyval(rUI,T_copy)*x + polyval(sUI,T_copy)*log(polyval(tUI,T_copy)*x + 1) - V_electro(i-1);
   EI_elecp(j) = (fzero(UI, 45))*2.5;
   x_elecp(j,:) = [ET_p(j),EI_elecp(j)];
  
% Importance weight calculation
   wp1(j) = (1/(det(R))^0.5)*exp( (-([MT(iT);MI_electro(i)] - H*x_elecp(j,:)')'*inv(R)*([MT(iT);MI_electro(i)] - H*x_elecp(j,:)'))/2 );
 
  end
 end 
   
  sum_wp = sum(wp1);
  for j = 1 : Np
   wp(j) = wp1(j)/sum_wp;
  end

%State estimate
 %%Regularized particle filter resampling
   mu = mean(x_elecp(:,1));
   
   S_r = zeros(1,1);  %%currnet prediction covariance 
   for j=1:Np
     S_r = S_r + (x_elecp(j,1)-mu)'*(x_elecp(j,1)-mu);
   end
   S_r = (1/(Np-1))*S_r;
   A_r = chol(S_r);     %%square root factorization
   
   v_1 = 2;  %% dimension of state vector is 2
   
   h = (1/2)*((8*5*(2*sqrt(pi))/(v_1))^(1/5))*Np^(-1/5); %%optimal kernal bandwidth(tuning parameter)
   
 %%random sampling
  w = abs(mu-MT(iT));  %%width
  if w>1
      w=1;
  end
  
  for j=1:5*Np    %%
    ET_pf(j) = (mu -0.5*w) + rand(1,1)*2*(0.5*w);
    prob(j) = pdf(ET_pf(j),Np,wp,ET_p,A_r,h);
  end  
   
   prob = prob./sum(prob);
   
   ET_p = randsample(ET_pf,Np,true,prob);
   
   for j = 1:Np
     UI = @(x)polyval(UrevUI,ET_p(j)) + polyval(rUI,ET_p(j))*x + polyval(sUI,ET_p(j))*log(polyval(tUI,ET_p(j))*x + 1) - V_electro(i);
     EI_elecp(j) = (fzero(UI, 45))*2.5;
     x_elecp(j,:) = [ET_p(j),EI_elecp(j)];
   end  
 
  ET(i) = mean(ET_p);
  UI = @(x)polyval(UrevUI,ET(i)) + polyval(rUI,ET(i))*x + polyval(sUI,ET(i))*log(polyval(tUI,ET(i))*x + 1) - V_electro(i);
  EI_elec(i) = (fzero(UI, 45))*2.5;
  
end

 EP_elec(i) = EI_elec(i)*V_electro(i)/1000;
 TP_elec(i) = 21*EP_elec(i);

end

%% -------------MSE calculations ----------------%%
%---------generating the temp. data for comparison
ET_modi = zeros(m_T,1);
iT = 0;
for i = 1:m
 if (i==1)
   ET_modi(i) = ET(i);
 elseif(mod((i-1),mRR) ==0)
   iT = ((i-1)/mRR) + 1;
   ET_modi(iT) = ET(i);
 end
end
RMSE_TM = (immse(T,MT))^0.5
RMSE_TE = (immse(T,ET_modi))^0.5
RMSE_IMelec = (immse(I_electro,MI_electro))^0.5
RMSE_IEelec = (immse(I_electro,EI_elec))^0.5
RMSE_PEelec = (immse(P_electro,EP_elec'))^0.5

elapsed_time = toc()
Computation_time = cputime()