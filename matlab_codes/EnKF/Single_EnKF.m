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
load('Ts_2.mat')

%%----GENERATING MEASURED DATA-------------%%
%%assuming no error in electrolyzer's voltage
MT = T + 2*randn(m,1) ;%% normal noise of standard deviation 0.2
MI_electro = I_electro + 3*randn(m,1); %% normal noise of standard deviation 0.1
MP_electro = V_electro(1:m).*MI_electro/1000; %%(kW)
%%----------PARAMETERS FOR STATE ESTIMATION-------%%
                   %--ELECTROLYZER

Ts = 2/3600;%%sampling time(in hrs)
H = [1,0;0,1];
P = [0.01 , 0 ; 0, 0.01];
Q = 0.05;
R = [0.1,0;0,0.1];
EI_elec = zeros(m,1);
ET = zeros(m,1);
Np = 20

tic();
%%---------STATE ESTIMATION--------------%%

for i=1:m
             %%-----Hydrogen generation Electrolyzer(PHOEBUS)-----%%
%% assume that we are getting correct values of voltage
 if(i==1)
  ET(i) = MT(i);
  UI = @(x)polyval(UrevUI,ET(i)) + polyval(rUI,ET(i))*x + polyval(sUI,ET(i))*log(polyval(tUI,ET(i))*x + 1) - V_electro(i);
  EI_elec(i) = fzero(UI,45)*2.5;
  EP_elec(i) = EI_elec(i)*V_electro(i)/1000;
  x_elec = [ET(i);EI_elec(i)];
  S = chol(P);
  
  % Initial Particles
  for j=1:Np
   yk = mvnrnd([0;0],eye(2));
   x_elecp = x_elec + S*yk';
   ET_p(j) = x_elecp(1);
   EI_elecp(j) = x_elecp(2);
  end 
  
 else
  i
  v1 = ET_p;
  v2 = EI_elecp;
  
                 %Forecast step%  
  parfor j=1:Np
   %yk = mvnrnd([0;0],eye(2));
   %x_elecp = x_elec + S*yk';
   %EI_elecp(j) = x_elecp(2);
   %ET_p(j) = x_elecp(1);

   %%old method
   %fx = -1/29 - 1000000*0.6*(0.004186/625)*(1 - exp(-0.00000036*(9.955+0.0091935*x_elecp(2))/(0.6*0.004186)));
   %fz = 3.6*21*(V_electro(i-1) -1.482)/625 - 3.6*(0.0091935/625)*(x_elecp(1) - 14.5)*exp(-0.00000036*(9.955 + 0.0091935*x_elec(2))/(0.6*0.004186));
   %Xf = exp(fx*Ts);
   %Zf = (Xf-1)*fz*inv(fx);
   %ET_pf(j) = Xf*ET_p(j) + Zf*EI_elecp(j) +  Q*randn();
   %UI = @(x)polyval(UrevUI,ET_pf(j)) + polyval(rUI,ET_pf(j))*x + polyval(sUI,ET_pf(j))*log(polyval(tUI,ET_pf(j))*x + 1) - V_electro(i-1);
   %EI_elecpf(j) = (fzero(UI, 45))*2.5;
   
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
  end
  ET(i) = (1/Np)*sum(ET_p); 
  EI_elec(i) = (1/Np)*sum(EI_elecp); 
  x_elec = [ET(i);EI_elec(i)];
  
  P = zeros(2,2);
  for j=1:Np
   e = [ET_p(j);EI_elecp(j)] - x_elec;
   P = P + e*e';
  end
  P = P./(Np-1);

  %%kalman Gain
   L = P*H'*inv(H*P*H' + R);

  %%update
  
   for j=1:Np
    x_elecp = [ET_p(j);EI_elecp(j)] + L*( [MT(i);MI_electro(i)] +  mvnrnd([0;0],R)' - H*[ET_p(j);EI_elecp(j)]) ;
    ET_p(j) = x_elecp(1);
    UI = @(x)polyval(UrevUI,ET_p(j)) + polyval(rUI,ET_p(j))*x + polyval(sUI,ET_p(j))*log(polyval(tUI,ET_p(j))*x + 1) - V_electro(i);
    EI_elecp(j) = (fzero(UI, 45))*2.5;
    %xel0 = [ET(i); EI_elec(i)];
    %y_elec = H*[ET_pf(j);EI_elecpf(j)] + mvnrnd([0;0],R)';
    %phi = @(xel)([xel(1);xel(2)] - [ET_pf(j);EI_elecpf(j)])'*inv(P)*([xel(1);xel(2)] - [ET_pf(j);EI_elecpf(j)])  +  ([MT(i);MI_electro(i)] - y_elec)'*inv(R)*([MT(i);MI_electro(i)] - y_elec);
    %constr = @(xel)polyval(UrevUI,xel(1)) + polyval(rUI,xel(1))*xel(2) + polyval(sUI,xel(1))*log(polyval(tUI,xel(1))*xel(2) + 1) - V_electro(i-1);   
    
    %[xel, obj, info, iter, nf, laambda] = sqp(xel0, phi, constr, [],[],[],500,10e-5);
    %ET_p(j) = xel(1);
    %EI_elecp(j) = xel(2);
   end
   ET(i) = (1/Np)*sum(ET_p);
   UI = @(x)polyval(UrevUI,ET(i)) + polyval(rUI,ET(i))*x + polyval(sUI,ET(i))*log(polyval(tUI,ET(i))*x + 1) - V_electro(i);
   EI_elec(i) = (fzero(UI, 45))*2.5;
   x_elec = [ET(i);EI_elec(i)];
   
   P = zeros(2,2);
   for j=1:Np
    e = [ET_p(j);EI_elecp(j)] - x_elec;
    P = P + e*e';
   end
   P = P./(Np-1);
  
   EP_elec(i) = EI_elec(i)*V_electro(i)/1000;
   TP_elec(i) = 21*EP_elec(i);
 end

end

%% -------------MSE calculations ----------------%%
RMSE_TM = (immse(T,MT))^0.5
RMSE_TE = (immse(T,ET))^0.5
RMSE_IMelec = (immse(I_electro,MI_electro))^0.5
RMSE_IEelec = (immse(I_electro,EI_elec))^0.5
%RMSE_PEelec = (immse(P_electro,EP_elec'))^0.5

elapsed_time = toc()
Computation_time = cputime()
