%%-----------------------SOLAR PANEL SPECS----------------%%
Pmp_ref = 49;
Imp_ref = 2.88;
Vmp_ref = 17;
Isc_ref = 3.11;
Voc = 21.8;
a_ref = 1.2;
G_ref = 1000;
Iph_ref = Isc_ref;
Io_ref = Isc_ref*exp(-Voc/a_ref);
Rp_ref = 3244.6;%%calculated from the previous code
Rs_ref = 0.571;%%calculatedfrom the previous code

%%-------------Hydrogen generation Electrolyzer(PHOEBUS) SPECS-----------
global Rt Ct Ta Tcw T_in eta_c Ccw F tao_t U_tn rUI sUI tUI UrevUI
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

%%------------FUEL CELL SPECS-----%%
U_FCo = 33.18; %V
e1 = -0.013; %V.C^*1
e2 = -1.57; %V.C^-1
I_FCo = 8.798; %A
R_FC = -2.04; %ohm.C^-1
N_FC = 35;
C_H2 = 2.39; %Ah/l_conversion coefficient for hydrogen
eta_FC = 0.7; %utilization factor


%%--------------DATA-----------%%
load('ts_2.mat');

%%----GENERATING MEASURED DATA-------------%%
%%assuming no error in electrolyzer's voltage
MT = T + 2*randn(m,1) ;%% normal noise of standard deviation 2
MI_electro = I_electro + 3*randn(m,1); %% normal noise of standard deviation 3
MP_electro = V_electro(1:m).*MI_electro/1000; %%(kW)

%%----------PARAMETERS FOR STATE ESTIMATION-------%%


                   %--ELECTROLYZER
Ts = 2/3600;%%sampling time for current and voltage
T_start = 1  %in hr
T_end = 2    % in hr
m_start = floor(T_start/Ts);
m_end = floor(T_end/Ts);
P = [0.01 , 0 ; 0, 0.01];
Q = 0.1;
R = [0.1,0;0,0.1];
R_f = 0.5;
H_f = [0,1];
H = [1,0;0,1];
EI_elec = zeros(m,1);
ET = zeros(m,1);

tic();
%%---------STATE ESTIMATION--------------%%

for i=1:m
i
             %%-----Hydrogen generation Electrolyzer(PHOEBUS)-----%%
%% assume that we are getting correct values of voltage
 if(i==1)
  ET(i) = MT(i);
  UI = @(x)polyval(UrevUI,ET(i)) + polyval(rUI,ET(i))*x + polyval(sUI,ET(i))*log(polyval(tUI,ET(i))*x + 1) - V_electro(i);
  EI_elec(i) = fzero(UI,45)*2.5;
  x_elec = [ET(i);EI_elec(i)];
  EP_elec(i) = EI_elec(i)*V_electro(i)/1000;

 else 
  fx = -1/29 - 1000000*0.6*(0.004186/625)*(1 - exp(-0.00000036*(9.955+0.0091935*EI_elec(i-1))/(0.6*0.004186)));
  fz = 3.6*21*(V_electro(i-1) -1.482)/625 - 3.6*(0.0091935/625)*(ET(i-1) - 14.5)*exp(-0.00000036*(9.955 + 0.0091935*EI_elec(i-1))/(0.6*0.004186));
  gx = UrevUI(1) + rUI(1)*(EI_elec(i-1)/2.5) + (polyval(sUI,ET(i-1))/(polyval(tUI,ET(i-1))*(EI_elec(i-1)/2.5) + 1))*(2*tUI(1)*ET(i-1) + tUI(2))*(EI_elec(i-1)/2.5) + sUI(1)*log(polyval(tUI,ET(i-1))*(EI_elec(i-1)/2.5) + 1);
  gz = polyval(rUI,ET(i-1))/2.5 + polyval(sUI,ET(i-1))*((polyval(tUI,ET(i-1)))/(polyval(tUI,ET(i-1))*EI_elec(i-1) + 2.5));
  
  gx_ = -inv(gz)*gx*fx;
  gz_ = -inv(gz)*gx*fz;
  
  %decoupled prediction
  %phy = [exp(fx*Ts), (exp(fx*Ts)-1)*fz*inv(fx); (exp(gz_*Ts)-1)*gx_*inv(gz_), exp(gz_*Ts)]; 

  %coupled prediction
  %A_aug = [fx,fz ; gx_,gz_];
  %phy = exp(A_aug*Ts);
  
  %using C2D transformation
  A_aug = [fx,fz ; gx_,gz_];
  sys = ss(A_aug,[0;0],H,0);
  sysd = c2d(sys,Ts);
  [phy,B_,C_,D_] = ssdata(sysd);
  %%results of decoupled and C2D transformation are same and better than coupled
  
  gamma = [1; -inv(gz)*gx];
  
  %%Predication Step
  M = diag([1,0]);
  y0 = zeros(2,1);
  y0(1) = ET(i-1);
  y0(2) = EI_elec(i-1);
  options=odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-2);
  tspan = [t(i-1) t(i)];
  [t_,y]=ode23t(@(t,y) ElecDae2(t,y,V_electro(i-1)),tspan,y0,options);
  
  ET(i) = y(end,1);
  EI_elec(i) = y(end,2);
  x_elec = [ET(i);EI_elec(i)];
   
  P = phy*P*(phy') + gamma*Q*(gamma');
  
  if(i>=m_start && i <=m_end )         %shutdown time
      
    L_f = P*H_f'*inv(H_f*P*H_f' + R_f);
    x_elec = x_elec + L_f*([MI_electro(i)] - H_f*x_elec);
    P = (eye(2) - L_f*H_f)*P;
    
  else            
   L = P*H'*inv(H*P*H' + R);
   x_elec = x_elec + L*([MT(i);MI_electro(i)] - H*x_elec);
   P = (eye(2) - L*H)*P;
  
  end
  %%update
  
  ET(i) = x_elec(1);  
  UI = @(x)polyval(UrevUI,ET(i)) + polyval(rUI,ET(i))*x + polyval(sUI,ET(i))*log(polyval(tUI,ET(i))*x + 1) - V_electro(i);
  EI_elec(i) = fzero(UI, 45)*2.5;
  
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

%plot(t,T,t,MT,'r*',t,ET,'-g')