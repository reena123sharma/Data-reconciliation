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
% t = xlsread('data_central','Sheet3','A3:A1143');
% G = xlsread('data_central','Sheet3','B3:B1143');
% V_solar = xlsread('data_central','Sheet3','C3:C1143');
% I_solar = xlsread('data_central','Sheet3','D3:D1143');
% P_solar = V_solar.*I_solar;
% 
% V_electro = xlsread('data_central','Sheet3','E3:E1143');
% I_electro = xlsread('data_central','Sheet3','F3:F1143');
% T = xlsread('data_central','Sheet3','G3:G1143');
% P_electro = V_electro.*I_electro/1000;
% 
% m = length(t);
% 
% U_FC = xlsread('data_central','Sheet3','H3:H1143');
% I_FC = xlsread('data_central','Sheet3','I3:I1143');
% Hyd = 0.1*ones(m,1); % kmol/hr
% P_FC = U_FC.*I_FC/1000;


%%----GENERATING MEASURED DATA-------------%%
MV_solar = V_solar + 0.5*randn(m,1) ;%% normal noise of standard deviation 0.2
MI_solar = I_solar + 0.3*randn(m,1); %% normal noise of standard deviation 0.1
MP_solar = MV_solar.*MI_solar;

%%assuming no error in electrolyzer's voltage
MT = T + 2*randn(m,1) ;%% normal noise of standard deviation 0.2
MI_electro = I_electro + 3*randn(m,1); %% normal noise of standard deviation 0.1
MP_electro = V_electro(1:m).*MI_electro/1000; %%(kW)

MHyd = Hyd +  0.01*randn(m,1);
MU_FC = U_FC + randn(m,1); 
MP_FC = I_FC.*MU_FC/1000;  %%(kW) 

%%----------PARAMETERS FOR STATE ESTIMATION-------%%

                    %--SOLAR

Var_solar = [0.25,0;0 , 0.01];
A_solar = zeros(2,1);
Iph = G.*Iph_ref/G_ref; %%Assume the cell temp. = ref cell temp. = 25C
Io = Io_ref;
a = a_ref;
Rs = Rs_ref;%%assumption
Rp = Rp_ref*G_ref./G;

                   %--ELECTROLYZER

Ts = 2/3600;%%sampling time
H = [1,0;0,1];
P = [0.01 , 0 ; 0, 0.05];
Q = 0.05;
R = [1,0;0,2];
EI_elec = zeros(m,1);
ET = zeros(m,1);
Np = 20
options = optimoptions('fsolve','MaxIterations',1000);
           %--FUEL CELL

T_FC = 298;                       
Var_FC = [0.0001,0;0, 1];
A_FC = zeros(2,1);

tic();
%%---------STATE ESTIMATION--------------%%

for i=1:m

               %%----SOLAR------%%
I = MI_solar(i);
V = MV_solar(i);
y_cap = [I;V];

 for j=1:4
  V = y_cap(2);
  I = y_cap(1);
  A_solar = [ -(Io*Rs/a)*exp((V+I*Rs)/a)-Rs/Rp(i)-1 ,  -(Io/a)*exp((V + I*Rs)/a)-1/Rp(i) ];
  f = Iph(i) - Io*(exp((V+I*Rs)/a) -1) -(V+I*Rs)/Rp(i) - I;
  y = y_cap - Var_solar*(A_solar')*(inv(A_solar*Var_solar*A_solar'))*f;
  y_cap = y;
 end
EP_solar(i) = y(1)*y(2);
EI_solar(i) = y(1);
EV_solar(i) = y(2);

             %%-----Hydrogen generation Electrolyzer(PHOEBUS)-----%%
%% assume that we are getting correct values of voltage
 if(i==1)
  x_elec = [MT(1);MI_electro(1)];
  EI_elec(i) = x_elec(2);
  ET(i) = x_elec(1);
  EP_elec(i) = EI_elec(i)*V_electro(i)/1000;
  
 else
  i   
  x_elec = [ET(i-1);EI_elec(i-1)];
  S = chol(P);

                 %Forecast step%  
   j = 0;              
  while(j<=Np)
   j = j+1;   
   yk = mvnrnd([0;0],eye(2));
   x_elecp = x_elec + S*yk';
   ET_p(j) = x_elecp(1);
   if((ET_p(j)<=49) || (ET_p(j)>=65))
       j= j-1
   continue;
   end    
   UI = @(x)polyval(UrevUI,ET_p(j)) + polyval(rUI,ET_p(j))*x + polyval(sUI,ET_p(j))*log(polyval(tUI,ET_p(j))*x + 1) - V_electro(i-1);
   EI_elecp(j) = (fsolve(UI,45,options))*2.5;
  end
  v1 = ET_p;
  v2 = EI_elecp;
  
  for j=1:Np
   A = -1/29 - 1000000*0.6*(0.004186/625)*(1 - exp(-0.00000036*(9.955+0.0091935*v2(j))/(0.6*0.004186)));
   B = 3.6*21*(V_electro(i-1) -1.482)/625 - 3.6*(0.0091935/625)*(v1(j) - 14.5)*exp(-0.00000036*(9.955 + 0.0091935*v2(j))/(0.6*0.004186));
     
   Gf = exp(A*Ts);
   Hf = (Gf -1)*B*inv(A);
   ET_pf(j) = Gf*v1(j) + Hf*v2(j) + Q*randn();
   T_c = ET_pf(j);
   UI = @(x)polyval(UrevUI,T_c) + polyval(rUI,T_c)*x + polyval(sUI,T_c)*log(polyval(tUI,T_c)*x + 1) - V_electro(i-1);
   EI_elecpf(j) = (fsolve(UI, 45,options))*2.5;
  end
  ET(i) = mean(ET_pf); 
  EI_elec(i) = mean(EI_elecpf); 
  x_elec = [ET(i);EI_elec(i)];
  
  P = zeros(2,2);
  for j=1:Np
   e = [ET_p(j);EI_elecp(j)] - [ET(i);EI_elec(i)];
   P = P + e*e';
  end
  P = P./(Np-1);

  %%kalman Gain
   L = P*H'*inv(H*P*H' + R);

  %%update
   for j=1:Np
    x_elecp = [ET_pf(j);EI_elecpf(j)] + L*( [MT(i);MI_electro(i)] +  mvnrnd([0;0],R)' - [ET_pf(j);EI_elecpf(j)]) ;
    ET_p(j) = x_elecp(1);
   %ET_p(j) = ET_pf(j) + L(1)*(MT(i) + 0.5*randn() - ET_pf(j)); 
    UI = @(x)polyval(UrevUI,ET_p(j)) + polyval(rUI,ET_p(j))*x + polyval(sUI,ET_p(j))*log(polyval(tUI,ET_p(j))*x + 1) - V_electro(i);  
    EI_elecp(j) = (fsolve(UI, 45,options))*2.5;
    %xel0 = [ET(i); EI_elec(i)];
    %y_elec = H*[ET_pf(j);EI_elecpf(j)] + mvnrnd([0;0],R)';
    %phi = @(xel)([xel(1);xel(2)] - [ET_pf(j);EI_elecpf(j)])'*inv(P)*([xel(1);xel(2)] - [ET_pf(j);EI_elecpf(j)])  +  ([MT(i);MI_electro(i)] - y_elec)'*inv(R)*([MT(i);MI_electro(i)] - y_elec);
    %constr = @(xel)polyval(UrevUI,xel(1)) + polyval(rUI,xel(1))*xel(2) + polyval(sUI,xel(1))*log(polyval(tUI,xel(1))*xel(2) + 1) - V_electro(i-1);   
    
    %[xel, obj, info, iter, nf, laambda] = sqp(xel0, phi, constr, [],[],[],500,10e-5);
    %ET_p(j) = xel(1);
    %EI_elecp(j) = xel(2);
   end
   ET(i) = mean(ET_p);
   UI = @(x)polyval(UrevUI,ET(i)) + polyval(rUI,ET(i))*x + polyval(sUI,ET(i))*log(polyval(tUI,ET(i))*x + 1) - V_electro(i);
   EI_elec(i) = (fsolve(UI, 45,options))*2.5;
   
   P = zeros(2,2);
   for j=1:Np
    e = [ET_p(j);EI_elecp(j)] - [ET(i);EI_elec(i)];
    P = P + e*e';
   end
   P = P./(Np-1);
  
   EP_elec(i) = EI_elec(i)*V_electro(i)/1000;
   TP_elec(i) = 21*EP_elec(i);
 end


             %%-------FUEL CELL--------%%

 y_cap = [MHyd(i); MU_FC(i)];

 for j=1:20
  x= y_cap;
  A_FC = [1 , 0; 0, -1];
  f = [x(1) - eta_FC*N_FC*I_FC(i)*0.09/(2000*C_H2); U_FCo + e1*T_FC + e2*log(I_FC(i)/I_FCo) + R_FC*I_FC(i)/T_FC - x(2)];
  y = y_cap - Var_FC*(A_FC')*(inv(A_FC*Var_FC*A_FC'))*f;
  y_cap = y;
 end
 EU_FC(i) = y(2);
 EHyd(i) = y(1);
 EI_FC(i) = 2000*C_H2*EHyd(i)/(N_FC*0.09*eta_FC);
 EP_FC(i) = EI_FC(i)*EU_FC(i)/1000; 

end

%% -------------MSE calculations ----------------%%
RMSE_PEsolar = (immse(P_solar,EP_solar'))^0.5
RMSE_PMsolar = (immse(P_solar,MP_solar))^0.5
RMSE_IMsolar = (immse(I_solar,MI_solar))^0.5
RMSE_IEsolar = (immse(I_solar,EI_solar'))^0.5
RMSE_VMsolar = (immse(V_solar,MV_solar))^0.5
RMSE_VEsolar = (immse(V_solar,EV_solar'))^0.5

RMSE_TM = (immse(T,MT))^0.5
RMSE_TE = (immse(T,ET))^0.5
RMSE_IMelec = (immse(I_electro,MI_electro))^0.5
RMSE_IEelec = (immse(I_electro,EI_elec))^0.5
RMSE_PEelec = (immse(P_electro,EP_elec'))^0.5

RMSE_PMFC = (immse(P_FC, MP_FC))^0.5
RMSE_PEFC = (immse(P_FC,EP_FC'))^0.5
RMSE_UMFC = (immse(U_FC,MU_FC))^0.5
RMSE_UEFC = (immse(U_FC,EU_FC'))^0.5
RMSE_VMHyd = (immse(Hyd,MHyd))^0.5
RMSE_VEHyd = (immse(Hyd,EHyd'))^0.5

elapsed_time = toc()
Computation_time = cputime()

figure(1)
plot(t,P_solar,t,MP_solar,'r*',t,EP_solar,'-g')
legend('True','Measured','Estimated')
xlabel('Time(hour)')
ylabel('Electricity(W)')
title('Maximum power output from Solar panel')

figure(2)
plot(t,I_solar,t,MI_solar,'r*',t,EI_solar,'g-');
legend('True','Measured','Estimated')
title('Solar')
xlabel('Time(hour)');
ylabel('Voltage/cell')

figure(3)
plot(V_solar,I_solar,'-b',MV_solar,MI_solar,'r*',EV_solar,EI_solar,'g*');
legend('True','Measured','Estimated')
title('solar')
ylabel('Current');
xlabel('Voltage/cell')

figure(4)
plot(t,T,'-b',t,MT,'r*',t,ET,'-g')
xlabel('time(hrs)')
ylabel('temp.(C)')
legend('true','measured','estimated')
title('Electrolyzer')

figure(5)
plot(t,I_electro,'-b',t,MI_electro,'r*',t,EI_elec,'-g')
xlabel('time(hrs)')
ylabel('Current(A)')
legend('true','measured','estimated')
title('Electrolyzer')

figure(6)
plot(t,P_electro,'-b',t,MP_electro,'r*',t,EP_elec,'-g')
xlabel('time(hrs)')
ylabel('Power(KW)')
legend('true','measured','estimated')
title('Electrolyzer/Cell')

figure(8)
plot(t,P_FC,'-b',t,MP_FC,'r*',t,EP_FC,'-g')
xlabel('time(hrs)')
ylabel('Power(kW)')
title('Fuel Cell')

figure(9)
plot(t,U_FC,'-b',t,MU_FC,'r*',t,EU_FC,'-g')
legend('true','Measured','Estimated')
xlabel('time(hrs)')
ylabel('Voltage')

figure(10)
plot(t,I_FC,'-b',t,EI_FC,'-g')
legend('true','Estimated')
xlabel('time(hrs)')
ylabel('Current')

figure(11)
plot(t,Hyd,'-b',t,MHyd,'r*',t,EHyd,'-g')
title('Fuel cell')
legend('true','measured','estimated')
xlabel('time(hrs)')
ylabel('Hydrogen flow rate(kmol/hr)')