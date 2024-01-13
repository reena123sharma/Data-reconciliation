I_electro = xlsread("data_central","Sheet3","F3:F1143");
t = xlsread("data_central","Sheet3","A3:A1143");
%%----------Fuel Cell
U_FCo = 33.18; %V
e1 = -0.013; %V.C^*1
e2 = -1.57; %V.C^-1
I_FCo = 8.798; %A
R_FC = -2.04; %ohm.C^-1
N_FC = 35;
C_H2 = 2.39; %Ah/l_conversion coefficient for hydrogen
eta_FC = 0.7; %utilization factor
T_FC = 25; % as we are not considering any temperature dynamics, we assume it to be constant
%--assuming constant Hydrogen supply
%Hyd_supp = 0.1*ones(m,1); % kmol/hr
%assuming no storage whatever electrolyzer produce, fuel cell consumes
Hyd_gen = 3.6*0.98*21*I_electro./(2*98485);%kmol/hr
% figure(4)
% plot(t,Hyd_gen)
% xlabel('time(hrs)')
% ylabel('Hydrogen flow rate(kmol/hr)')

I_FC = 2000*C_H2*Hyd_gen./(N_FC*0.09*eta_FC);
U_FC = U_FCo + e1*T_FC + e2*log(I_FC./I_FCo) + R_FC*I_FC./T_FC;
P_FC = U_FC.*I_FC/1000;

%%----GENERATING MEASURED DATA-------------%%
m = length(t);
MHyd = Hyd_gen +  0.01*randn(m,1);
MU_FC = U_FC + randn(m,1); 
MI_FC = I_FC + 2.0*randn(m,1);
MP_FC = MI_FC.*MU_FC/1000;  %%(kW) 

Var_FC = [0.5,0,0;0 ,0.5,0;0,0,0.01];

for i=1:m
 %%-------FUEL CELL--------%%
 y = [MI_FC(i); MU_FC(i); MHyd(i)];
 y_cap = y;
 for j=1:100
  I = y_cap(1);
  U = y_cap(2);
  H = y_cap(3);
  A_FC = [e2/I + R_FC/T_FC , -1, 0;  eta_FC*N_FC*0.09/(2000*C_H2), 0, -1];
  f = [U_FCo + e1*T_FC + e2*log(I/I_FCo) + R_FC*I/T_FC - U; eta_FC*N_FC*I*0.09/(2000*C_H2) - H];
  b = A_FC*y_cap - f;
  y_cap = y - Var_FC*(A_FC')*(inv(A_FC*Var_FC*A_FC'))*(A_FC*y - b);
 end
 EI_FC(i) = y_cap(1);
 EU_FC(i) = y_cap(2);
 EHyd(i) = y_cap(3);
 %EI_FC(i) = 2000*C_H2*EHyd(i)/(N_FC*0.09*eta_FC);
 EP_FC(i) = EI_FC(i)*EU_FC(i)/1000; 

end

RMSE_PMFC = (immse(P_FC, MP_FC))^0.5
RMSE_PEFC = (immse(P_FC,EP_FC'))^0.5
RMSE_UMFC = (immse(U_FC,MU_FC))^0.5
RMSE_UEFC = (immse(U_FC,EU_FC'))^0.5
RMSE_IMFC = (immse(I_FC,MI_FC))^0.5
RMSE_IEFC = (immse(I_FC,EI_FC'))^0.5
RMSE_VMHyd = (immse(Hyd_gen,MHyd))^0.5
RMSE_VEHyd = (immse(Hyd_gen,EHyd'))^0.5