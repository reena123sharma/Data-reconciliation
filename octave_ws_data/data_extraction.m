
t = xlsread('data_central','Sheet3','A3:A1143');
G = xlsread('data_central','Sheet3','B3:B1143');
V_solar = xlsread('data_central','Sheet3','C3:C1143');
I_solar = xlsread('data_central','Sheet3','D3:D1143');
P_solar = V_solar.*I_solar;

V_electro = xlsread('data_central','Sheet3','E3:E1143');
I_electro = xlsread('data_central','Sheet3','F3:F1143');
T = xlsread('data_central','Sheet3','G3:G1143');
P_electro = V_electro.*I_electro/1000;

m = length(t);

U_FC = xlsread('data_central','Sheet3','H3:H1143');
I_FC = xlsread('data_central','Sheet3','I3:I1143');
Hyd = 0.1*ones(m,1); % kmol/hr
P_FC = U_FC.*I_FC/1000;

