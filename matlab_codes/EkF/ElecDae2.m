function[dydt]=ElecDae2(t,y,V_elec)
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

dydt = zeros(2,1);
dydt(1) = (1/Ct)*(3.6*eta_c*(V_elec-U_tn)*y(2) - 3.6*(1/Rt)*(y(1)-Ta) - 1000000*F*Ccw*(-Tcw + y(1))*(1-exp(-0.00000036*( 9.557 + 0.0091935*y(2))/(F*Ccw))));
dydt(2) = 1.249-0.0008182*y(1) + polyval(rUI,y(1))*y(2)/2.5 + (polyval(sUI,y(1)))*log(polyval(tUI,y(1))*y(2)/2.5 + 1) - V_elec;
end