Ts = 0.5/3600; %sampling time in seconds

%%-----------------------SOLAR PANEL----------------%%
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

Tim = xlsread('irridanceVstime','Sheet1','C2:C171');
S = xlsread('irridanceVstime','Sheet1','B2:B171');
G_prof = splinefit(Tim,S,40);
Time = 0:Ts:9.5;

##fid = fopen("ascasc.txt","w");
##fdisp(fid,T);
##fclose(fid);

G = ppval(G_prof,Time);
maxIter = length(Time);

figure(11)
plot(Tim,S,'-r',Time,G,'-g')


Iph = G.*Iph_ref/G_ref; %%Assume the cell temp. = ref cell temp. = 25C
Io = Io_ref;
a = a_ref;
Rs = Rs_ref;%%assumption
Rp = Rp_ref*G_ref./G;

%%---solve for maximum power for different G
for i=1:1:maxIter
P_prev =0;
 for V=0:0.5:25
 y = @(I)Iph(i) - Io*(exp((V + Rs*I)/a) - 1) - (V + I*Rs)/Rp(i) - I;
 curr = fzero(y,2.);
 P = V*curr;
  if (P<P_prev)
  break;
  endif
 P_prev = P;
 endfor  
P_solar(i) = P;
I_solar(i) = curr;
V_solar(i) =V;
end

figure(1)
plot(Time,P_solar)
xlabel('Time(hour)')
ylabel('Electricity(W)')
title('Maximum power output from Solar panel')

figure(2)
plot(Time,I_solar);
title('Voltage/cell Vs time')
xlabel('Time(hour)');
ylabel('Voltage/cell')

figure(3)
plot(V_solar,I_solar,'b*');
title('Voltage/cell Vs time')
ylabel('Current');
xlabel('Voltage/cell')

TV_solar = 2*V_solar;%%Assume 2 PWX 500PV in series

%%-------------Hydrogen generation Electrolyzer(PHOEBUS)-----------
V_electro = TV_solar/21; %% assume that we are getting correct values of voltage
m = length(V_electro);
%Hydr = 3600*0.98*21*Current./(2*98485);%mol/hr

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

m = maxIter;
T = zeros(m,1);
T(1) = T_in;
for i=2:m
 func = @(I) 1.249-0.0008182*T(i-1) + polyval(rUI,T(i-1))*(I) + (polyval(sUI,T(i-1)))*log(polyval(tUI,T(i-1))*(I) + 1) - V_electro(i);
 I_electro(i) = (fzero(func,40))*2.5;
 
 UA_t = 9.557 + 0.0091935*I_electro(i);
 fthe = @(x) T(i-1) + (Ts/Ct)*(3.6*eta_c*(V_electro(i)-U_tn)*I_electro(i) - 3.6*(1/Rt)*(x-Ta) - 1000000*F*Ccw*(-Tcw + x)*(1-exp(-0.00000036*UA_t/(F*Ccw))) ) - x;
 T(i) = fzero(fthe,50); 
end
I_electro(1) = I_electro(2);
t = Time;
P_electro = V_electro(1:m).*I_electro/1000;

figure(4)
plot(t,T,'-b')
xlabel('time(hrs)')
ylabel('temp.(C)')
title('Electrolyzer')

figure(5)
plot(t,I_electro,'-b')
xlabel('time(hrs)')
ylabel('Current(A)')
title('Electrolyzer')

figure(6)
plot(t,P_electro,'-b')
xlabel('time(hrs)')
ylabel('Power(KW)')
title('Electrolyzer/Cell')

%%----------Fuel Cel

U_FCo = 33.18; %V
e1 = -0.013; %V.C^*1
e2 = -1.57; %V.C^-1
I_FCo = 8.798; %A
R_FC = -2.04; %ohm.C^-1
N_FC = 35;
C_H2 = 2.39; %Ah/l_conversion coefficient for hydrogen
eta_FC = 0.7; %utilization factor
T_FC = 298;

%--assuming constant Hydrogen supply
Hyd_supp = 0.1*ones(m,1); % kmol/hr
I_FC = 2000*C_H2*Hyd_supp./(N_FC*0.09*eta_FC);
U_FC = U_FCo + e1*T_FC + e2*log(I_FC./I_FCo) + R_FC*I_FC./T_FC;
P_FC = U_FC.*I_FC/1000;

figure(8)
plot(t,P_FC,'-b')
xlabel('time(hrs)')
ylabel('Power(kW)')
title('Fuel Cell')

figure(9)
plot(t,U_FC,'-b')
xlabel('time(hrs)')
ylabel('Voltage')

figure(10)
plot(t,I_FC,'-b')
xlabel('time(hrs)')
ylabel('Current')
