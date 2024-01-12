load('Np50_Ts6s_RPF','ET','EI_elec','T','I_electro','t','m');
ET_RPF = ET;
EI_RPF = EI_elec;
t_6 = t;
m_6 = m;
clear ET;
clear EI_elec;
clear t;
clear m;
load('Np20_ts6s_EnKF.mat','ET','EI_elec');
ET_EnKF = ET;
EI_EnKF = EI_elec;
clear ET;
clear EI_elec;
load('Ts_1s_EKF.mat','ET','EI_elec','Ts')
ET_EKF = ET;
EI_EKF = EI_elec;
t_1 = (0:Ts:9.5)';
m_1 = length(t_1);
clear ET;
clear EI_elec;
clear Ts;
j = 1;
for i =1:10:m_6
    t_6_m(j) = t_6(i);
    ET_RPF_m(j) = ET_RPF(i);
    ET_EnKF_m(j) = ET_EnKF(i);
    T_m(j) = T(i);
    I_m(j) = I_electro(i);
    EI_RPF_m(j) = EI_RPF(i);
    EI_EnKF_m(j) = EI_EnKF(i);
    j = j+1;
end 
j=1;
for i=1:30:m_1
    t_1_m(j) = t_1(i);
    ET_EKF_m(j) = ET_EKF(i);
    EI_EKF_m(j) = EI_EKF(i);
    j = j+1;
end

figure(1)
p = plot( t_1_m,ET_EKF_m,'-y', t_6_m, ET_EnKF_m,'-c',t_6_m, ET_RPF_m,'-r',t_6_m, T_m, '-k')
xlabel('Time(hrs)')
ylabel('Temp.(C)')
legend('DAE-EKF','DAE-EnKF','RPF','True')
title('Electrolyzer')
p(4).LineWidth = 1.5;
p(1).LineWidth = 1.5;
p(2).LineWidth = 1.5;
p(3).LineWidth = 2;
set(gca,'FontSize',20);

figure(2)
p2 = plot( t_1_m,EI_EKF_m,'-y', t_6_m, EI_EnKF_m,'-c',t_6_m, EI_RPF_m,'-r',t_6_m, I_m, '-k')
xlabel('Time(hrs)')
ylabel('Current(A)')
legend('DAE-EKF','DAE-EnKF','RPF','True')
title('Electrolyzer')
title('Electrolyzer')
p2(4).LineWidth = 1.5;
p2(1).LineWidth = 1.5;
p2(2).LineWidth = 1.5;
p2(3).LineWidth = 2;
set(gca,'FontSize',20);
