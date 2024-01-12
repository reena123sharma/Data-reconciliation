load('ts_6.mat','T','I_electro','t')

load('Np50_Ts6s_RPF.mat','ET','EI_elec','m','MT', 'MI_electro');
ET_RPF = ET;
EI_RPF = EI_elec;
clear ET;
clear EI_elec;

load('Single_EnKF_20_6_res1.mat','ET','EI_elec');
ET_EnKF = ET;
EI_EnKF = EI_elec;
clear ET;
clear EI_elec;

load('SingleEKF_6s_res.mat','ET','EI_elec','Ts')
ET_EKF = ET;
EI_EKF = EI_elec;
clear ET;
clear EI_elec;

fig = figure('units','inch','position',[0,0,8,14]);
ax1 = subplot(2,1,1);
plot(t, MT, ':go','MarkerIndices',1:10:length(t),'MarkerSize',3,'MarkerFaceColor','g','LineWidth',0.01)
hold on;
plot(t,ET_EKF,':o','Color',[0 0.4470 0.7410],'MarkerIndices',1:10:length(t),'LineWidth',0.1,'MarkerSize',3,'MarkerFaceColor','b')
hold on;
plot(t, ET_EnKF,':kd','MarkerIndices',1:10:length(t),'LineWidth',0.1,'MarkerSize',3.5,'MarkerFaceColor','k')
hold on;
plot(t, ET_RPF,'-r','LineWidth',0.9)
hold on;
plot(t, T, '-y','LineWidth',0.7)
hold off;
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel('Temperature(\circ C)','FontSize',12,'FontWeight','bold')
[h,icons] = legend('Measured','DAE-EKF','DAE-EnKF','DAE-RPF','True','FontSize',10,'Location','southeast')
title('(a) Estimated values of differential state x')
set(gca, 'FontSize',10,'LineWidth',0.8);
ylim([45 72])

icons = findobj(icons,'Type','line');
set(icons,'LineWidth',1.3);
icons = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize',5);

ax2 = subplot(2,1,2)
plot(t, MI_electro, ':go','MarkerIndices',1:10:length(t),'MarkerSize',3,'MarkerFaceColor','g','LineWidth',0.01)
hold on;
plot(t,EI_EKF,':o','Color',[0 0.4470 0.7410],'MarkerIndices',1:10:length(t),'LineWidth',0.1,'MarkerSize',3,'MarkerFaceColor','b')
hold on;
plot(t, EI_EnKF,':kd','MarkerIndices',1:10:length(t),'LineWidth',0.1,'MarkerSize',3.5,'MarkerFaceColor','k')
hold on;
plot(t, EI_RPF,'-r','LineWidth',0.9)
hold on;
plot(t, I_electro, '-y','LineWidth',0.7)
hold off;
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel('Current(Amp.)','FontSize',12,'FontWeight','bold')
%legend('DAE-EKF','DAE-EnKF','DAE-RPF','True','FontSize',10,'Location','southeast')
title('(b) Estimated values of algebraic state z')
set(gca, 'FontSize',10,'LineWidth',0.8);
ylim([70 190])

