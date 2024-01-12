m2 = 4.0; %in hrs

%%%%  loading temperature failure data  %%%%%%%
load('shutdownT_RPF_2s.mat','ET','EI_elec','m','MT', 'MI_electro','m_start','m_end','Ts');
ET_RPFt = ET(1:m2/Ts);
EI_RPFt = EI_elec(1:m2/Ts);
MT_elect = MT(1:m2/Ts);
MI_elect = MI_electro(1:m2/Ts);
clear ET;
clear EI_elec;
clear MT;
clear MI_electro;

load('shutT_ENKF_2_res.mat','ET','EI_elec');
ET_EnKFt = ET(1:m2/Ts);
EI_EnKFt = EI_elec(1:m2/Ts);
clear ET;
clear EI_elec;

load('shutT_EKF_2_res.mat','ET','EI_elec')
ET_EKFt = ET(1:m2/Ts);
EI_EKFt = EI_elec(1:m2/Ts);
clear ET;
clear EI_elec;

%%%%  loading Current failure data  %%%%%%%
load('shutdownI_RPF_2s.mat','ET','EI_elec','MT', 'MI_electro');
ET_RPFi = ET(1:m2/Ts);
EI_RPFi = EI_elec(1:m2/Ts);
MT_eleci = MT(1:m2/Ts);
MI_eleci = MI_electro(1:m2/Ts);
clear ET;
clear EI_elec;
clear MT;
clear MI_electro;

load('shutI_ENKF_2_res.mat','ET','EI_elec');
ET_EnKFi = ET(1:m2/Ts);
EI_EnKFi = EI_elec(1:m2/Ts);
clear ET;
clear EI_elec;

load('shutI_EKF_2_res.mat','ET','EI_elec')
ET_EKFi = ET(1:m2/Ts);
EI_EKFi = EI_elec(1:m2/Ts);
clear ET;
clear EI_elec;

load('Ts_2.mat','T','I_electro','t')
t_e = t(1:m2/Ts);
T_elec = T(1:m2/Ts);
I_elec = I_electro(1:m2/Ts);
clear T; clear I_electro;

for i=1:m
    if (i>=m_start) && (i<=m_end)
        MT_elect(i) = nan;
        MI_eleci(i) = nan;
    end
end

%%%%%%%%%%%%%%            plots        %%%%%%%%%%%%%%%%%%%%
tbars = [1 2];
fig = figure('units','inch','position',[0,0,14,14]);

subplot(2,2,1)
plot(t_e, MT_elect, ':go','MarkerIndices',1:10:length(t_e),'MarkerSize',3,'MarkerFaceColor','g','LineWidth',0.01)
hold on;
plot(t_e,ET_EKFt,':o','Color',[0 0.4470 0.7410],'MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3,'MarkerFaceColor','b')
hold on;
plot(t_e, ET_EnKFt,':kd','MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3.5,'MarkerFaceColor','k')
hold on;
plot(t_e, ET_RPFt,'-r','LineWidth',0.9)
hold on;
plot(t_e, T_elec, '-y','LineWidth',0.7)
hold on;
p1 = patch( [tbars(1) tbars(1), tbars(2) tbars(2)],[45 65 65 45], [0.8 0.8 0.8]);
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel('Temperature(\circ C)','FontSize',12,'FontWeight','bold')
title('(a) Temperature sensor failure')
set(gca, 'FontSize',10,'LineWidth',0.8);
alpha(0.25)
p1.EdgeColor = 'none';

subplot(2,2,3);
plot(t_e, MI_elect, ':go','MarkerIndices',1:10:length(t_e),'MarkerSize',3,'MarkerFaceColor','g','LineWidth',0.01)
hold on;
plot(t_e,EI_EKFt,':o','Color',[0 0.4470 0.7410],'MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3,'MarkerFaceColor','b')
hold on;
plot(t_e, EI_EnKFt,':kd','MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3.5,'MarkerFaceColor','k')
hold on;
plot(t_e, EI_RPFt,'-r','LineWidth',0.9)
hold on;
plot(t_e, I_elec, '-y','LineWidth',0.7)
hold off;
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel('Current(Amp.)','FontSize',12,'FontWeight','bold')
[h,icons] = legend({'Measured','DAE-EKF','DAE-EnKF','DAE-RPF','True'},'FontSize',10,'Location','southeast','NumColumns',3);
set(gca, 'FontSize',10,'LineWidth',0.8);
ylim([70 165])

icons = findobj(icons,'Type','line');
set(icons,'LineWidth',1.3);
icons = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize',5);

%%%%%%%%%%%%%%  Current failure   %%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)  
plot(t_e, MT_eleci, ':go','MarkerIndices',1:10:length(t_e),'MarkerSize',3,'MarkerFaceColor','g','LineWidth',0.01)
hold on;
plot(t_e,ET_EKFi,':o','Color',[0 0.4470 0.7410],'MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3,'MarkerFaceColor','b')
hold on;
plot(t_e, ET_EnKFi,':kd','MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3.5,'MarkerFaceColor','k')
hold on;
plot(t_e, ET_RPFi,'-r','LineWidth',0.9)
hold on;
plot(t_e, T_elec, '-y','LineWidth',0.7)
hold off;
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel('Temperature(\circ C)','FontSize',12,'FontWeight','bold')
title('(b) Current sensor failure')
set(gca, 'FontSize',10,'LineWidth',0.8);

subplot(2,2,4);
plot(t_e, MI_eleci, ':go','MarkerIndices',1:10:length(t_e),'MarkerSize',3,'MarkerFaceColor','g','LineWidth',0.01)
hold on;
plot(t_e,EI_EKFi,':o','Color',[0 0.4470 0.7410],'MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3,'MarkerFaceColor','b')
hold on;
plot(t_e, EI_EnKFi,':kd','MarkerIndices',1:10:length(t_e),'LineWidth',0.1,'MarkerSize',3.5,'MarkerFaceColor','k')
hold on;
plot(t_e, EI_RPFi,'-r','LineWidth',0.9)
hold on;
plot(t_e, I_elec, '-y','LineWidth',0.7)
hold on;
p1 = patch( [tbars(1) tbars(1), tbars(2) tbars(2)],[60 170 170 60], [0.8 0.8 0.8]);
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel('Current(Amp.)','FontSize',12,'FontWeight','bold')
set(gca, 'FontSize',10,'LineWidth',0.8);
ylim([70 160])
alpha(0.25)
p1.EdgeColor = 'none';
