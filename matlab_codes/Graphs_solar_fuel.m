load('solar_6s_res.mat','t','MI_solar','EI_solar','I_solar')

%%%%%%%%%%%%%%%%%%%%%%%%Solar Cell%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fig = figure('units','inch','position',[0,0,12,5])
%subplot(1,2,1)
plot(t, MI_solar,':o','Color',[0.8500 0.3250 0.0980],'MarkerIndices',1:5:length(t),'MarkerSize',3,'LineWidth',0.8,'MarkerFaceColor',[0.8500 0.3250 0.0980])
hold on
plot(t, EI_solar,'-d','Color',[0 0.4470 0.7410],'MarkerIndices',1:5:length(EI_solar),'LineWidth',0.5,'MarkerSize',3,'MarkerFaceColor','b')
hold on
plot(t, I_solar,'-y','LineWidth',0.7)
hold off
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel( 'Current/Cell(A)','FontSize',12,'FontWeight','bold')
[h,icons] = legend('Measured', 'Estimated','True','FontSize',12)
set(gca, 'FontSize',10,'LineWidth',0.8)
ylim([0,3.5])
%to change the sizes in legend only 
icons = findobj(icons,'Type','line');
set(icons,'LineWidth',1.3);
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',5);
% subplot(1,2,2)
% plot(t, MV_solar,':*','Color',[0.8500 0.3250 0.0980],'MarkerIndices',1:100:length(MV_solar))
% hold on
% plot(t, EV_solar,':d','Color',[0 0.4470 0.7410],'MarkerIndices',1:100:length(EV_solar),'LineWidth',0.2,'MarkerSize',7 )
% hold off
% xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
% ylabel( 'Voltage/Cell(V)','FontSize',12,'FontWeight','bold')
% legend('Measured', 'Estimated','FontSize',12)
% set(gca, 'FontSize',12,'LineWidth',1.2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fuel Cell%%%%%%%%%%%%%%%%%%%%%%%%
load('fuel_6s_res.mat','t','U_FC','MU_FC','EU_FC','Hyd_gen','MHyd','EHyd')

fig = figure('units','inch','position',[0,0,7,12])
subplot(2,1,1)
plot(t, MU_FC,':o','Color',[0.8500 0.3250 0.0980],'MarkerIndices',1:5:length(t),'MarkerSize',3,'LineWidth',0.8,'MarkerFaceColor',[0.8500 0.3250 0.0980])
hold on
plot(t, EU_FC,'-d','Color',[0 0.4470 0.7410],'MarkerIndices',1:5:length(t),'LineWidth',0.5,'MarkerSize',3,'MarkerFaceColor','b')
hold on
plot(t, U_FC,'-y','LineWidth',0.9)
hold off
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel( 'Voltage(V)','FontSize',12,'FontWeight','bold')
%[h,icons] = legend('Measured', 'Estimated','True','FontSize',12)
title('(a) Estimated values of Voltage(V)')
set(gca, 'FontSize',10,'LineWidth',0.8)

subplot(2,1,2)
plot(t, MHyd,':o','Color',[0.8500 0.3250 0.0980],'MarkerIndices',1:5:length(t),'MarkerSize',3,'LineWidth',0.8,'MarkerFaceColor',[0.8500 0.3250 0.0980])
hold on
plot(t, EHyd,'-d','Color',[0 0.4470 0.7410],'MarkerIndices',1:5:length(t),'LineWidth',0.5,'MarkerSize',3,'MarkerFaceColor','b')
hold on
plot(t, Hyd_gen,'-y','LineWidth',0.9)
hold off
xlabel('Time(hrs)','FontSize',12,'FontWeight','bold')
ylabel( 'Hydrogen flow rate(kmol/hr)','FontSize',12,'FontWeight','bold')
title('(b) Estimated values of Hydrogen flow rate')
[h,icons] = legend('Measured', 'Estimated','True','FontSize',12,'Location','southeast')
set(gca, 'FontSize',10,'LineWidth',0.8)

icons = findobj(icons,'Type','line');
set(icons,'LineWidth',1.3);
icons = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize',5);

