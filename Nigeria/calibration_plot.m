% First run 'Show_model_fts.m' and then run this code to highlight only the
% required plot for Global Plan report 

% --- Comparing model outputs with simulations
figure;
subplot(1,2,1); hold on;
% Plot data
plt  = [data.inc2000; data.inc2022; data.inc_h1; data.noti2022; data.noti_cu/10]'; 
hilo = diff(plt,1); md = plt(2,:);
xpts = [1:length(plt)]-0.1;

plot(xpts, md, 'r.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'Color', 'r');

% Plot simulations
plt = prctile(sims(:,1:5),[2.5,50,97.5],1);
hilo = diff(plt,1); md = plt(2,:);
xpts = [1:size(hilo,2)]+0.1;

plot(xpts, md, 'b.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'Color', 'b');
% xlim([0.5 5.5]);
ylabel('Rate per 100,000 population');
set(gca,'fontsize',fs,'XTick',1:size(plt,2),'XTickLabel',{'Incidence 2000','Incidence 2022', 'Incidence(HIV TB)2022', 'Notification 2022','Cumulative notification/10'});
xtickangle(45); 



subplot(1,2,2); hold on; 
% Plot data
plt  = [data.mort_H02000;data.mort_H02022; data.mort_H12022; data.ART_covg; data.HIV_prev; data.mdr2022; data.mdriniTX; 100*data.sym]'; 
hilo = diff(plt,1); md = plt(2,:);
xpts = [1:length(plt)]-0.1;

pl1=plot(xpts, md, 'r.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'Color', 'r');

% Plot simulations
plt = prctile(sims,[2.5,50,97.5],1);
hilo = diff(plt(:,6:end),1); md = plt(2,6:end);
xpts = [1:size(hilo,2)]+0.1;

pl2=plot(xpts, md, 'b.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'Color', 'b');
% xlim([0.5 3.5]);
legend([pl1(1),pl2(1)],'Data','Model')
set(gca,'fontsize',fs,'XTick',1:size(plt,2),'XTickLabel',{'Mortality(HIV-ve)2000','Mortality(HIV-ve)2022', 'Mortality(HIV+ve)2022', 'ART coverage', 'HIV prevalence','MDR 2022','MDR Tx-initiation','% Symptomatic'});
xtickangle(45);

sgtitle('Nigeria')
