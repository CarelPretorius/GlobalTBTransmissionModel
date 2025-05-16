clear all; load calibres2.mat  % 

tic
obj = @(x) get_livessaved(x, prm, ref, sel, agg, gps,gps1,gps2,lhd);

load Model_setup;
xsto =xsto3; 

ix0 = 3e4; nx = 15; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);

pct_xs= prctile(xs,[2.5,50,97.5])

% inds = find(outsto==max(outsto));
% xs = xsto(inds(1),:);

pop2003_2022 = data.pop00_22(4:end);
pop2000_2022 = data.pop00_22;

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sims(ii,:) = [aux.inc2000,aux.inc2022, aux.inc_h1, aux.noti, aux.noti_cu/10, aux.mort_H02000,aux.mort_H02022, aux.mort_H1, aux.ART_covg, aux.HIV_prev, aux.mdr2022, aux.mdriniTX,100*aux.sym];
    
    pop = sum(aux.soln(1:end-1,1:i.nstates),2);  %1997:2022 (1:end -> 1997:2023)
    inct(:,:,ii) = (diff(aux.soln(:,i.aux.inc),1)./pop)*1e5;  %1997:2022

    % Lives saved calculation (below are already normalized with the population)
    incd(:,ii)  = aux.incd;   % Current
    incdb(:,ii) = aux.incdb;  % Null coverage from 2003 
    incdc(:,ii) = aux.incdc;  % Constant coverage from 2003 
    incdd(:,ii) = aux.incdd;  % Null coverage from 2000 
    incde(:,ii) = aux.incde;  % Constant coverage from 2000 
    incdf(:,ii) = aux.incdf;  % Null coverage from 2022 
    incdg(:,ii) = aux.incdg;  % Constant coverage from 2022 

    mor(:,ii)  = aux.mort;   % Current 1997 - 2022
    morb(:,ii) = aux.mortb;  % Null scenario from 2000
    morc(:,ii) = aux.mortc;  % Constant coverage at 2000
    mord(:,ii) = aux.mortd;  % Null coverage at 2003
    more(:,ii) = aux.morte;  % Constant coverage at 2003
    morf(:,ii) = aux.mortf;  % Null coverage at 2022
    morg(:,ii) = aux.mortg;  % Constant coverage at 2022
     
    mor0(:,ii)= ((aux.mort(4:end))./1e5).*pop2000_2022'; 
    mor1(:,ii)= ((aux.mort(7:end))./1e5).*pop2003_2022'; 
    morb1(:,ii) = ((aux.mortb)./1e5).*pop2000_2022';
    morc1(:,ii) = ((aux.mortc)./1e5).*pop2000_2022';
    mord1(:,ii) = ((aux.mortd)./1e5).*pop2003_2022';
    more1(:,ii) = ((aux.morte)./1e5).*pop2003_2022';
    
    % Lives saved by null coverge (LS1) and by cont coverage (LS2) from 2000
    LS1(:,ii) = (morb1(:,ii) - mor0(:,ii)); % 2000-2022
    LS2(:,ii) = (morc1(:,ii) - mor0(:,ii)); % 2000-2022
    
    % Lives saved by null coverge (LS3) and by cont coverage (LS4) from 2003
    LS3(:,ii) = (mord1(:,ii) - mor1(:,ii)); % 2003-2022
    LS4(:,ii) = (more1(:,ii) - mor1(:,ii)); % 2003-2022
    
end
fprintf('\n');

inc_pct = permute(prctile(inct,[2.5,50,97.5],3),[3,1,2]);                  % 1. Lo/Md/Hi, 2.Time, 3.All/HIV only
     

incd_pct = prctile(incd,[2.5,50,97.5],2);    %1997-2019 % Current
incdb_pct = prctile(incdb,[2.5,50,97.5],2);  %2003-2019 % Null coverage
incdc_pct = prctile(incdc,[2.5,50,97.5],2);  %2003-2019 % Constant Coverage
incdd_pct = prctile(incdd,[2.5,50,97.5],2);  %2000-2019 % Null coverage
incde_pct = prctile(incde,[2.5,50,97.5],2);  %2000-2019 % Constant Coverage
incdf_pct = prctile(incdf,[2.5,50,97.5],2);  %2022-2035 % Null Coverage
incdg_pct = prctile(incdg,[2.5,50,97.5],2);  %2022-2035 % Constant Coverage

mor_pct = prctile(mor,[2.5,50,97.5],2);    %1997-2019 % Current
morb_pct = prctile(morb,[2.5,50,97.5],2);  %2003-2019 % Null coverage
morc_pct = prctile(morc,[2.5,50,97.5],2);  %2003-2019 % Constant Coverage
mord_pct = prctile(mord,[2.5,50,97.5],2);  %2000-2019 % Null coverage
more_pct = prctile(more,[2.5,50,97.5],2);  %2000-2019 % Constant Coverage
morf_pct = prctile(morf,[2.5,50,97.5],2);  %2022-2035 % Null Coverage
morg_pct = prctile(morg,[2.5,50,97.5],2);  %2022-2035 % Constant Coverage

LS1_pct = prctile(LS1,[2.5,50,97.5],2); % 2000-2022
LS2_pct = prctile(LS2,[2.5,50,97.5],2); % 2000-2022
LS3_pct = prctile(LS3,[2.5,50,97.5],2); % 2003-2022
LS4_pct = prctile(LS4,[2.5,50,97.5],2); % 2003-2022

% --- Show all on a plot --------------------------------------------------

figure; fs = 14;

% --- Incidence plot
% subplot(2,2,1);
% 
% plt = inc_pct(:,:,1);  % All TB  %sum(inc_pct,3); 
% plot(1997:2022,plt(2,:)); hold on;
% jbfill(1997:2022,plt(3,:),plt(1,:),'b','None',1,0.3); hold on; %1:size(plt,2)
% yl = ylim; yl(1) = 0; ylim(yl);
% %xlim([2011 2022]);
% 
% plt = inc_pct(:,:,2); % HIV-TB %sum(inc_pct(:,:,[2,3]),3);  
% plot(1997:2022,plt(2,:)); hold on;
% jbfill(1997:2022,plt(3,:),plt(1,:),'r','None',1,0.3); hold on;
% yl = ylim; yl(1) = 0; ylim(yl);
% %xlim([2011 2022]);
% 
% ylabel('Incidence');
% set(gca,'fontsize',fs);

% subplot(2,2,2)
% plt = inc_pct(:,:,3); % MDR-TB %sum(inc_pct(:,:,[2,3]),3);  
% plot(1997:2022,plt(2,:)); hold on;
% jbfill(1997:2022,plt(3,:),plt(1,:),'r','None',1,0.3); hold on;
% yl = ylim; yl(1) = 0; ylim(yl);
% %xlim([2011 2022]);
% 
% ylabel('MDR Incidence');
%set(gca,'fontsize',fs);

% --- Comparing model outputs with simulations
%subplot(2,2,3); hold on;
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
ylabel('Rate per 100,000 population','fontsize',fs);
set(gca,'fontsize',fs,'XTick',1:size(plt,2),'XTickLabel',{'Incidence 2000','Incidence 2022', 'Incidence(HIV TB)2022', 'Notification 2022','Cumulative notification/10'});
xtickangle(45); 

%subplot(2,2,4); hold on; 
subplot(1,2,2); hold on; 
% Plot data
plt  = [data.mort_H02000;data.mort_H02022; data.mort_H12022; data.ART_covg; data.HIV_prev; data.mdr2022; data.mdriniTX; 100*data.sym]'; 
hilo = diff(plt,1); md = plt(2,:);
xpts = [1:length(plt)]-0.1;

plot(xpts, md, 'r.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'Color', 'r');

% Plot simulations
plt = prctile(sims,[2.5,50,97.5],1);
hilo = diff(plt(:,6:end),1); md = plt(2,6:end);
xpts = [1:size(hilo,2)]+0.1;

plot(xpts, md, 'b.', 'markersize', 24);
errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'Color', 'b');
% xlim([0.5 3.5]);

set(gca,'fontsize',fs,'XTick',1:size(plt,2),'XTickLabel',{'Mortality(HIV-ve)2000','Mortality(HIV-ve)2022', 'Mortality(HIV+ve)2022', 'ART coverage', 'HIV prevalence','MDR 2022','MDR Tx-initiation','% Symptomatic'});
xtickangle(45);
sgtitle('Indonesia')








% Plot of mortality rate for all three scenarios
figure;  
plt = mor_pct(1:end,:)';  %current
lg(1,:) = plot(1997:2022, plt(2,:),'r','LineWidth',1.5);  %Mortality 1997:2019
jbfill(1997:2022, plt(3,:), plt(1,:), 'r', 'None', 1, 0.3); hold on

plt = morb_pct(1:end,:)'; % Null from 2000
lg(2,:) = plot(2000:2022, plt(2,:),'g','LineWidth',1.5);  %Mortality 2000:2022
jbfill(2000:2022, plt(3,:), plt(1,:), 'g', 'None', 1, 0.3); hold on

plt = morc_pct(1:end,:)'; % Constant coverage from 2000
lg(3,:) = plot(2000:2022, plt(2,:),'b','LineWidth',1.5);  %Mortality 2000:2022
jbfill(2000:2022, plt(3,:), plt(1,:), 'b', 'None', 1, 0.3); hold on

plt = mord_pct(1:end,:)'; % Null coverage from 2003
lg(4,:) = plot(2003:2022, plt(2,:),'y','LineWidth',1.5);  %Mortality 2005:2022
jbfill(2003:2022, plt(3,:), plt(1,:), 'y', 'None', 1, 0.3); hold on

plt = more_pct(1:end,:)'; % Constant coverage from 2005
lg(5,:) = plot(2003:2022, plt(2,:),'m','LineWidth',1.5);  %Mortality 2005:2022
jbfill(2003:2022, plt(3,:), plt(1,:), 'm', 'None', 1, 0.3); hold on

plt = morf_pct(1:end,:)'; % Null coverage from 2022
lg(6,:) = plot(2022:2035, plt(2,:),'c','LineWidth',1.5);  %Mortality 2022:2035
jbfill(2022:2035, plt(3,:), plt(1,:), 'c', 'None', 1, 0.3); hold on

plt = morg_pct(1:end,:)'; % Constant coverage from 2005
lg(7,:) = plot(2022:2035, plt(2,:),'m','LineWidth',1.5);  %Mortality 2022:2035
jbfill(2022:2035, plt(3,:), plt(1,:), 'm', 'None', 1, 0.3);

set(gca,'fontsize',fs);
ylabel('Rate per 100,000','fontsize',fs);
yl = ylim; yl(1) = 0; ylim(yl);
xlim([2000 2035]);
legend(lg(1:7,:),'Current trend','Null services from 2000','Constant coverage from 2000','Null services from 2003','Constant coverage from 2003','Null services from 2022','Constant coverage from 2022')
title('Mortality rate')


% Plot of incidence rate for all three scenarios
figure;
plt = incd_pct(1:end,:)';  %end-14
lg(1,:) = plot(1997:2022, plt(2,:),'r','LineWidth',1.5);  % incidence  2005:2019
jbfill(1997:2022, plt(3,:), plt(1,:), 'r', 'None', 1, 0.3); hold on

plt = incdb_pct(1:end,:)'; 
lg(2,:) = plot(2000:2022, plt(2,:),'g','LineWidth',1.5);  %incidence2000:2022
jbfill(2000:2022, plt(3,:), plt(1,:), 'g', 'None', 1, 0.3); hold on

plt = incdc_pct(1:end,:)'; 
lg(3,:) = plot(2000:2022, plt(2,:),'b','LineWidth',1.5);  %incidence2000:2022
jbfill(2000:2022, plt(3,:), plt(1,:), 'b', 'None', 1, 0.3); hold on

plt = incdd_pct(1:end,:)'; 
lg(4,:) = plot(2003:2022, plt(2,:),'y','LineWidth',1.5);  %incidence2005:2022
jbfill(2003:2022, plt(3,:), plt(1,:), 'y', 'None', 1, 0.3); hold on

plt = incde_pct(1:end,:)'; 
lg(5,:) = plot(2003:2022, plt(2,:),'m');  %incidence2005:2022
jbfill(2003:2022, plt(3,:), plt(1,:), 'm', 'None', 1, 0.3); hold on

plt = incdf_pct(1:end,:)'; % Null coverage from 2022
lg(6,:) = plot(2022:2035, plt(2,:),'c','LineWidth',1.5);  %Incidence 2022:2035
jbfill(2022:2035, plt(3,:), plt(1,:), 'c', 'None', 1, 0.3); hold on

plt = incdg_pct(1:end,:)'; % Constant coverage from 2005
lg(7,:) = plot(2022:2035, plt(2,:),'m','LineWidth',1.5);  %Incidence 2022:2035
jbfill(2022:2035, plt(3,:), plt(1,:), 'm', 'None', 1, 0.3);



set(gca,'fontsize',fs);
ylabel('Rate per 100,000','fontsize',fs);
yl = ylim; yl(1) = 0; ylim(yl);
xlim([2000 2035]);
legend(lg(1:7,:),'Current trend','Null services from 2000','Constant coverage from 2000','Null services from 2003','Constant coverage from 2003','Null services from 2022','Constant coverage from 2022')
title('Incidence rate')







figure; %Lives saved

plt = LS1_pct(1:end,:)'/1e3;
lg(1,:) = plot(2000:2022, plt(2,:),'y','LineWidth',1.5);  % LS from 2000:2022 (null covegare)
jbfill(2000:2022, plt(3,:), plt(1,:), 'y', 'None', 1, 0.3); hold on

plt = LS2_pct(1:end,:)'/1e3;
lg(2,:) = plot(2000:2022, plt(2,:),'m','LineWidth',1.5);  % LS from 2000:2022 (constant covegare)
jbfill(2000:2022, plt(3,:), plt(1,:), 'm', 'None', 1, 0.3);  hold on

plt = LS3_pct(1:end,:)'/1e3;  %end-14
lg(3,:) = plot(2003:2022, plt(2,:),'g','LineWidth',1.5);  % LS from 2005:2022 (null covegare)
jbfill(2003:2022, plt(3,:), plt(1,:), 'g', 'None', 1, 0.3); hold on

plt = LS4_pct(1:end,:)'/1e3;
lg(4,:) = plot(2003:2022, plt(2,:),'b','LineWidth',1.5);  % LS from 2005:2022 (constant covegare)
jbfill(2003:2022, plt(3,:), plt(1,:), 'b', 'None', 1, 0.3); 


set(gca,'fontsize',12);
ylabel('Annual lives saved (in thousands)','fontsize',14);
xlabel('Year','fontsize',12);
legend(lg(1:4,:), 'Relative to null services(from 2000)','Relative to constant coverage (from 2000)','Relative to null service(from 2003)','Relative to constant coverage (from 2003)','Location','NorthEast');
yl = ylim; yl(1) = 0; ylim(yl);
xlim([2000,2022])



%%%%%%%%%%%%%%%
% Create JSON data file 
s = struct('data',data, 'r', r, 'p',p, 'prm',prm, 'xs',xs);
jsonencode(s);

tmp = jsonencode(s);
fid = fopen('JSON_data.json', 'w');
fprintf(fid, '%s', tmp);
fclose(fid);


