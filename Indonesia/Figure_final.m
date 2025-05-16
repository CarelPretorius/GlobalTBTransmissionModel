clear all; 
load res_forward2; % Load outputs from the, 'Simulate_forward.m'

fs = 14; ms = 24; lw = 1.6; ts = 0.1;

years = {'2022','2023','2024','2025','2026','2027','2028','2029','2030','2031','2032','2033','2034','2035'};

% --- Incidence and mortality projections ---------------------------------
ff = figure;
skip = 12; 
cols = linspecer(8);
pmat = cat(4,inc_pct,mrt_pct); 
pmat = permute(pmat,[3,2,1,4]);    % By 'permute' we can reorganise the matrix-order as required

tis  = {'Annual incidence','Annual mortality'}; 
sis  = [1,2]; % 1;
target1 = [(370)*(1-0.8),48*(1-0.9)]; % [incd_pct(end-4,2),mor_pct(end-4,2)];2030 target relative to 2015 
target2 = [(370)*(1-0.9),48*(1-0.95)]; %[incd_pct(end-4,2),mor_pct(end-4,2)]: 2035 target relative to 2015 
for is = 1:2  %2
    subplot(1,2,sis(is));
    for ii = 1:7                % All intervetions
        plt = pmat(:,:,ii,is); 
        pl1(ii,:) = plot(plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
        jbfill(1:size(plt,2),plt(3,:),plt(1,:),cols(ii,:),'None',1,ts); hold on; 
       % plot(1:size(plt,2),repmat(target1(is),1,size(plt,2)),'-.r','linewidth',lw); hold on;    
        plot(1:size(plt,2),repmat(target2(is),1,size(plt,2)),'-.k','linewidth',lw); hold on;    
    end
    
     xlim([1 size(plt,2)]);
     %ylim([0,50]);
     yl = ylim; yl(1) = 0; ylim(yl);
     % set(gca, 'YTick', [])
     ylabel('Annual rate per 100 000','fontsize',12);
     xlabel('Year','fontsize',12);
     legend(pl1,'Baseline','+ PPM','+ Improved diagnostics, routine TB services','+ Upstream case-finding','+ Detecting asymptomatic TB','+ Preventive therapy (risk group)','+ Vaccine','location','SouthEast','fontsize',9);
  
    skip = 2;
    set(gca,'fontsize',12,'XTick',1:skip:size(plt,2),'XTickLabel',years(1:skip:size(plt,2)));
   % xlim([1 8])
    xtickangle(60) %45
    title(tis{is});
    %title('Indonesia')
    txt1 = {'2035 End-TB targets'};
    scaleposition = [0.7, 1];
    text(1.5,target1(is)*scaleposition(is),txt1,'fontsize',12) 
    
end

set(ff,'Position',[451   308   640   538]);

sgtitle('Global Plan Scenario')


% TB incidence in 2023 and in 2035 (mean, lower, upper)
% Dim: = (intervention_stage, year, [lo, mid, hi])
int_layer = 3;
inc2023 = [inc_pct(1,2,1),inc_pct(1,2,2),inc_pct(1,2,3)];  % [mean, lower, upper]
inc2035 = [inc_pct(int_layer,end,1),inc_pct(int_layer,end,2),inc_pct(int_layer,end,3)];  % [mean, lower, upper]

mrt2023 = [mrt_pct(1,2,1),mrt_pct(1,2,2),mrt_pct(1,2,3)];  % [mean, lower, upper]
mrt2035 = [mrt_pct(int_layer,end,1),mrt_pct(int_layer,end,2),mrt_pct(int_layer,end,3)];  % [mean, lower, upper]

Get_reduction(inc2023, inc2035)
Get_reduction(mrt2023, mrt2035)

