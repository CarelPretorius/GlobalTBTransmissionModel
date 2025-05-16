load res_forward


%%%% 4.1 GENERAL EPIDEMIOLOGICAL BENCHMARKS

% Estimation of 'Cum_inc_adult_0_5'
tem1 = sum((LTBI_recent(1:5,1)).*r.progression(1),1);
tem2 = LTBI_recent(1,1);
tem3 = sum((LTBI_remote(1:5,1)).*r.reactivation(1),1);
tem4 = LTBI_remote(1,1);

Cum_inc_adult_0_5Lf = prctile(tem1*100./tem2,[2.5,50,97.5])
Cum_inc_adult_0_5Ls = prctile(tem3*100./tem4,[2.5,50,97.5])
Cum_inc_adult_0_5 = prctile(((tem1+tem3)*100./(tem2+tem4)),[2.5,50,97.5])


% Estimation of 'Ann_inc_adult_5p'
tem1 = LTBI_recent(1,1).*r.progression(1);
tem2 = LTBI_recent(1,1);
tem3 = LTBI_remote(1,1).*r.reactivation(1);
tem4 = LTBI_remote(1,1);

Ann_inc_adult_5p = prctile(((tem1+tem3)*100./(tem2+tem4)),[2.5,50,97.5])

% CFR: BRR_Rpt_TB_mort_no_Tx
Rpt_TB_mort_no_Tx = prctile(xs(xi.r_mort_TB(1))*100./(xs(xi.r_mort_TB(1))+xs(xi.r_self_cure)+r.mort),[2.5,50,97.5])

% % % Duration of Active TB in absense of treatment (year)
% Rpt_TB_dur_no_Tx1 = prctile(1./(xs(xi.r_cs)),[2.5,50,97.5])   % sx(:,6) -> r.cs
% or the below????
Rpt_TB_dur_no_Tx = prctile(1./(xs(:,xi.r_cs)+xs(:,xi.r_mort_TB(1))+xs(:,xi.r_self_cure)+r.mort),[2.5,50,97.5]) % xs(:,8)--> mu_TB; xs(:,18) --> self cure        

% Rpt_Part_imm_prior_inf
Rpt_Part_imm_prior_inf = 100*p.imm(1)    %prctile((Sus2(1,1,:)./Sus1(1,1,:)),[2.5,50,97.5])



%%%% 4.2 COUNTRY-SPECIFIC EPIDEMIOLOGICAL BENCHMARKS

% TB incidence rate 2022-2035: 'tb_incid_level'
tb_incid_level = [tab.incy_pct(1,:)]

% tb_incid_trend in last 5 years (i.e in 2026 relative to 2022 at baseline
tb_incid_trend = (tab.incy_pct(1,1)-tab.incy_pct(1,5))*100./(5*tab.incy_pct(1,1)) % or, %((tab.incy_pct(1,5,2)/tab.incy_pct(1,1,2)).^0.2)*100-100

% TB mortality level 'tb_mort_level_h0' (hiv-negative)
tb_mort_level_h0 = [tab.mrty_pct(1,:)]

% TB mortality level 'tb_mort_level_h1' (hiv+negative)
tb_mort_level_h1 = [tab.mrthivy_pct(1,:)]

% tb_mort_trend (hiv-negative) in last 5 years (i.e in 2023 relative to 2019)
tb_mort_trend_h0 = (tab.mrty_pct(1,1)-tab.mrty_pct(1,5))*100./(5*tab.mrty_pct(1,1)) % or, %((tab.mrty_pct(1,5,2)/tab.mrty_pct(1,1,2)).^0.2)*100-100

% tb_mort_trend (hiv-positive) in last 5 years (i.e in 2023 relative to 2019)
tb_mort_trend_h1 = (tab.mrthivy_pct(1,1)-tab.mrthivy_pct(1,5))*100./(5*tab.mrthivy_pct(1,1)) % or, %((tab.mrthivy_pct(1,5,2)/tab.mrthivy_pct(1,1,2)).^0.2)*100-100

% tb prevalence 
tb_prev_level = [tab.prvty_pct(1,:)]

% MDR-TB prevalence
mdr_naiv_level = [tab.prvtmdry_pct(1,:)]

% Prevalence of MDR-TB among treatment experienced notified TB cases, in the most recent available year mdr_expd_level
mdr_expd_level = [tab.notimdry_pct(1,:,1)]

% Percent of incident TB cases among HIV positive individuals, in the most recent available year 
tb_incid_hiv = [tab.inchivy_pct(1,:)]

% TB-HIV prevalence
TBhiv_prev = [tab.prvthivy_pct(1,:)] 

% TB-HIV prevalence
hiv_prev = data.HIV_prev*100;


%%%% 4.4 ADDITIONAL STANDARD OUTPUTS
% Percentage of total population infected with latent M.tb infection (LTBI), in most recent year (%)
ltbi_prev = [tab.LTBI_pct(1,:)]./1e3  % Total population 1e5


% Percent of incident TB cases due to recent infection (M.tb infection or reinfection within the last 5 years),
tem1 = LTBI_recent(1,1).*r.progression(1);
tem3 = LTBI_remote(1,1).*r.reactivation(1);
recent_pct= prctile((tem1*100./(tem1+tem3)),[2.5,50,97.5])

% % ARTI
%arti = [tab.ARTI(1,:,1);tab.ARTI(1,:,2);tab.ARTI(1,:,3)]


% Average number of new M.tb infections / reinfections produced by an infectious case, in most recent year.
eff_cont_case = pct_xs(:,xi.r_beta)'

% Duration of active TB (ie to death, self-cure,or treatment), in most recent year
duration_tb = prctile(1./xs(:,xi.r_cs)+(xs(:,xi.r_mort_TB(1))+xs(:,xi.r_self_cure)+r.mort),[2.5,50,97.5]) % sx(:,6) -> r.cs; xs(:,8)--> mu_TB; xs(:,18) --> self cure        

%%%%  Care cascade [only include public sector]
% Percent of all incident TB cases that will access TB diagnosis,for current year
access_pct = 100*prctile(xs(:,xi.p_pu),[2.5, 50, 97.5])  %xs(:,5)--> p.p_pu

% Percent of all incident TB cases that will receive a positive TB diagnosis, for current year
diag_pct =  100*prctile(xs(:,xi.p_Dx(1)),[2.5, 50, 97.5])   %xs(:,10) --> p_Dx

% Percent of all incident TB cases that will be notified, for current year
notif_pct = [tab.noti_pct(1,:)]

% Percent of all incident TB cases will initiate treatment, for current year
init_pct = [tab.noti_pct(1,:)]

% Percent of all incident TB cases will complete a treatment regimen, for current year
complete_pct = 100*prctile(xs(:,xi.p_Tx_complete(1)),[2.5, 50, 97.5])    % xs(:,12) --> p_Tx_complete

% Percent of all incident TB cases that will be cured via treatment, for current year
r.default = r.Tx*(1-xs(:,xi.p_Tx_complete(1)))./xs(:,xi.p_Tx_complete(1));
cure_pct =  100*prctile(r.Tx./(r.default+r.Tx),[2.5,50,97.5])    

% Percent of notifications without TB (ie false positive fraction), for most recent year
false_pos_pct = -1   % ???

% Percent of TB cases treated in private sector, for most recent year
tx_private_pct = 100*prctile(1-xs(:,xi.p_pu),[2.5, 50, 97.5])


%----------------------------
%General Epidemiology
BRR_Rpt_Cum_inc_adult_0_5  = 1;
BRR_Rpt_Ann_inc_adult_5p   = 2;
BRR_Rpt_TB_mort_no_Tx      = 3;
BRR_Rpt_TB_dur_no_Tx       = 4;
BRR_Rpt_Part_imm_prior_inf = 5;
%----------------------------
%Country-Specific Epidemiological Benchmarks
BRR_Rpt_TB_incid_level = 6;
BRR_Rpt_TB_incid_trend  = 7;
BRR_Rpt_TB_mort_level   = 8;
BRR_Rpt_TB_mort_trend   = 9;
BRR_Rpt_TB_prev_level   = 10;
BRR_Rpt_MDR_naiv_level  = 11;
BRR_Rpt_MDR_expd_level  = 12;
BRR_Rpt_TB_incid_hiv    = 13;
BRR_Rpt_HIV_prev        = 14;
%---------------------------
%Country-specific Economic Benchmarks
BRR_Rpt_TB_spending = 15;
BRR_Rpt_TB_cost_FLD = 16;
%----------------------------
%Additional standard outputs
%-Epidemiology
BRR_Rept_LTBI_prev    = 17;
BRR_Rpt_Recent        = 18;
BRR_Rpt_ARI           = 19;
BRR_Rpt_Eff_cont_rate = 20;
BRR_Rpt_TB_Dur        = 21;
%----------------------------
%-Care Cascade
BRR_Rpt_Access_Num     = 22;
BRR_Rpt_Diag_Num       = 23;
BRR_Rpt_Notif_Num      = 24;
BRR_Rpt_Init_Num       = 25
BRR_Rpt_Complete_Num   = 26;
BRR_Rpt_TxCure_Num     = 27;
BRR_Rpt_FP_Num         = 28;
BRR_Rpt_Tx_private_Num = 29;

BRR_Rpt_Access_Pct_Inc     = 22;
BRR_Rpt_Diag_Pct_Inc       = 23;
BRR_Rpt_Notif_Pct_Inc      = 24;
BRR_Rpt_Init_Pct_Inc       = 25
BRR_Rpt_Complete_Pct_Inc   = 26;
BRR_Rpt_TxCure_Pct_Inc     = 27;
BRR_Rpt_FP_Pct_Inc         = 28;
BRR_Rpt_Tx_private_Pct_Inc = 29;

BRR_Rpt_Access_Pct_Risk     = 22;
BRR_Rpt_Diag_Pct_Risk       = 23;
BRR_Rpt_Notif_Pct_Risk      = 24;
BRR_Rpt_Init_Pct_Risk       = 25
BRR_Rpt_Complete_Pct_Risk   = 26;
BRR_Rpt_TxCure_Pct_Risk     = 27;
BRR_Rpt_FP_Pct_Risk         = 28;
BRR_Rpt_Tx_private_Pct_Risk = 29;

%-----------------------------
%Policy Projections
BRR_Rpt_Fut_incid_trend_sq  = 30;
BRR_Rpt_Fut_mort_trend_sq   = 31;
BRR_Rpt_Fut_incid_trend_max = 32;
BRR_Rpt_Fut_mort_trend_max  = 33;

BRR_Rpt_max                 = 33;
%------------------------------

BRRSheetNames{1}='LB';
BRRSheetNames{2}='Mean';
BRRSheetNames{3}='UB';

for I=1:3   % WHAT IS I here?

num_years=11;
BRR_Tab  = zeros(BRR_Rpt_max,num_years);   
BRR_Rows = cell(BRR_Rpt_max,1);   
    
%General Epidemiology
BRR_Tab(BRR_Rpt_Cum_inc_adult_0_5,1)=Cum_inc_adult_0_5(I);
BRR_Rows{BRR_Rpt_Cum_inc_adult_0_5,1}='Cum_inc_adult_0_5';
BRR_Tab(BRR_Rpt_Ann_inc_adult_5p,1)=Ann_inc_adult_5p(I);
BRR_Rows{BRR_Rpt_Ann_inc_adult_5p,1}='Ann_inc_adult_5p';
BRR_Tab(BRR_Rpt_TB_mort_no_Tx,1)=Rpt_TB_mort_no_Tx(I);
BRR_Rows{BRR_Rpt_TB_mort_no_Tx,1}='Rpt_TB_mort_no_Tx';
BRR_Tab(BRR_Rpt_TB_dur_no_Tx,1)= Rpt_TB_dur_no_Tx(I);
BRR_Rows{BRR_Rpt_TB_dur_no_Tx,1}='Rpt_TB_dur_no_Tx';
BRR_Tab(BRR_Rpt_Part_imm_prior_inf,1)=Rpt_Part_imm_prior_inf(1);
BRR_Rows{BRR_Rpt_Part_imm_prior_inf,1}='Rpt_Part_imm_prior_inf';

%Country-Specific Epidemiological Benchmarks
BRR_Tab(BRR_Rpt_TB_incid_level,1:num_years)=tb_incid_level(I,:);
BRR_Rows{BRR_Rpt_TB_incid_level,1}='tb_incid_level';
BRR_Tab(BRR_Rpt_TB_incid_trend,1)=tb_incid_trend(1);
BRR_Rows{BRR_Rpt_TB_incid_trend,1}='tb_incid_trend';
BRR_Tab(BRR_Rpt_TB_mort_level,1:num_years)=tb_mort_level_h1(I,:);
BRR_Rows{BRR_Rpt_TB_mort_level,1}='tb_mort_level_h1';
BRR_Tab(BRR_Rpt_TB_mort_trend,1)=tb_mort_trend_h1(1);
BRR_Rows{BRR_Rpt_TB_mort_trend,1}='tb_mort_trend_h1';
BRR_Tab(BRR_Rpt_TB_prev_level,1:num_years)=tb_prev_level(I,:);
BRR_Rows{BRR_Rpt_TB_prev_level,1}='tb_prev_level';
BRR_Tab(BRR_Rpt_MDR_naiv_level,1:num_years)=mdr_naiv_level(I,:);
BRR_Rows{BRR_Rpt_MDR_naiv_level,1}='mdr_naiv_level';
BRR_Tab(BRR_Rpt_MDR_expd_level,1:num_years)=mdr_expd_level(I,:);
BRR_Rows{BRR_Rpt_MDR_expd_level,1}='mdr_expd_level';
BRR_Tab(BRR_Rpt_TB_incid_hiv,1:num_years)=tb_incid_hiv(I,:);
BRR_Rows{BRR_Rpt_TB_incid_hiv,1}='tb_incid_hiv';
BRR_Tab(BRR_Rpt_HIV_prev,1)=hiv_prev(I);
BRR_Rows{BRR_Rpt_HIV_prev,1}='hiv_prev';

%Country-specific Economic Benchmarks
%TBD

%Additional standard outputs
%-Epidemiology
BRR_Tab(BRR_Rept_LTBI_prev,1:num_years)=ltbi_prev(1:num_years);
BRR_Rows{BRR_Rept_LTBI_prev,1}='ltbi_prev';
BRR_Tab(BRR_Rpt_Recent,1)=recent_pct(I);
BRR_Rows{BRR_Rpt_Recent,1}='recent_pct';
%BRR_Tab(BRR_Rpt_ARI,1)=arti(I);
BRR_Rows{BRR_Rpt_ARI,1}='arti';
BRR_Tab(BRR_Rpt_Eff_cont_rate,1)=eff_cont_case(I);
BRR_Rows{BRR_Rpt_Eff_cont_rate,1}='eff_cont_case';
BRR_Tab(BRR_Rpt_TB_Dur,1)=duration_tb(I);
BRR_Rows{BRR_Rpt_TB_Dur,1}='duration_tb';

%-Care Cascade
BRR_Tab(BRR_Rpt_Access_Pct_Inc,1)=access_pct(I);
BRR_Rows{BRR_Rpt_Access_Pct_Inc,1}='access_pct';
BRR_Tab(BRR_Rpt_Diag_Pct_Inc,1)=diag_pct(I);
BRR_Rows{BRR_Rpt_Diag_Pct_Inc,1}='diag_pct';
BRR_Tab(BRR_Rpt_Notif_Pct_Inc,1:num_years)=notif_pct(I,:);
BRR_Rows{BRR_Rpt_Notif_Pct_Inc,1}='notif_pct';
BRR_Tab(BRR_Rpt_Init_Pct_Inc,1:num_years)=init_pct(I,:);
BRR_Rows{BRR_Rpt_Init_Pct_Inc,1}='init_pct';
BRR_Tab(BRR_Rpt_Complete_Pct_Inc,1)=complete_pct(I);
BRR_Rows{BRR_Rpt_Complete_Pct_Inc,1}='complete_pct';
BRR_Tab(BRR_Rpt_TxCure_Pct_Inc,1)=cure_pct(I);
BRR_Rows{BRR_Rpt_TxCure_Pct_Inc,1}='cure_pct';
%BRR_Tab(BRR_Rpt_FP,1:num_years)=0;
BRR_Tab(BRR_Rpt_Tx_private_Pct_Inc,1)=tx_private_pct(I);
BRR_Rows{BRR_Rpt_Tx_private_Pct_Inc,1}='tx_private_pct';

xlswrite('BRRTableKenya.xls',BRR_Tab,BRRSheetNames{I});
if(I==1)
xlswrite('BRRTableKenya.xls',BRR_Rows,'BRRRows');
end

end

'table written'




