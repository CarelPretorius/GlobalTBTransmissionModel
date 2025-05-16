function [ BRR_Table] = Make_BRR_Table(params) %,indicators)

%Construct the BRR table, defined in Section 4 of:
%Supplementary Materials for “A new approach to improve the quality of mathematical
%modelling for country-level TB decision-making - development and piloting”.

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
BRR_Rpt_Access     = 22;
BRR_Rpt_Diag       = 23;
BRR_Rpt_Notif      = 24;
BRR_Rpt_Complete   = 25;
BRR_Rpt_FP         = 26;
BRR_Rpt_Tx_private = 27;
%-----------------------------
%Policy Projections
BRR_Rpt_Fut_incid_trend_sq  = 30;
BRR_Rpt_Fut_mort_trend_sq   = 31;
BRR_Rpt_Fut_incid_trend_max = 32;
BRR_Rpt_Fut_mort_trend_max  = 33;

BRR_Rpt_max                 = 27;
%------------------------------

prm=params;

% prm.p:
% progression: [0.0826 0.8260 0.1652]
% LTBI_stabil: [0.8720 0 0.8720]
% reactivation: [6.0000e-04 0.0600 0.0120]
% self_cure: [0.1667 0 0.1667]
% mort_TB: 0.1667
% relapse: [0.0320 0.1400 0.0015]
% mort: 0.0152
% Dx: 52
% Tx: 2
% Tx2: 0.5000
% default: 0.8169
% default2: [0.5417 5.0000e+05]
% pt: 0
% access: 0
% cs3: 0
% vacc: 0
% waning: 0
% 
% 
% prm.r: 
% imm: [0.8000 0 0.8000]
% MDR_rec2015: [1.0000e-03 0]
% Tx_init2: [0.8800 0]
% SL_trans: [0.8800 0]
% Tx_init: [1 1]
% default: 0.2900
% cure: [1 1]
% tsrsl: [0.4800 1.0000e-06]
% cure2: [0.5000 0]
% PT_PLHIV: [1 1 1]
% VE: [0 0]


pop_19_30 = [52573973,53771296,54985698,56215221,57458897,58714804,59981310,61257803,52544014,63838878,65141206,66449654]; 

num_years = size(pop_19_30,2);
BRR_Table=zeros(BRR_Rpt_max,num_years);

for I=1:BRR_Rpt_max %do    
for t=1:num_years
 
pop_t=pop_19_30(t);
inc_t=inc(t);
mort_t=mort(t);
prev_t=prev(t);

pop_t_5=0;
inc_t_5=0;
mort_t_5=0;

if(t>5)
    pop_t_5=pop(t-5); 
    inc_t_5=inc(t-5);
    mort_t_5=inc(t-5);
end    

value=0;
switch I
    
    %General Epidemiology
    case BRR_Rpt_Cum_inc_adult_0_5    
        value = prm.p.reactivation;
          
    case BRR_Rpt_Ann_inc_adult_5p    
        value = prm.p.reactivation;
    
    case BRR_Rpt_TB_mort_no_Tx
        value = prm.p.r_mort_TB(1);
   
    case BRR_Rpt_TB_dur_no_Tx
        value = 1/(prm.p.mort_TB+prm.p.self_cure);
    
    case BRR_Rpt_Part_imm_prior_inf
        value = -1;
    
    %Country-Specific Epidemiological Benchmarks
    case BRR_Rpt_TB_incid_level
        value=10^5*inc_t/pop_t;
        
    case BRR_Rpt_TB_incid_trend
        if( (t>5) && (inc_t_5>0) ) 
            value=( (inc_t/pop_t)/(inc_t_5/pop_t_5)-1 )^1/5;
        end
        
    case BRR_Rpt_TB_mort_level
        value=10^5*mort_t/pop_t;
        
    case BRR_Rpt_TB_mort_trend
        if( (t>5) && (mort_t_5>0) ) 
            value=( (mort_t/pop_t)/(mort_t_5/pop_t_5) -1 )^1/5;
        end 
        
    case BRR_Rpt_TB_prev_level
        value=10^5*prev_t/pop_t; 
                
    case BRR_Rpt_MDR_naiv_level
        value=100*mdr_prev_t/prev_tb_no_tx;

    case BRR_Rpt_MDR_expd_level
        value=100*mdr_prev_t/prev_tb_had_tx;

    case BRR_Rpt_TB_incid_hiv
        value=100*inc_hiv_t/inc_t;   
        
    case BRR_Rpt_HIV_prev
        value=100*pop_hiv_t/pop_t; 
        
            
        
  otherwise
    value=-1;    
end

BRR_Table(BRR_Rpt_Cum_inc_adult_0_5,t)=value;

end



end

