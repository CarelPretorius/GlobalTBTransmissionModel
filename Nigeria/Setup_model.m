clear all;

load('../Collect_data/NGA.mat')

% Age-group of Population in each compartment [0-4,5-9,10-14,15-65,65+]   
prm.N = [0.2 0.2 0.2 0.2 0.2];                     % Initial population (assumption)
prm.cm = data.contact_mat;                         % Contact matrix correcting with population weight
prm.gth = [0   0.03   0.035    0.044];            %[0,data.gth];  % population growth rate 
r.mort =  data.age_mort';                           % Age specific mortality rate 
r.ageing = [0.1762    0.1702    0.1606    0.0050];  %data.ageing';  % Ageing rate  

% --- Get calibration targets ---------------------------------------------
popn              = data.pop00_22(end);           % 2022
data.sym          = [0.30 0.50 0.70];             % proportion of prevalent TB that has symptoms % 0.5[0.36, 0.8]subclinical
data.inc2000      = data.inc_all(1,:); 
data.inc2015      = data.inc_all(end-7,:); 
data.inc2022      = data.inc_all(end,:); 
data.mort_H02000  = data.mort_H0(1,:);
data.mort_H02022  = data.mort_H0(end,:);
data.noti_cu      = sum(data.noti,1)*[0.8 1 1.2];             % from 2000 to 2022
data.noti2022     = [0.8 1 1.2]*data.noti(end);
data.mdr2015      = data.mdr2015; 
data.mdr2022      = data.mdr2022;
data.mdriniTX     = [0.8 1 1.2]*data.mdriniTX(end); 
data.inc_h1       = data.inc_h1;
data.mort_H12022  = data.mort_H1;

% Get the HIV data                                                         
prm.ART_start = data.ART_start;                     %ART_start;
data.ART_covg =  [0.85 1 1.15]*data.ART_covg;       %ARTcovg_2022;
data.HIV_prev = [0.85 1 1.15]*(data.HIV_prev./popn);%HIVprev_2022 (proportion of population)
HIV_incd = data.HIV_incd;                           %HIV incidence (proportion of population) 1980-2022


% Extrapolate HIV incidence by 20 years
ys1 = HIV_incd;                                     %1980-2022
xs1 = 1:length(ys1);
xs2 = 1:length(ys1)+20;                             % 2022-2042
ys2 = interp1(xs1,ys1,xs2,'linear','extrap');       % linier trend for the next 21 years
prm.rHIV = ys2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gps1.age  = {'a1','a2','a3','a4','a5'}; % Simplest model
gps1.vacc = {'v0'}; 
gps1.hiv  = {'h0'}; 
gps1.strains = {'ds'}; 
gps1.provs = {'pr'};  

gps2.age  = {'a1','a2','a3','a4','a5'};  % With MDR
gps2.vacc = {'v0'};  
gps2.hiv  = {'h0'}; 
gps2.strains = {'ds','mdr'};
gps2.provs = {'pr'}; 

gps.age  = {'a1','a2','a3','a4','a5'};  % With MDR, pu/pr and HIV
gps.vacc = {'v0'};                
gps.hiv  = {'h0','h1','hart'};
gps.strains = {'ds','mdr'};
gps.provs = {'pu','pr'};

gps0.age  = {'a1','a2','a3','a4','a5'};  % Full model including vaccine
gps0.vacc = {'v0','v1','v2'};            % After model calibration we can add these stages
gps0.hiv  = {'h0','h1','hart'};
gps0.strains = {'ds','mdr'};
gps0.provs = {'pu','pr'};

states0 = {'U'};                              % structured by hiv and vaccine status 
states1 = {'Lf','Ls'};                        % structured by hiv, vaccine and strain 
states2 = {'Isc','I','E','Rlo','Rhi','R'};    % structured by hiv and strain 
states3 = {'Dx','Tx','Tx2'};                  % structured by hiv, strain and provider type

% For equilibrium model without mdr/hiv/vaccine/public sector 
[i1, s1, d1, lim1] = get_addresses({states0, gps1.age, gps1.hiv, gps1.vacc}, [], [], [], 0);
[i1, s1, d1, lim1] = get_addresses({states1, gps1.age, gps1.hiv, gps1.vacc, gps1.strains}, i1, s1, d1, lim1);
[i1, s1, d1, lim1] = get_addresses({states2, gps1.age, gps1.hiv, gps1.strains}, i1, s1, d1, lim1);
[i1, s1, d1, lim1] = get_addresses({states3, gps1.age, gps1.hiv, gps1.strains, gps1.provs}, i1, s1, d1, lim1);
d1 = char(d1);

% Indtroducing mdr with the previous small model
[i2, s2, d2, lim2] = get_addresses({states0, gps2.age, gps2.hiv, gps2.vacc}, [], [], [], 0);
[i2, s2, d2, lim2] = get_addresses({states1, gps2.age, gps2.hiv, gps2.vacc, gps2.strains}, i2, s2, d2, lim2);
[i2, s2, d2, lim2] = get_addresses({states2, gps2.age, gps2.hiv, gps2.strains}, i2, s2, d2, lim2);
[i2, s2, d2, lim2] = get_addresses({states3, gps2.age, gps2.hiv, gps2.strains, gps2.provs}, i2, s2, d2, lim2);
d2 = char(d2);

% For the full model without vaccination
[i, s, d, lim] = get_addresses({states0, gps.age, gps.hiv, gps.vacc}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states1, gps.age, gps.hiv, gps.vacc, gps.strains}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states2, gps.age, gps.hiv, gps.strains}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states3, gps.age, gps.hiv, gps.strains, gps.provs}, i, s, d, lim);
d = char(d);

% For the full model with vaccination
[i0, s0, d0, lim0] = get_addresses({states0, gps0.age, gps0.hiv, gps0.vacc}, [], [], [], 0);
[i0, s0, d0, lim0] = get_addresses({states1, gps0.age, gps0.hiv, gps0.vacc, gps.strains}, i0, s0, d0, lim0);
[i0, s0, d0, lim0] = get_addresses({states2, gps0.age, gps0.hiv, gps0.strains}, i0, s0, d0, lim0);
[i0, s0, d0, lim0] = get_addresses({states3, gps0.age, gps0.hiv, gps0.strains, gps.provs}, i0, s0, d0, lim0);
d0 = char(d0);

s1.infectious      = [s1.Isc, s1.I, s1.E, s1.Dx];
s1.infectious_wosc = [s1.I, s1.E, s1.Dx];    % Infectious compartment without subclinical (contributing on TB death)
s1.prevalent       = unique([s1.infectious, [s1.Tx,s1.Tx2]]);

s2.infectious      = [s2.Isc, s2.I, s2.E, s2.Dx, intersect(s2.Tx, s2.mdr)];
s2.infectious_wosc = [s2.I, s2.E, s2.Dx, intersect(s2.Tx, s2.mdr)];    % Infectious compartment without subclinical (contributing on TB death)
s2.prevalent       = unique([s2.infectious, [s2.Tx,s2.Tx2]]);

s.infectious      = [s.Isc, s.I, s.E, s.Dx, intersect(s.Tx, s.mdr)];
s.infectious_wosc = [s.I, s.E, s.Dx, intersect(s.Tx, s.mdr)];    % Infectious compartment without subclinical (contributing on TB death)
s.prevalent       = unique([s.infectious, [s.Tx,s.Tx2]]);

s0.infectious      = [s0.Isc, s0.I, s0.E, s0.Dx, intersect(s0.Tx, s0.mdr)];
s0.infectious_wosc = [s0.I, s0.E, s0.Dx, intersect(s0.Tx, s0.mdr)];    % Infectious compartment without subclinical (contributing on TB death)
s0.prevalent       = unique([s0.infectious, [s0.Tx,s0.Tx2]]);


% --- Include the auxiliaries ---------------------------------------------
names = {'inc','noti','noti2','mort'};
lgths =     [3,     3,      1,     2];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
  
for ii = 1:length(names)
    inds0 = lim0 + [1:lgths(ii)];
    i0.aux.(names{ii}) = inds0;
    lim0 = inds0(end);
end

i1.nx = lim1;
i2.nx = lim2;
i.nx  = lim;
i0.nx = lim0;

% --- Make aggregators and selectors --------------------------------------

% Selectors for the incidence
tmp = zeros(3,i.nstates);
tmp(1,s.Isc) = 1;                                            % All
tmp(2,[intersect(s.Isc,s.h1),intersect(s.Isc,s.hart)]) = 1;  % hiv-TB
tmp(3,intersect(s.Isc,s.mdr)) = 1;                           % mdr-TB 
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Isc,:) = 1;
tmp(s.h1,s.h0) = 0; tmp(s.hart,s.h1) = 0; tmp(s.h0,s.h1) = 0; tmp(s.h1,s.hart) = 0;
tmp(s.ds,s.mdr)= 0; tmp(s.mdr,s.ds) = 0;   % ADD
sel.inc = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(intersect(s.Tx,s.mdr),intersect(s.Tx,s.ds)) = 1;
sel.acqu = sparse(tmp - diag(diag(tmp)));

% Selectors for public sector notifications
tmp = zeros(3,i.nstates);
tmp(1,intersect(intersect([s.Tx,s.Tx2],s.h0),s.pu)) = 1;  % public sector notification
tmp(2,intersect(intersect([s.Tx,s.Tx2],s.h1),s.pu)) = 1;
tmp(3,intersect(intersect([s.Tx,s.Tx2],s.hart),s.pu)) = 1;
agg.noti = sparse(tmp);

tmp = zeros(i.nstates);
tmp([s.Tx,s.Tx2],:) = 1;
tmp(s.h1,s.h0) = 0; tmp(s.hart,s.h1) = 0; tmp(s.h0,s.h1) = 0; tmp(s.h1,s.hart) = 0;
tmp(s.ds,s.mdr)= 0; tmp(s.mdr,s.ds) = 0;
sel.noti = tmp - diag(diag(tmp));

% Selectors for second-line notifications
tmp = zeros(i.nstates);
tmp(s.Tx2,:) = 1;
tmp(s.h1,s.h0) = 0; tmp(s.hart,s.h1) = 0; tmp(s.h0,s.h1) = 0; tmp(s.h1,s.hart) = 0;
tmp(s.ds,s.mdr)= 0; tmp(s.mdr,s.ds) = 0;
sel.noti2 = tmp - diag(diag(tmp));

% Selectors for the incidence from the full model with vaccine
tmp = zeros(3,i0.nstates);
tmp(1,s0.Isc) = 1;                                            % All
tmp(2,[intersect(s0.Isc,s0.h1),intersect(s0.Isc,s0.hart)]) = 1;  % hiv-TB
tmp(3,intersect(s0.Isc,s0.mdr)) = 1;                           % mdr-TB 
agg0.inc = sparse(tmp);

tmp = zeros(i0.nstates);
tmp(s0.Isc,:) = 1;
tmp(s0.h1,s0.h0) = 0; tmp(s0.hart,s0.h1) = 0; tmp(s0.h0,s0.h1) = 0; tmp(s0.h1,s0.hart) = 0;
tmp(s0.ds,s0.mdr)= 0; tmp(s0.mdr,s0.ds) = 0;   % ADD
sel0.inc = tmp - diag(diag(tmp));

tmp = zeros(i0.nstates);
tmp(intersect(s0.Tx,s0.mdr),intersect(s0.Tx,s0.ds)) = 1;
sel0.acqu = sparse(tmp - diag(diag(tmp)));

% Selectors for public sector notifications
tmp = zeros(3,i0.nstates);
tmp(1,intersect(intersect([s0.Tx,s0.Tx2],s0.h0),s0.pu)) = 1;  % public sector notification
tmp(2,intersect(intersect([s0.Tx,s0.Tx2],s0.h1),s0.pu)) = 1;
tmp(3,intersect(intersect([s0.Tx,s0.Tx2],s0.hart),s0.pu)) = 1;
agg0.noti = sparse(tmp);

tmp = zeros(i0.nstates);
tmp([s0.Tx,s0.Tx2],:) = 1;
tmp(s0.h1,s0.h0) = 0; tmp(s0.hart,s0.h1) = 0; tmp(s0.h0,s0.h1) = 0; tmp(s0.h1,s0.hart) = 0;
tmp(s0.ds,s0.mdr)= 0; tmp(s0.mdr,s0.ds) = 0;
sel0.noti = tmp - diag(diag(tmp));

% Selectors for second-line notifications
tmp = zeros(i0.nstates);
tmp(s0.Tx2,:) = 1;
tmp(s0.h1,s0.h0) = 0; tmp(s0.hart,s0.h1) = 0; tmp(s0.h0,s0.h1) = 0; tmp(s0.h1,s0.hart) = 0;
tmp(s0.ds,s0.mdr)= 0; tmp(s0.mdr,s0.ds) = 0;
sel0.noti2 = tmp - diag(diag(tmp));

% --- Define variables and ranges -----------------------------------------
names = {'r_beta','rfbeta_mdr','rfbeta_hiv','r_sym','p_pu','r_cs','rf_cs2','r_mort_TB','p_Dx','p_Tx_complete','p_MDRrec2022', 'r_MDR_acqu','r_ART_init','r_HIV_mort','r_self_cure','p_HIV_relrate', 'rf_p_imm', 'rf_cs2022', 'r_betared'};   %'rf_LTBI_cure', <--- Updated
lgths =        [1,           1,           1,     1,     1,     1,      1,          2,     2,              2,             1,            1,            1,           1,            1,               1,           1,          1,         1];

xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
%xi.calib = xi.r_HIV_mort;
xi.nx = lim;

bds = [];
bds(xi.r_beta,:)        = [0 40]; 
bds(xi.rfbeta_mdr,:)    = [0 1];
bds(xi.rfbeta_hiv,:)    = [0 1];
bds(xi.r_sym,:)         = [0.1 100];
bds(xi.p_pu,:)          = [0 1];
bds(xi.r_cs,:)          = [0.1 10]; 
bds(xi.rf_cs2,:)        = [1 40]; 
%bds(xi.r_Tx_init,:)    = [0.1 10];
bds(xi.r_mort_TB,:)     = 1/6*[0 2; 0 100];
bds(xi.p_Dx,:)          = [0.75 0.9; 0.1  0.7]; 
bds(xi.p_Tx_complete,:) = [0.75 0.95; 0.4 0.8];
bds(xi.p_MDRrec2022,:)  = [0  1];
bds(xi.r_MDR_acqu,:)    = [0  0.06];
bds(xi.r_ART_init,:)    = [0 10];
bds(xi.r_HIV_mort,:)    = [0 10];
bds(xi.r_self_cure,:)   = 1/6*[0.85 1.15];
bds(xi.p_HIV_relrate,:) = [1 100];
bds(xi.rf_p_imm,:)      = [0.75, 1.25]; %[0.1 0.75];
bds(xi.rf_cs2022,:)     = [1 10];
%bds(xi.rf_LTBI_cure,:)  = [0.75, 1.25];
bds(xi.r_betared,:)     = [0.8 1];

prm.bounds = bds';

% --- Define baseline parameter values ------------------------------------
% Natural history (h0, h1, hart)
r.progression  = 0.0826*[1 10 10*0.2];
r.LTBI_stabil  = 0.872*[1 0 1];
r.LTBI_cure    = 0.028; 
r.reactivation = 0.0006*[1 100 100*0.2];
r.relapse      = [0.032 0.14 0.0015];
p.imm          = [0.8 0 0.8];

% Diagnosis stage
r.Dx             = 52;

p.MDR_rec2015 = [0.001 0]; % CHECK for each country [pu pr]
p.Tx_init2    = [0.88 0];  % CHECK for each country [pu pr]
p.SL_trans    = [0.88 0];  % [pu pr]

% Treatment stage
p.Tx_init     = [1 1];
r.Tx      = 2;
r.Tx2     = 0.5;   
% p.default = 1 - 0.71;
% r.default = r.Tx*p.default./(1-p.default);
p.cure    = [1 1];

p.tsrsl     = [0.48 1e-6];                % Second-line treatment success rate in pu and pr
r.default2  = r.Tx2*(1-p.tsrsl)./p.tsrsl;
p.cure2     = [0.5 0];                 % Should I consider a different value for those with HIV?

% Intervention parameter
% r.pt          = 0;
p.PT_PLHIV    = [1 1 1]; % With PT intervenion among PLHIv [1 1 (1-0.6)]
r.access      = 0;
r.cs3         = 0;          % case finding rate among subclinical TB
p.VE          = [0 0];
r.vacc        = [0 0 0 0 0]; 
r.waning      = 0;


% Bring them all together
prm.p = p; prm.r = r;
ref.i = i;   ref.s = s;   ref.d = d;  ref.xi = xi;
ref.i1 = i1; ref.s1 = s1; ref.d1 = d1; 
ref.i2 = i2; ref.s2 = s2; ref.d2 = d2; 


show = 0;
f0 = get_distribution_fns(data.sym,         'lognorm', show);
f1 = get_distribution_fns(data.inc2000,     'lognorm', show);  
f2 = get_distribution_fns(data.inc2022,     'lognorm', show);
f3 = get_distribution_fns(data.mort_H02000, 'lognorm', show);
f4 = get_distribution_fns(data.mort_H02022, 'lognorm', show);
f5 = get_distribution_fns(data.noti_cu,     'lognorm', show);
f6 = get_distribution_fns(data.noti2022,    'lognorm', show);
f7 = get_distribution_fns(data.mdr2015,     'lognorm', show);
f8 = get_distribution_fns(data.mdr2022,     'lognorm', show);
f9 = get_distribution_fns(data.mdriniTX,    'lognorm', show); 
f10 = get_distribution_fns(data.HIV_prev,   'lognorm', show);
f11 = get_distribution_fns(data.inc_h1,     'lognorm', show);
f12 = get_distribution_fns(data.ART_covg,   'lognorm', show);
f13 = get_distribution_fns(data.mort_H12022,'lognorm', show);

lhd.fn  = @(sym,inc2000, inc2022, mort_H02000,mort_H02022,noti_cu, noti2022, mdr2015, mdr2022, mdriniTX,HIV_prev,inc_h1,ART_covg, mort_H12022)f0(sym)+f1(inc2000)+ 10*f2(inc2022)+f3(mort_H02000)+ 10*f4(mort_H02022)+ f5(noti_cu)+10*f6(noti2022)+f7(mdr2015)+f8(mdr2022)+f9(mdriniTX)+f10(HIV_prev)+f11(inc_h1)+f12(ART_covg)+f13(mort_H12022);

lhd.sgn = -Inf;


save Model_setup;

