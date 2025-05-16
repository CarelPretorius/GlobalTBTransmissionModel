clear all; %load calibres2; xsto = xsto4;
         
load calibres2; xsto = xsto3; %
load Model_setup;

ix0 = 2e4; nx = 150; dx = round((size(xsto,1)-ix0)/nx); 
xs = xsto(ix0:dx:end,:,1);
pct_xs= prctile(xs,[2.5,50,97.5])  % Estimated parameters

% % --- Get the Regional notifications
tend1   = 2024;  % Intervention starts except vaccination  
tend1b  = 2028;  % Vaccination starts
tend2   = 2035+1;  % End date for simulation

% --- Do the simulations 
opts = odeset('NonNegative',[1:i0.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);

r0 = r; p0 = p;
inct = [];
noti = [];

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    
    [r,p] = alloc_parameters(xs(ii,:),prm.r,prm.p,xi);
    % Get the initial conditions
    [out, aux] = obj(xs(ii,:));
    init = aux.soln(end,1:end);

    % Extend initial condition for the full model with vaccine compartments  
    inds1 = intersect(s0.U, s0.v0);
    inds2 = intersect([s0.Lf, s0.Ls], s0.v0);
    inds3 = [s0.Isc, s0.I, s0.E, s0.Rlo, s0.Rhi, s0.R];
    inds4 = [s0.Dx,s0.Tx,s0.Tx2];
    inds5 =  s0.Tx2(end)+[1:9];  % Auxiliary terms
    ind = [inds1, inds2, inds3,inds4,inds5];
   
    % To allocate results from a model without vaccine to input vector for full model with vaccine
    init0 = zeros(1,i0.nx);
    init0(ind) = init;   

    tref = [2022:1:tend2];
   
    % Interventions
%    PPM; Imp Dx+routine TB service; Upstream case finding; Casefinding subclinical;
%    PT; Upfront DST; Shorter duration of 2ndline; Imp 2ndline; Vaccination
    
% --- Now simulate without intervention (baseline)
     p.pu = p.pu_2022; r.cs = r.cs_2022; r.cs2 = r.cs2_2022; p.MDR_rec = p.MDR_rec2022;
     prm.growth = prm.gth(4);
     M0 = make_model3(p, r, i0, s0, gps0, prm); 
 
     geq = @(t,in) goveqs_scaleup(t, in, r, M0, M0, [2022 2023], i0, prm, sel0, agg0);
    [t0, soln0] = ode15s(geq, [2022:tend2], init0, odeset('NonNegative',[1:i0.nstates]));
       
    soln  = soln0;
    t     = t0;
    soln0 = interp1(t, soln, tref);
     
    % Find the initial condition at the starting point for the interventions
    geq = @(t,in) goveqs_scaleup(t, in, r, M0, M0, [2022 2023], i0, prm, sel0, agg0);
    [ta, solna] = ode15s(geq, [2022:tend1], init0, odeset('NonNegative',[1:i0.nstates]));
    initb = solna(end,:);
   
    % --- Next simulate with interventions; PPM
    r1 = r; p1 = p;
    p1.pu  = p.pu + (1- p.pu)*0.35; 
     
    M1 = make_model3(p1, r1, i0, s0, gps0,prm);
      
    geq = @(t,in) goveqs_scaleup(t, in, r, M0, M1, tend1 + [0  5], i0, prm, sel0, agg0);
    [tb, solnb] = ode15s(geq, [tend1 tend2], initb, opts);
  
    soln  = [solna; solnb(2:end,:)]; 
    t     = [ta;    tb(2:end)];
    soln1 = interp1(t, soln, tref); 
    
    %%%%%%
        % --- Next simulate with interventions; PPM+Imp Dx+ Routine TB service
    r1 = r; p1 = p;
    p1.pu  = p.pu + (1- p.pu)*0.35; 
    p1.Dx(1) = 0.95;         
    p1.MDR_rec = [0.5,0]; 
    r1.Tx2 = 12/9;                             % Shorter Tx2 regimn (9 months)
    p1.cure2 = [0.7, 0]; 
    del = 0.1;                                   % Percent increase of treatment completion. 
    p_Tx_complete(1) = min(xs(12)*(1+del),1);    % Estimated treatment completion in public sector
    r1.default(1) = r.Tx*(1-p_Tx_complete(1))./p_Tx_complete(1);

    M1 = make_model3(p1, r1, i0, s0, gps0,prm);
      
    geq = @(t,in) goveqs_scaleup(t, in, r, M0, M1, tend1 + [0  5], i0, prm, sel0, agg0);
    [tb, solnb] = ode15s(geq, [tend1 tend2], initb, opts);
  
    soln  = [solna; solnb(2:end,:)]; 
    t     = [ta;    tb(2:end)];
    soln2 = interp1(t, soln, tref);
   
   % --- Next simulate with interventions; + Upstream case finding (symptomatic)
    r2 = r1; p2 = p1;
    age_gp = 1:5;  % choose specific age group or age-groups
    csfac = 1.3;   % Choose any positive number to increase or decrease of care-seeking rate 
    r2.cs(age_gp)  = max(r.cs_1997(age_gp),r.cs(age_gp)*csfac);
    r2.cs2(age_gp) = max(r.cs2_1997(age_gp),r.cs2(age_gp)*csfac);

    M1 = make_model3(p2, r2, i0, s0, gps0,prm);
      
    geq = @(t,in) goveqs_scaleup(t, in, r, M0, M1, tend1 + [0  5], i0, prm, sel0, agg0);
    [tb, solnb] = ode15s(geq, [tend1 tend2], initb, opts);
  
    soln  = [solna; solnb(2:end,:)]; 
    t     = [ta;    tb(2:end)];
    soln3 = interp1(t, soln, tref);
    
       % --- Next simulate with interventions; + Upstream case finding (asymptomatic)
    r2 = r1; p2 = p1;
    age_gp = 1:5;  % choose specific age group or age-groups
    csfac = 1.3;   % Choose any positive number to increase or decrease of care-seeking rate 
    r2.cs(age_gp)  = max(r.cs_1997(age_gp),r.cs(age_gp)*csfac);
    r2.cs2(age_gp) = max(r.cs2_1997(age_gp),r.cs2(age_gp)*csfac);
    r2.cs3(age_gp) = max(r.cs_1997(age_gp),r.cs(age_gp)*csfac);
   
    M1 = make_model3(p2, r2, i0, s0, gps0,prm);
      
    geq = @(t,in) goveqs_scaleup(t, in, r, M0, M1, tend1 + [0  5], i0, prm, sel0, agg0);
    [tb, solnb] = ode15s(geq, [tend1 tend2], initb, opts);
  
    soln  = [solna; solnb(2:end,:)]; 
    t     = [ta;    tb(2:end)];
    soln4 = interp1(t, soln, tref);
    
     % --- Next simulate with interventions; + Preventive therapy
    r2 = r1; p2 = p1;
    age_gp = 3:5;   % choose specific age group or age-groups
    csfac = 1.3;    % Choose any positive number to increase or decrease of care-seeking rate 
    r2.cs(age_gp)  = max(r.cs_1997(age_gp),r.cs(age_gp)*csfac);
    r2.cs2(age_gp) = max(r.cs2_1997(age_gp),r.cs2(age_gp)*csfac);
    r2.cs3(age_gp) = max(r.cs_1997(age_gp),r.cs(age_gp)*csfac);
    qq = 0.05;                                                     % Preventive therapy among household contacts independent of hiv status
    r2.progression(age_gp,:)  = r.progression(age_gp,:)*(1-qq);    % Replace by this  --> p.PT_PLHIV = [1 1 (1-0.6)]
    r2.reactivation(age_gp,:) = r.reactivation(age_gp,:)*(1-qq);
    p2.PT_PLHIV = [1 1 (1-0.6)];                                    % PT among PLHIV reduces rate by 60%

    M1 = make_model3(p2, r2, i0, s0, gps0,prm);
      
    geq = @(t,in) goveqs_scaleup(t, in, r, M0, M1, tend1 + [0  5], i0, prm, sel0, agg0);
    [tb, solnb] = ode15s(geq, [tend1 tend2], initb, opts);
  
    soln  = [solna; solnb(2:end,:)]; 
    t     = [ta;    tb(2:end)];
    soln5 = interp1(t, soln, tref);
   
    % --- Next simulate with vaccination 
    % Find vaccination starting point
    [tv, solnv] = ode15s(geq, [tend1 tend1b], initb, opts);
    initv = solnv(end,:);
   
    % Choose suitable rate of vaccination and vaccine efficacy
    r3 = r2; p3 = p2;
    p3.VE = [0.6  0.6];                    % Efficacy--> [infection prevention, disease prevention]
    r3.vacc = [0 0 3 3 3];                 % Vaccination rate (chosen)
    M3 = make_model3(p3, r3, i0, s0, gps0, prm); % With vaccine
      
    geq = @(t,in) goveqs_scaleup(t, in,r, M1, M3, tend1b + [0  5], i0, prm, sel0, agg0); 
    [tc, solnc] = ode15s(geq, [tend1b tend2], initv, opts);
  
    soln  = [solna; solnv(2:end,:); solnc(2:end,:)];
    t     = [ta;    tv(2:end);      tc(2:end)];
     
    soln6 = interp1(t, soln, tref);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    % --- Bring them all together (baseline, interventions other than vaccine, interventions including vaccine)
   allmat = cat(3, soln0,soln1, soln2,soln3,soln4,soln5,soln6); 
   
    pop = sum(allmat(1:end-1,1:i0.nstates,1),2);  %2022:2035 (1:end -> 2022:2036)
    
    dallmat = diff(allmat,[],1);
    
    inct(:,:,ii)   = squeeze(dallmat(:,i0.aux.inc(1),:)./pop)*1e5;      % Total incidence
    %inch1(:,:,ii)  = squeeze(dallmat(:,i0.aux.inc(2),:)./pop)*1e5;      % HIV incidence
    mdr(:,:,ii)    = squeeze(dallmat(:,i0.aux.inc(3),:)./pop)*1e5;      % MDR incidence
    morth0(:,:,ii) = squeeze(dallmat(:,i0.aux.mort(1),:)./pop)*1e5;     % hiv-negative mortality
    %morth1(:,:,ii) = squeeze(dallmat(:,i0.aux.mort(2),:)./pop)*1e5;     % hiv-positive mortality
    noti(:,:,ii)   = squeeze(sum(dallmat(:,i0.aux.noti,:),2)./pop)*1e5; % Total notification
    notimdr(:,:,ii)= squeeze(dallmat(:,i0.aux.noti2,:)./pop)*1e5;       % MDR Treatment initiation

    
end
fprintf('\n'); 
    
inc_pct  = permute(prctile(inct,[2.5,50,97.5],3),[2,1,3]);  %2022 - 2035
mdr_pct  = permute(prctile(mdr*1.64,[2.5,50,97.5],3),[2,1,3]);
mrt_pct  = permute(prctile(morth0,[2.5,50,97.5],3),[2,1,3]);
noti_pct = permute(prctile(noti,[2.5,50,97.5],3),[3,1,2,4]);              % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario
notimdr_pct = permute(prctile(notimdr,[2.5,50,97.5],3),[3,1,2,4]);        % Dims: 1.Lo/Md/Hi, 2.Month, 3.Pu/Pr, 4.Scenario
 
 save res_forward2;   %Save output

 Figure_final
 

 
