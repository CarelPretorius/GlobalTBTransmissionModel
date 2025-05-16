function [out, aux] = get_objective(x, prm, ref, sel, agg, gps,gps1,gps2,calfn)     %India

r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;
%i0 = ref.i0; s0 = ref.s0; 
i1 = ref.i1; s1 = ref.s1; i2 = ref.i2; s2 = ref.s2;

mat = [prm.bounds(2,1:length(x))-x; x-prm.bounds(1,1:length(x))];

if min(mat(:)) < 0
    out = Inf*calfn.sgn;
    aux = [];
else
    
    [r,p] = alloc_parameters(x,r,p,xi);
    
    % --- Set up the necessary models -------------------------------------
    
    % Equilibrium model - no MDR, no HIV
    p0 = p; r0 = r;
    p0.pu = 0; r0.cs  = r.cs_1997; r0.cs2 = r.cs2_1997;
    p0.rfbeta_mdr = 0; r0.MDR_acqu = 0;  p0.MDR_rec = [0 0]; 
    p0.rfbeta_hiv = 0; r0.ART_init = 0;  
    M0 = make_smallmodel(p0, r0, i1, s1, gps1,prm);   % 60 compartments
   
    % Introduction of MDR, no ART (1970-1997),
    p1 = p; r1 = r;
    p1.pu = 0; r1.cs  = r.cs_1997; r1.cs2 = r.cs2_1997;
    p1.MDR_rec = [0 0];
    p1.rfbeta_hiv = 0; r1.ART_init = 0;               % p1.rfbeta_hiv not introduced
    M1a = make_smallmodel(p1, r1, i2, s2, gps2,prm);  % 115 compartments
    M1 = make_model2(p1, r1, i, s, gps,prm);          % 435 compartments (To transfer it into full matrix)
  
    % Start of RNTCP to start of ART  (prm.ART_start(here 2004 for India) (RNTCP scale up 1997 - 2007)
    p2 = p; r2 = r;
    p2.pu = 0 + (p.pu_2022 - 0)*((prm.ART_start-1997)/25);  % Assuming linear scaleup from 1997 to 2022
    r2.cs  = r.cs_1997 + (r.cs_2022 - r.cs_1997)*((prm.ART_start-1997)/25);
    r2.cs2 = r.cs2_1997 + (r.cs2_2022 - r.cs2_1997)*((prm.ART_start-1997)/25); % Assumed to be same from 1997 to 2022
    p2.MDR_rec = [0 0]; r2.ART_init = 0; %p.MDR_rec2015;         
    M2 = make_model2(p2, r2, i, s, gps,prm);          % 435 compartments 

    % Start of ART to start of DST in 2015 (scaleup 2003 - 2022)
    p3 = p; r3 = r;
    p3.pu = 0 + (p.pu_2022 - 0)*((2015-1997)/25);  % Assuming linear scaleup from 1997 to 2022
    r3.cs  = r.cs_1997 + (r.cs_2022 - r.cs_1997)*((2015-1997)/25);
    r3.cs2 = r.cs2_1997 + (r.cs2_2022 - r.cs2_1997)*((2015-1997)/25); % Assumed to be same from 1997 to 2022
    p3.MDR_rec = p.MDR_rec2015;   %p.MDR_rec2015       
    M3 = make_model2(p3, r3, i, s, gps,prm);          % 435 compartments 
     
    % Scale-up of DST (2015-2022)
    p4 = p; r4 = r;
    p4.pu = p.pu_2022; r4.cs  = r.cs_2022; r4.cs2 = r.cs2_2022; 
    p4.MDR_rec = p.MDR_rec2022;
    M4 = make_model2(p4, r4, i, s, gps,prm);         % 435 compartments 
%tic

    % --- Simulate the models
    % Equilibrium model from the simplified small model
    init1 = zeros(1,i1.nx); seed = 1e-6;
    init1(i1.U.a1.h0.v0) = prm.N(1)/sum(prm.N); init1(i1.U.a2.h0.v0) = prm.N(2)/sum(prm.N); init1(i1.U.a3.h0.v0) = prm.N(3)/sum(prm.N); 
    init1(i1.U.a4.h0.v0) = (prm.N(4)/sum(prm.N))-seed; init1(i1.U.a5.h0.v0) = prm.N(5)/sum(prm.N); init1(i1.I.a4.h0.ds) = seed;

    prm.growth = 0;
    geq = @(t,in) goveqs_basissmall(t, in, r, M0, i1, prm);
    [t0, soln0] = ode15s(geq, [0:5e3], init1, odeset('NonNegative',[1:i1.nstates]));
%keyboard;
    inds1 = intersect(intersect(s2.U, s2.h0),s2.v0);
    inds2 = intersect(intersect(intersect([s2.Lf, s2.Ls], s2.h0), s2.ds),s2.v0);
    inds3 = intersect(intersect([s2.Isc, s2.I, s2.E, s2.Rlo, s2.Rhi, s2.R], s2.h0), s2.ds);
    inds4 = intersect(intersect(intersect([s2.Dx,s2.Tx,s2.Tx2],s2.ds),s2.h0),s2.pr);
    inds0  = [inds1, inds2, inds3,inds4];

    % To allocate results from a small model to input vector for a bigger
    % model (60 to 115 compartments)
    init_smallmodel =  soln0(end,:);
    init2 = zeros(1,i2.nx);             %i.nx  115
    init2(inds0) = init_smallmodel;
    
     % Introduction of MDR  
    prm.growth = prm.gth(1);
    geq = @(t,in)goveqs_basissmall(t, in, r, M1a, i2, prm);  
    [t1, soln1] = ode15s(geq, [1970 1997], init2, odeset('NonNegative',[1:i2.nstates]));

    % To initialise with MDR components
    inds1 = intersect(intersect(s.U, s.h0),s.v0);
    inds2 = intersect(intersect([s.Lf, s.Ls], s.h0),s.v0);
    inds3 = intersect([s.Isc, s.I, s.E, s.Rlo, s.Rhi, s.R], s.h0);
    inds4 = intersect(intersect([s.Dx,s.Tx,s.Tx2],s.h0),s.pr);
    inds_full  = [inds1, inds2, inds3,inds4];
    
    %init1 = init(inds_full);
    % To allocate results from a small model to input vector for a bigger model
    init_smallmodel2 =  soln1(end,:);
    init = zeros(1,i.nx);              %444     
    init(inds_full) = init_smallmodel2;
    
%      % Introduction of RNTCP % 
%     %init = init_smallmodel2; %soln1(end,:);
%     prm.growth = prm.gth(2);
%     geq = @(t,in) goveqs_scaleup(t, in, r, M1, M2, [1997 prm.ART_start], i, prm, sel, agg);
%     [t2, soln2] = ode15s(geq, [1997:1:prm.ART_start], init, odeset('NonNegative',[1:i.nstates]));

      % Before Introduction of RNTCP % 
    prm.growth = 0; %prm.gth(2);  %2
    %geq = @(t,in) goveqs_scaleup(t, in, r, M1, M2, [1997 prm.ART_start], i, prm, sel, agg);
    geq = @(t,in) goveqs_basis2(t, in, r, M1, i, prm, sel, agg); %1e3
    [t2a, soln2a] = ode15s(geq, [0:600], init, odeset('NonNegative',[1:i.nstates]));
   
      % Introduction of RNTCP %
    prm.growth = prm.gth(2);  %2
    geq = @(t,in) goveqs_scaleup(t, in, r, M1, M2, [1997 prm.ART_start], i, prm, sel, agg);
    %geq = @(t,in) goveqs_basis2(t, in, r, M1, i, prm, sel, agg); 
    [t2, soln2] = ode15s(geq, [1997:1:prm.ART_start], soln2a(end,:), odeset('NonNegative',[1:i.nstates]));
  

    % Introduction of ART (prm.ART_start(here 2004) onwards 2003-2015)
    init = soln2(end,:); prm.growth = prm.gth(3);
    geq = @(t,in) goveqs_scaleup(t, in, r, M2, M3, [prm.ART_start 2015], i, prm, sel, agg);
    [t3, soln3] = ode15s(geq, [prm.ART_start:1:2015], init, odeset('NonNegative',[1:i.nstates]));

    % DST scale-up (Assuming from 2015 to 2022)
    init = soln3(end,:); prm.growth = prm.gth(4); 
    %[t4, soln4] = ode15s(@(t,in) goveqs_scaleup(t, in,r, M3, M4, [2015 2022], i, prm, sel, agg), [2015 2023], init, odeset('NonNegative',[1:i.nstates]));
    geq = @(t,in) goveqs_scaleup(t, in, r, M3, M4, [2015 2022], i, prm, sel, agg);
    [t4, soln4] = ode15s(geq, [2015:2023], init, odeset('NonNegative',[1:i.nstates]));

  keyboard;
    allt   = [t2; t3(2:end); t4(2:end)]; % ;   % 1997:2023
    allsol = [soln2; soln3(2:end,:);soln4(2:end,:)];

    soln = interp1(allt, allsol, 1997:2023);
    dsol   = diff(soln,[],1);
    sfin  = soln(end,:);
    
    pop = sum(soln(1:end-1,1:i.nstates),2);  %1997:2022 (1:end -> 1997:2023)
%     %%%%%%%%%%%%%%%%%% To estimate Lives saved  %%%%%%%%%%
%     % Initial condition in 2000
%     init2000 = allsol(4,:); %soln2(4,:);
%     % Initial condition in 2005
%     init2005 = allsol(9,:); %soln3(end-10,:);
% 
%     % Condition in 1999
%     p2 = p; r2 = r;
%     p2.pu = 0 + (p.pu_2022 - 0)*((1999-1997)/25);  % Assuming linear scaleup from 1997 to 2022
%     r2.cs  = r.cs_1997 + (r.cs_2022 - r.cs_1997)*((1999-1997)/25);
%     r2.cs2 = r.cs2_1997 + (r.cs2_2022 - r.cs2_1997)*((1999-1997)/25); % Assumed to be same from 1997 to 2022
%     r2.ART_init = 0; p2.MDR_rec = [0 0];         
%     M2a = make_model2(p2, r2, i, s, gps,prm);
% 
%     % Scenario1: In absense of program (removal of program in 2000)
%     p5 = p; r5 = r;
%     p5.pu = 0;  
%     r5.cs  = r.cs_1997;%+ (r.cs_2022 - r.cs_1997)*((2000-1997)/25);
%     r5.cs2 = r.cs2_1997;% + (r.cs2_2022 - r.cs2_1997)*((2000-1997)/25); % Assumed to be same from 1997 to 2022
%     r5.ART_init = 0; p5.MDR_rec = [0 0];          
%     MN2000 = make_model2(p5, r5, i, s, gps,prm); % Null coverage
% 
%   
%   % Fixing the intervention parameters in 2000
%     p6 = p; r6 = r;
%     p6.pu = 0 + (p.pu_2022 - 0)*((2000-1997)/25);  
%     r6.cs  = r.cs_1997 + (r.cs_2022 - r.cs_1997)*((2000-1997)/25);
%     r6.cs2 = r.cs2_1997 + (r.cs2_2022 - r.cs2_1997)*((2000-1997)/25); % Assumed to be same from 1997 to 2022
%     r6.ART_init = 0; p6.MDR_rec = [0 0];          
%     MC2000 = make_model2(p6, r6, i, s, gps,prm); % constant coverage
% 
%      % Condition in 2004
%     p2 = p; r2 = r;
%     p2.pu = 0 + (p.pu_2022 - 0)*((2004-1997)/25);  % Assuming linear scaleup from 1997 to 2022
%     r2.cs  = r.cs_1997 + (r.cs_2022 - r.cs_1997)*((2004-1997)/25);
%     r2.cs2 = r.cs2_1997 + (r.cs2_2022 - r.cs2_1997)*((2004-1997)/25); % Assumed to be same from 1997 to 2022
%     p2.MDR_rec = [0 0]; %r2.ART_init = 0;        
%     M2b = make_model2(p2, r2, i, s, gps,prm);
% 
%    % In absense of program (removal of program in 2005 and back to 1997 condition)
%     p7 = p; r7 = r;
%     p7.pu = 0;  
%     r7.cs  = r.cs_1997;% + (r.cs_2022 - r.cs_1997)*((2005-1997)/25);
%     r7.cs2 = r.cs2_1997;% + (r.cs2_2022 - r.cs2_1997)*((2005-1997)/25); % Assumed to be same from 1997 to 2022
%     p7.MDR_rec = [0 0]; r7.ART_init = 0;           
%     MN2005 = make_model2(p7, r7, i, s, gps,prm); % Null coverage
% 
%     % Fixing the intervention parameters in 2005
%     p8 = p; r8 = r;
%     p8.pu = 0 + (p.pu_2022 - 0)*((2005-1997)/25);  % Assuming linear scaleup from 1997 to 2022
%     r8.cs  = r.cs_1997 + (r.cs_2022 - r.cs_1997)*((2005-1997)/25);
%     r8.cs2 = r.cs2_1997 + (r.cs2_2022 - r.cs2_1997)*((2005-1997)/25); % Assumed to be same from 1997 to 2022
%     p8.MDR_rec = [0 0]; %r8.ART_init = 0;         
%     MC2005 = make_model2(p8, r8, i, s, gps,prm); % constant coverage
% 
%     % Without program from 2000
%     prm.growth = prm.gth(2);
%    %geq = @(t,in)goveqs_basis2(t, in, r5, MN2000, i, prm, sel, agg);  
%     geq = @(t,in) goveqs_scaleup(t, in, r5, M2a, MN2000, [2000 2023], i, prm, sel, agg);
%    [tb, solnb] = ode15s(geq, [2000:1:2023], init2000, odeset('NonNegative',[1:i.nstates]));
% 
%    % With constant coverage as in 2000
%     prm.growth = prm.gth(2);
%     %geq = @(t,in)goveqs_basis2(t, in, r6, MC2000, i, prm, sel, agg); 
%     geq = @(t,in) goveqs_scaleup(t, in, r6, M2a, MC2000, [2000 2023], i, prm, sel, agg);
%     [tc, solnc] = ode15s(geq, [2000:1:2023], init2000, odeset('NonNegative',[1:i.nstates]));
%         
%      % Without program from 2005
%      prm.growth = prm.gth(3);
%      %geq = @(t,in)goveqs_basis2(t, in, r7, MN2005, i, prm, sel, agg);  
%      geq = @(t,in) goveqs_scaleup(t, in, r7, M2b, MN2005, [2005 2023], i, prm, sel, agg); % Arbitary scale up
%     [td, solnd] = ode15s(geq, [2005:1:2023], init2005, odeset('NonNegative',[1:i.nstates]));
%      
%     % Baseline with constant coverage as in 2005
%     prm.growth = prm.gth(3);
%     %geq = @(t,in)goveqs_basis2(t, in, r8, MC2005, i, prm, sel, agg);  
%     geq = @(t,in) goveqs_scaleup(t, in, r8, M2b, MC2005, [2005 2023], i, prm, sel, agg); % Arbitary scale up
%     [te, solne] = ode15s(geq, [2005:1:2023], init2005, odeset('NonNegative',[1:i.nstates]));
%
%     soln_b = interp1(tb, solnb, 2000:2023);  % With null coverage from 2000
%     dsolb  = diff(soln_b,[],1);
% 
%     soln_c = interp1(tc, solnc, 2000:2023); % With constant coverage from 2000
%     dsolc  = diff(soln_c,[],1);
% 
%     soln_d = interp1(td, solnd, 2005:2023); % With null coverage from 2005
%     dsold  = diff(soln_d,[],1);
% 
%     soln_e = interp1(te, solne, 2005:2023); % With constant coverage from 2005
%     dsole  = diff(soln_e,[],1);
%      
%       incd   = dsol(:,i.aux.inc(1))*1e5;  % 1997-2022
%       incd_b = dsolb(:,i.aux.inc(1))*1e5; % 2000-2022
%       incd_c = dsolc(:,i.aux.inc(1))*1e5; % 2000-2022
%       incd_d = dsold(:,i.aux.inc(1))*1e5; % 2005-2022
%       incd_e = dsole(:,i.aux.inc(1))*1e5; % 2005-2022
%       mort   = dsol(:,i.aux.mort(1))*1e5; % 1997-2022
%       mort_b = dsolb(:,i.aux.mort(1))*1e5;% 2000-2022
%       mort_c = dsolc(:,i.aux.mort(1))*1e5;% 2000-2022
%       mort_d = dsold(:,i.aux.mort(1))*1e5;% 2005-2022
%       mort_e = dsole(:,i.aux.mort(1))*1e5;% 2005-2022

%  toc
    % --- Get the objectives ----------------------------------------------
    inc2000      = (sum(dsol(end-22,i.aux.inc(1)))./pop(end-22))*1e5;
    inc2022      = (sum(dsol(end,i.aux.inc(1)))./pop(end))*1e5;
    inc_h1       = (dsol(end,i.aux.inc(2))./pop(end))*1e5;     %([2,3])
    mdr2015      = (dsol(end-7,i.aux.inc(3))./pop(end-7))*1e5; %end--> 2022
    mdr2022      = (dsol(end,i.aux.inc(3))./pop(end))*1e5;
    noti_cu      =  sum(sum((dsol(end-22:end,i.aux.noti)./pop(end-22:end))*1e5,2));  %cumulative notification from 2000 to 2022
    noti2022     = (sum(dsol(end,i.aux.noti))./pop(end))*1e5; 
    ART_covg     = sum(sfin(s.hart)/sum(sfin([s.h1,s.hart])));
    HIV_prev     = sum(sfin([s.h1, s.hart]))/sum(sfin(1:i.nstates));
    mort_H02000  = (dsol(end-22,i.aux.mort(1))./pop(end-22))*1e5;
    mort_H02022  = (dsol(end,i.aux.mort(1))./pop(end))*1e5;
    mort_H12022  = (dsol(end,i.aux.mort(2))./pop(end))*1e5;
    mdriniTX     = (dsol(end,i.aux.noti2)./pop(end))*1e5;  % MDR treatment initiation 2022
    sym_sim      =  sum(sfin(s.I))/sum(sfin(s.prevalent)); 
    inc2015      = (sum(dsol(end-7,i.aux.inc(1)))./pop(end-7))*1e5;
    %inc2003      = (sum(dsol(end-19,i.aux.inc(1)))./pop(end-19))*1e5;
    
    
%      % --- Get the objectives ----------------------------------------------
%     ainc2000      = sum(dsol(end-22,i.aux.inc(1)))*1e5;
%     ainc2022      = sum(dsol(end,i.aux.inc(1)))*1e5;
%     ainc_h1       = dsol(end,i.aux.inc(2))*1e5;     %([2,3])
%     amdr2015      = dsol(end-7,i.aux.inc(3))*1e5; %end--> 2022
%     amdr2022      = dsol(end,i.aux.inc(3))*1e5;
%     anoti_cu      =  sum(sum(dsol(end-22:end,i.aux.noti)*1e5,2));  %cumulative notification from 2000 to 2022
%     anoti2022     = sum(dsol(end,i.aux.noti))*1e5; 
%     aART_covg     = sum(sfin(s.hart)/sum(sfin([s.h1,s.hart])));
%     aHIV_prev     = sum(sfin([s.h1, s.hart]))/sum(sfin(1:i.nstates));
%     amort_H02000  = dsol(end-22,i.aux.mort(1))*1e5;
%     amort_H02022  = dsol(end,i.aux.mort(1))*1e5;
%     amort_H12022  = dsol(end,i.aux.mort(2))*1e5;
%     amdriniTX     = dsol(end,i.aux.noti2)*1e5;  % MDR treatment initiation 2019
%     asym_sim      =  sum(sfin(s.I))/sum(sfin(s.prevalent)); 
%     ainc2015      = sum(dsol(end-7,i.aux.inc(1)))*1e5;
%     %inc2003      = sum(dsol(end-19,i.aux.inc(1)))*1e5;

  
    if inc2022 < 1
        out = Inf*calfn.sgn;
        aux = [];
    else
        %out = calfn.fn(inc2003, inc2022, inc_h1, noti, ART_covg, HIV_prev, mort_H0, mort_H1, mdr2015, mdr2022, mdriniTX, sym_sim);

       out = calfn.fn(sym_sim,inc2000, inc2022, mort_H02000,mort_H02022,noti_cu, noti2022, mdr2015, mdr2022, mdriniTX,HIV_prev,inc_h1,ART_covg,mort_H12022);

        
        if imag(out)>0
            out = Inf*calfn.sgn;
            aux = [];
            
        else
            
            % --- Get additional outputs and package --------------------------
            aux.soln     = soln; %[soln, t];         % 2003:2023
            aux.allsol   = [soln2; soln3(2:end,:); soln4(2:end,:)];  % 1997:2023
            aux.inc2022  = inc2022;
            aux.inc2015 = inc2015;
            aux.inc2000 = inc2000;
            aux.inc_h1   = inc_h1;
            aux.noti     = noti2022;
            aux.noti_cu  = noti_cu;
            aux.ART_covg = ART_covg;
            aux.HIV_prev = HIV_prev;
            aux.mort_H02000  = mort_H02000;
            aux.mort_H02022  = mort_H02022;
            aux.mort_H1  = mort_H12022;
            aux.mdr2015  = mdr2015;
            aux.mdr2022  = mdr2022;
            aux.mdriniTX = mdriniTX;
            aux.sym      = sym_sim;
            aux.M1       = M1;
            
%             aux.incd  = incd;
%             aux.incdb = incd_b;
%             aux.incdc = incd_c;
%             aux.incdd = incd_d;
%             aux.incde = incd_e;
%             aux.mort  = mort;
%             aux.mortb = mort_b;
%             aux.mortc = mort_c;
%             aux.mortd = mort_d;
%             aux.morte = mort_e;

        end
        
    end
end
