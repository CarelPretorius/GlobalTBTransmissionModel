function [epi_outputs]=RunGlobalEpiModelb(PARS,IntvnScn, PathToModel, PathToScenarioModel)
cd(PathToModel)

epi_outputs=[];
%clear all; %load calibres2; xsto = xsto4;
%load calibres1; 

xsto2=PARS{1};
outsto2=PARS{2};

xsto = xsto2;
[M,I] = max(outsto2);

xi=PARS{3};
obj=PARS{4};

i0=PARS{5};
s0=PARS{6};
gps0=PARS{7};
sel0=PARS{8};
agg0=PARS{9};
r=PARS{10};
p=PARS{11};
prm=PARS{12};


%load Model_setup;

xs = xsto(I,:);  % Most probable parameter set 

% ix0 = 1e4; nx = 2; dx = round((size(xsto,1)-ix0)/nx); %2e4 nx = 150
% xs_all = xsto(ix0:dx:end,:,1);
% 
% pct_xs= prctile(xs_all,[2.5,50,97.5])


% % --- Get the Regional notifications
tend1   = IntvnScn.year1;%2024;  % Intervention starts except vaccination  
tend1b  = IntvnScn.year2;%2028;  % Vaccination starts
tend2   = IntvnScn.year3;%2036;  % End date for simulation

% --- Do the simulations 
opts = odeset('NonNegative',[1:i0.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);


r0 = r; p0 = p;

inct   = [];
inch1  = [];
mdr    = [];
morth0 = [];
morth1 = [];
noti   = [];
notimdr= [];


[r,p] = alloc_parameters(xs,prm.r,prm.p,xi);
% Get the initial conditions
[out, aux] = obj(xs);
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

 % --- Now simulate without intervention (baseline)
M0 = make_model3(p, r, i0, s0, gps0, prm); 
prm.growth = prm.gth(3);

% Find the initial condition at the starting point for the interventions
geq = @(t,in) goveqs_scaleup(t, in, r, M0, M0, [2022 2023], i0, prm, sel0, agg0);
[ta, solna] = ode15s(geq, [2022 tend1], init0, odeset('NonNegative',[1:i0.nstates]));
initb = solna(end,:);

% --- Next simulate with interventions
%r1 = r; p1 = p;
%Set intervention to model parameters
[p1,r1]=SetInterventionScenario(p,r,IntvnScn);

M1 = make_model3(p1, r1, i0, s0, gps0, prm);


% p1
% r1
% pause

geq = @(t,in) goveqs_scaleup(t, in, r1, M0, M1, [tend1 tend1b], i0, prm, sel0, agg0);
[tb, solnb] = ode15s(geq, [tend1 tend2], initb, opts);

soln  = [solna; solnb(2:end,:)];
t     = [ta;    tb(2:end)];
soln8 = interp1(t, soln, tref); 

allmat = soln8;
dallmat = diff(allmat,[],1);

inct(:,:) = squeeze(dallmat(:,i0.aux.inc(1),:))*1e5;    % Total incidence
inch1(:,:) = squeeze(dallmat(:,i0.aux.inc(2),:))*1e5;   % HIV incidence
mdr(:,:) = squeeze(dallmat(:,i0.aux.inc(3),:))*1e5;     % MDR incidence
morth0(:,:) = squeeze(dallmat(:,i0.aux.mort(1),:))*1e5; % hiv-negative mortality
morth1(:,:) = squeeze(dallmat(:,i0.aux.mort(2),:))*1e5; % hiv-positive mortality
noti(:,:) = squeeze(sum(dallmat(:,i0.aux.noti,:),2))*1e5; % Total notification
notimdr(:,:) = squeeze(dallmat(:,i0.aux.noti2,:))*1e5;  % MDR Treatment initiation


epi_outputs{1}=inct;
epi_outputs{2}=inch1;
epi_outputs{3}=mdr;
epi_outputs{4}=morth0;
epi_outputs{5}=morth1;
epi_outputs{6}=noti;
epi_outputs{7}=notimdr;

%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%




    
cd(PathToScenarioModel) 

end
 