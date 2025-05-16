clear all; load Model_setup;

tic
obj = @(x) get_objective(x, prm, ref, sel, agg, gps,gps1,gps2,lhd);
         
         
% --- Generate initial guesses

nsam = 1e3;
%nx   = xi.calib;
nx   = xi.nx;

lo = repmat(prm.bounds(1,1:nx),nsam,1);
hi = repmat(prm.bounds(2,1:nx),nsam,1);
sams = lo + (hi-lo).*lhsdesign(nsam,nx);  
%load sams;     % 

% out = zeros(1,nsam);
% mk  = round(nsam/25); %25
% for ii = 1:nsam
%     if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
%     out(ii) = obj(sams(ii,:));
% end
% fprintf('\n');
% 
% % Order them to get the best-fitting ones
% % out = out(~isnan(out));
% % mat = sortrows([out;1:length(out)]',-1);
% 
% tmp1 = [out;1:nsam]';
% tmp2 = tmp1(~isnan(out),:);
% mat  = sortrows(tmp2,-1);
% 
% x0   = sams(mat([1:3],2),:);
% 
% % First round of MCMC
% [xsto1, outsto1] = MCMC_adaptive2(obj, x0(1,:), 5e4, 1, [], [], [], true); 
                 
load xsto1
%xsto1(:,22)=1; %xsto1(:,23)=1;
load outsto1

% Second round of MCMC
[rr,cc] = find(outsto1==max(outsto1(:)));
xx = xsto1(cc(1),:,rr(1));
mat  = xsto1(3e4:end,:,1);
cov0 = cov(mat);
%cov0(22,:) = 0; cov0(:,22) = 0; cov0(22,22) = 1e-3;
[xsto2, outsto2] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 


%Next round of MCMC
[rr,cc] = find(outsto2==max(outsto2(:)));
xx = xsto2(cc(1),:,rr(1));
mat  = xsto2(3e4:end,:,1);
cov0 = cov(mat);
[xsto3, outsto3] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 

% % Next round of MCMC
% [rr,cc] = find(outsto3==max(outsto3(:)));
% xx = xsto3(cc(1),:,rr(1));
% mat  = xsto3(3e4:end,:,1);
% cov0 = cov(mat);
% [xsto4, outsto4] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 

% % Next round of MCMC
% [rr,cc] = find(outsto4==max(outsto4(:)));
% xx = xsto4(cc(1),:,rr(1));
% mat  = xsto4(3e4:end,:,1);
% cov0 = cov(mat);
% [xsto5, outsto5] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 

% % Next round of MCMC
% [rr,cc] = find(outsto5==max(outsto5(:)));
% xx = xsto5(cc(1),:,rr(1));
% 
% mat  = xsto5(3e4:end,:,1);
% cov0 = cov(mat);
% [xsto6, outsto6] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 

% % Second round of MCMC
% [rr,cc] = find(outsto6==max(outsto6(:)));
% xx = xsto6(cc(1),:,rr(1));
% 
% mat  = xsto6(3e4:end,:,1);
% cov0 = cov(mat);
% [xsto7, outsto7] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 
% 
% % Second round of MCMC
% [rr,cc] = find(outsto7==max(outsto7(:)));
% xx = xsto7(cc(1),:,rr(1));
% 
% mat  = xsto7(3e4:end,:,1);
% cov0 = cov(mat);
% [xsto8, outsto8] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 
% 
% % Second round of MCMC
% [rr,cc] = find(outsto8==max(outsto8(:)));
% xx = xsto8(cc(1),:,rr(1));
% 
% mat  = xsto8(3e4:end,:,1);
% cov0 = cov(mat);
% [xsto9, outsto9] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 
% 
% % Second round of MCMC
% [rr,cc] = find(outsto9==max(outsto9(:)));
% xx = xsto9(cc(1),:,rr(1));
% 
% mat  = xsto9(3e4:end,:,1);
% cov0 = cov(mat);
% [xsto10, outsto10] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 


% % Second round of MCMC
% [rr,cc] = find(outsto2==max(outsto2(:)));
% xx = xsto2(cc(1),:,rr(1));
% xx(1,8) = 1;     % Introduced 8th column for mortality factor
% 
% mat  = xsto2(3e4:end,:,1);
% cov0 = zeros(8);
% cov0(1:7,1:7) = cov(mat);
% cov0(8,8) = 1;
% 
% [xsto3, outsto3] = MCMC_adaptive2(obj, xx, 5e4, 1, [], [], cov0, true); 
 
save calibres2.mat;  % With population correction


% Alarm me once the program runnning has been finished
load handel  
sound(y,Fs)


