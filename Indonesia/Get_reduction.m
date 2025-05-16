%%% Define a function to estimate the percentage reduction from two distributions 
function  pct_reduction = Get_reduction(Nold, Nnew)

% Estimate standard deviations assuming symmetric normal distribution
% Using 95% CI: upper ≈ mean + 1.96*std
stdold = (Nold(3) - Nold(1)) / 1.96;
stdnew = (Nnew(3) - Nnew(1)) / 1.96;

% Number of Monte Carlo samples
N = 1e5;
% Sample from the normal distributions
samplesold = normrnd(Nold(1), stdold, [N, 1]);
samplesnew = normrnd(Nnew(1), stdnew, [N, 1]);

% Compute percentage reduction
perc_reduction = 100 * (samplesold - samplesnew) ./ samplesold;
pct_reduction = prctile(perc_reduction,[2.5,50,97.5]);
end
