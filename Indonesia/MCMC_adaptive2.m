% Adaptive MCMC, using Haario et al:
% https://link.springer.com/article/10.1007/s11222-008-9110-y
% and
% http://probability.ca/jeff/ftpdir/adaptex.pdf

% F:          Function giving log-posterior density for a parameter set x
% x0:         Initial value of parameter set x
% n:          Number of iterations
% cov0:       Initial covariance matrix
% fac:        Scaling factor for covariance matrix. Set fac = 1 for default
% fixinds:    Elements of x that should be held fixed. Set fixinds = [] for full MCMC
% blockinds:  Number of 'epi parameters' (e.g. beta, X2) in x, if we want to vary epi and non-epi parameters as independent 'blocks'. Set blockinds = [] if we want a full covariance matrix
% displ:      Structure with display options. Set displ = true to show progress

function [xsto, outsto, history, accept_rate] = MCMC_adaptive(F, x0, n, sigma, fixinds, blockind, cov0, displ)

d = length(x0); b = 0.05; sd = sigma*2.4^2/d;
if ~isempty(fixinds)
    inds = fixinds(1,:); %vals = fixinds(2,:);
    vals = x0(fixinds);
else
    inds = []; vals = [];
end

if isempty(cov0)
    % Checks on the initial covariance matrix
    cov0 = eye(d)*1e-3; cov0(inds,:) = 0; cov0(:,inds) = 0;
    cov0(1:blockind,blockind+1:end) = 0; cov0(blockind+1:end,1:blockind) = 0;
end
% Initiate the output matrices
xsto = zeros(d,n); outsto = zeros(1,n);
history = zeros(d+1,n);                                                    % Rows: 1:d Proposed values 4. Accept or reject

xsto(:,1) = x0(:); xbar = xsto;
%keyboard;
FX = F(x0); outsto(1) = FX;
acc = 0;

if displ; figure; end

% --- Start the MCMC loop -------------------------------------------------
for t = 2:n
    X = xsto(:,t-1);
    
    % --- Make a proposal from the distribution
    Y0 = mvnrnd(X,0.1^2*cov0*sigma/d);
    if t < 2*d
        Y = Y0;
        Y(inds) = vals;
        
        mu     = mean(xsto(:,1:t-1),2);
        covmat = cov(xsto(:,1:t-1)');
    else        
        mat = sd*covmat; mat = (mat+mat.')/2;
        try
            Y = max((1-b)*mvnrnd(X,mat) + b*Y0,0);
        catch
            keyboard;
        end
        Y(inds) = vals;
    end
    
    
    history(1:d,t) = Y;
    
    % --- Decide whether to accept or not
    FY = F(Y);
    if rand < exp(FY-FX)
        % Accept
        xsel = Y(:);
        FX = FY;
        acc = acc+1;
        history(end,t) = 1;
    else
        % Reject
        xsel = xsto(:,t-1);
    end
    xsto(:,t) = xsel;
    outsto(t) = FX;
    xbar(:,t) = (xbar(:,t-1)*(t-1) + xsel)/t;
    
    % Update means and covariance matrices
    if t>=2*d
        mu     = (t-1)/t*mu + xsel/t;
        covmat = (t-1)/t*covmat + 1/(t-1)*(xsel-mu)*(xsel-mu)';
    end
    
    % Display options
    if displ && (mod(t,round(n/25))==0); fprintf('%0.5g ', t/n*25); end
    if displ && (mod(t,200)==0)
        plot(xsto(setdiff(1:d,fixinds),1:t-1)'); xlim([0 n]); title(sprintf('%0.5g',acc/t));
        drawnow;
    end
    
end

fprintf('\n');
accept_rate = acc/n;
xsto = xsto';
history = history';