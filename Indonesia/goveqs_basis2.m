function [out, lam] = goveqs_basis2(t, in, r, M, i, prm, sel, agg)         % <--- Arguments changed to accommodate prm

invec = in(1:i.nstates);

if t<1980                                                                  % <--- NEW
    rHIV = 0;
else
    rHIV = interp1(0:length(prm.rHIV)-1, prm.rHIV, (t-1980)); 
end

fac = (r.betared)^max(t-2000, 0); % trend on activation and progression rate

try
% Normalise by populations
lam = (M.lambda*invec)/sum(invec); 
nlin_ds = (lam(1)*M.nlin.a1.ds+ lam(2)*M.nlin.a2.ds+ lam(3)*M.nlin.a3.ds+lam(4)*M.nlin.a4.ds+ lam(5)*M.nlin.a5.ds); 
nlin_mdr = lam(6)*M.nlin.a1.mdr+ lam(7)*M.nlin.a2.mdr+ lam(8)*M.nlin.a3.mdr +lam(9)*M.nlin.a4.mdr+ lam(10)*M.nlin.a5.mdr; 
allmat = M.lin + M.Dxlin + rHIV*M.linHIV + fac*M.linLTBI + nlin_ds + nlin_mdr;


out = allmat*invec;
catch
   keyboard; 
end

% Implement deaths
morts = sum(M.mortvec,2).*invec;
out = out - morts;
% Implement births
births = sum(morts) + prm.growth;
out(i.U.a1.h0.v0) = out(i.U.a1.h0.v0)+births;

% Get the auxiliaries
out(i.aux.inc)     = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.noti)    = agg.noti*(sel.noti.*allmat)*invec;
out(i.aux.noti2)   = sum((sel.noti2.*allmat)*invec);
out(i.aux.mort(1)) = sum(M.mortvec(:,2).*invec);
out(i.aux.mort(2)) = sum(M.mortvec(:,3).*invec);

