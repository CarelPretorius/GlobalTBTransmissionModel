function out = goveqs_scaleup(t, in,r, M0, M1, times, i, prm, sel, agg)

scale = max((t-times(1))/(times(2)-times(1)),0);                           % <--- NEW: removed min(,1), so that model can continue extrapolating ART coverage into future
Mt = M1; Mt.lin = M0.lin + scale*(M1.lin - M0.lin); 
Mt.Dxlin =  M0.Dxlin + scale*(M1.Dxlin - M0.Dxlin);
Mt.linHIV = M0.linHIV + scale*(M1.linHIV - M0.linHIV);
Mt.linLTBI = M0.linLTBI + scale*(M1.linLTBI - M0.linLTBI);
out = goveqs_basis2(t, in, r, Mt, i, prm, sel, agg);
