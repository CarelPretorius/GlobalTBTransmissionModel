function [r,p] = alloc_parameters(x,r,p,xi)

r.beta         = x(xi.r_beta);                                              % Transmission for HIV-ve, DS-TB
p.rfbeta_mdr   = x(xi.rfbeta_mdr);
p.rfbeta_hiv   = x(xi.rfbeta_hiv);
r.betared      = x(xi.r_betared);
r.sym          = x(xi.r_sym);  
p.pu_2022      = x(xi.p_pu);
r.cs_1997      = x(xi.r_cs).*[1,1,1,1,1];
r.cs_2022      = x(xi.r_cs)*x(xi.rf_cs2022).*[1,1,1,1,1];
r.cs2_1997     = x(xi.r_cs)*x(xi.rf_cs2).*[1,1,1,1,1];
r.cs2_2022     = x(xi.r_cs)*x(xi.rf_cs2).*[1,1,1,1,1];
r.cs3          = r.cs3.*[1,1,1,1,1];
r.mort_TB      = x(xi.r_mort_TB);
p.Dx           = x(xi.p_Dx);
r.default      = r.Tx*(1-x(xi.p_Tx_complete))./x(xi.p_Tx_complete);
r.ART_init     = x(xi.r_ART_init);
p.MDR_rec2022  = [x(xi.p_MDRrec2022),0];
r.MDR_acqu     = x(xi.r_MDR_acqu);
r.HIV_mort     = x(xi.r_HIV_mort);

r.progression([2,3])  = r.progression(1)*x(xi.p_HIV_relrate)*[1 0.4];
r.reactivation([2,3]) = r.reactivation(1)*x(xi.p_HIV_relrate)*[1 0.4];
r.progression   = r.progression.*[1 1 1 1 1]';    %[0.0826    3.4262    1.3705].*[1 1 1 1 1]'; %
r.reactivation  = r.reactivation.*[1 1 1 1 1]';   %[0.0006    0.0249    0.0100].*[1 1 1 1 1]'; %
%keyboard;

%if length(x) > xi.calib
   r.self_cure([1,3]) = x(xi.r_self_cure);
 %  r.LTBI_cure    = r.LTBI_cure*x(xi.rf_LTBI_cure);
%end