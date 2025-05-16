function M = make_smallmodel(p, r, i, s, gps,prm)  % Kenya

% --- Get the linear rates ------------------------------------------------
m  = zeros(i.nstates);
m2 = zeros(i.nstates);
m4 = zeros(i.nstates);

for ia = 1:length(gps.age)
    age  = gps.age{ia};

for istr = 1:length(gps.strains)
    strain = gps.strains{istr};
    
    for ih = 1:length(gps.hiv)
        hiv = gps.hiv{ih};
        
        getst = @(st) i.(st).(age).(hiv).(strain);
         
        Isc  = getst('Isc');     % Subclinical 
        I    = getst('I');       % clinical symptom 
        E    = getst('E');
        Dxs  = getst('Dx');
        Txs  = getst('Tx');
        Tx2s = getst('Tx2');
        Rlo  = getst('Rlo');
        Rhi  = getst('Rhi');
        R    = getst('R');
       
        % Group Selectors
        ismdr   = strcmp(strain, 'mdr');
        
        for iv = 1:length(gps.vacc)
            vacc = gps.vacc{iv};
            getst2 = @(st) i.(st).(age).(hiv).(vacc).(strain);
            Lf = getst2('Lf');
            Ls = getst2('Ls');
            
            % --- Fast progression and LTBI stabilisation
            source  = Lf;
            destins = [Isc,                                                   Ls];
            rates   = [r.progression(ia,ih).*p.PT_PLHIV(ih)*(1-(iv==2)*p.VE(2)), r.LTBI_stabil(ih)];     % p.PT_PLHIV = [1 1 (1-0.6)] during PT among PLHIV
            m4(destins, source) = m4(destins, source) + rates';
           
            % --- Reactivation (ia,ih)
            source = Ls; destin = Isc; rate = r.reactivation(ia,ih).*p.PT_PLHIV(ih)*(1-(iv==2)*p.VE(2));  % p.PT_PLHIV = [1 1 (1-0.6)] during PT among PLHIV
            m4(destin, source) = m4(destin, source) + rate;
            
        end
        
        % --- Sub-clinical to clinical Infection
        source = Isc;
        destin = I;
        rate   = r.sym;                                 % rate of symptom development (estimate this value)
        m(destin, source) = m(destin, source) + rate;

        
        % Introduce public/private sector
        % --- Primary careseeking, including access to public sector care
%         source = I; destin = Dxs.pu; rate = r.cs*p.pu + r.access;
%         m2(destin, source) = m2(destin, source) + rate;
        
        source = I; destin = Dxs.pr; rate = r.cs(ia)*(1-p.pu);
        m2(destin, source) = m2(destin, source) + rate;
        
        
%         % --- Secondary careseeking
%         source = E; destin = Dxs.pu; rate = r.cs2*p.pu;
%         m2(destin, source) = m2(destin, source) + rate;
        
        source = E; destin = Dxs.pr; rate = r.cs2(ia)*(1-p.pu);
        m2(destin, source) = m2(destin, source) + rate;
        
        
%         % --- Case finding from subclinical TB
%         source = Isc; destin = Dxs.pu; rate = r.cs3*p.pu ;
%         m2(destin, source) = m2(destin, source) + rate;
        
        source = Isc; destin = Dxs.pr; rate = r.cs3(ia)*(1-p.pu);
        m2(destin, source) = m2(destin, source) + rate;
  
        
        for ip = 1:length(gps.provs)
            prov = gps.provs{ip};
            
            Dx = Dxs.(prov); Tx = Txs.(prov); Tx2 = Tx2s.(prov);
            
            try
                pFLinit = p.Dx(ip)*p.Tx_init(ip)*(1 - ismdr*p.MDR_rec(ip)); %*p.Tx_init2(ip)
            catch
                keyboard
            end
            pSLinit = p.Dx(ip)*p.Tx_init2(ip)*ismdr*p.MDR_rec(ip);
            p_ltfu  = 1-(pFLinit+pSLinit); %p.Dx(ip)*p.Tx_init(ip);
            
            source  = Dx;
            destins =          [Tx,       Tx2,      E];
            rates   = r.Dx*[pFLinit,  pSLinit,  p_ltfu];
            m(destins, source) = m(destins, source) + rates';
            
            % --- FL Treatment
            pFLcure = p.cure(ip)*(1-ismdr);
            pSLtran = p.SL_trans(ip)*ismdr;
            rMDRacq = r.MDR_acqu*(1-ismdr);
            
            source  = Tx;
            destins = [Rlo            Tx2                        E                              Rhi]; %,             i.Tx.(age).(hiv).mdr.(prov)];
            rates   = [r.Tx*pFLcure,  r.Tx*(1-pFLcure)*pSLtran,  r.Tx*(1-pFLcure)*(1-pSLtran),  r.default(ip)];%,   rMDRacq];
            m(destins, source) = m(destins, source) + rates';
            
            % --- SL Treatment
            source  = Tx2;
            destins = [Rlo                 E];          %Rhi];
            rates   = [r.Tx2*p.cure2(ip),  r.Tx2*(1-p.cure2(ip))+ r.default2(ip)];
            m(destins, source) = m(destins, source) + rates';
            
            
        end  % provider loop
        
%         % --- LTBI cure
%         sources = [s.Lf, s.Ls];
%         destin  = i.U; %i.U.(age).(hiv).(vacc); %
%         rates   = r.LTBI_cure;
%         m(destin, sources) = m(destin, sources) + rates;

        % --- Relapse
        sources = [Rlo, Rhi, R];
        destin  = Isc;
        rates   = r.relapse;
        m(destin, sources) = m(destin, sources) + rates;
        
        sources = [Rlo, Rhi];
        destin  = R;
        rates   = 0.5;
        m(destin, sources) = m(destin, sources) + rates;
        
        % --- Self cure
        sources = intersect(intersect(intersect(s.infectious,s.(hiv)), s.(strain)),s.(age));
        destin  = Rhi;
        rates   = r.self_cure(ih);
        m(destin, sources) = m(destin, sources) + rates;
        
    end  % hiv loop
end % MDR loop



%   for ih = 1:length(gps.hiv)
%     hiv = gps.hiv{ih};
%     
%     source = i.U.(age).(hiv).v0;
%     destin = i.U.(age).(hiv).v1;
%     rate   = r.vacc;
%     m(destin,source) = m(destin,source) + rate;
%     
%     source = i.U.(age).(hiv).v1;
%     destin = i.U.(age).(hiv).v2;
%     rate   = r.waning;
%     m(destin,source) = m(destin,source) + rate;
%     
%     for istr = 1:length(gps.strains)
%         strain = gps.strains{istr};
%         
%         source = i.Lf.(age).(hiv).v0.(strain);
%         destin = i.Lf.(age).(hiv).v1.(strain);
%         rate   = r.vacc;
%         m(destin,source) = m(destin,source) + rate;
%         
%         source = i.Lf.(age).(hiv).v1.(strain);
%         destin = i.Lf.(age).(hiv).v2.(strain);
%         rate   = r.waning;
%         m(destin,source) = m(destin,source) + rate;
%         
%     end
%   end
end

% --- Introduction of Ageing rate                                           
sources = s.a1;
destins = s.a2;
inds    = sub2ind(size(m), destins, sources);
rate   = r.ageing(1);
m(inds) = m(inds) + rate;

sources = s.a2;
destins = s.a3;
inds    = sub2ind(size(m), destins, sources);
rate   = r.ageing(2);
m(inds) = m(inds) + rate;

sources = s.a3;
destins = s.a4;
inds    = sub2ind(size(m), destins, sources);
rate   = r.ageing(3);
m(inds) = m(inds) + rate;

sources = s.a4;
destins = s.a5;
inds    = sub2ind(size(m), destins, sources);
rate   = r.ageing(4);
m(inds) = m(inds) + rate;

% --- HIV acquisition                                                      % <--- NEW
m3 = zeros(i.nstates);
% sources = s.h0;
% destins = s.h1;
% inds    = sub2ind(size(m), destins, sources);
% rates   = 1;
% m3(inds) = m3(inds) + rates;
% 
% % --- ART initiation
% sources = s.h1;
% destins = s.hart;
% inds    = sub2ind(size(m), destins, sources);
% rates   = r.ART_init;
% m(inds) = m(inds) + rates;

% --- Bring them together
M.lin    = sparse(m - diag(sum(m,1)));
M.Dxlin  = sparse(m2 - diag(sum(m2,1)));
M.linHIV = sparse(m3 - diag(sum(m3,1)));                                   % <--- NEW, now HIV acquisition is separate matrix made in lines 62-67
M.linLTBI = sparse(m4 - diag(sum(m4,1)));

% --- Get the nonlinear rates ---------------------------------------------

% --- Allocating transitions
for istr = 1:length(gps.strains)
    strain = gps.strains{istr};
      
  for ia = 1:length(gps.age)
      age  = gps.age{ia};
    
    m = zeros(i.nstates);
    for ih = 1:length(gps.hiv)
        hiv = gps.hiv{ih};
        
        for iv = 1:length(gps.vacc)
            vacc = gps.vacc{iv};
            
            sources = [i.U.(age).(hiv).(vacc), intersect(intersect(intersect(intersect([s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(age)),s.(hiv)),s.(strain)),s.(vacc))];
            destin  = i.Lf.(age).(hiv).(vacc).(strain);
            m(destin, sources) = m(destin, sources) + (1-(iv==2)*p.VE(1));
            
            %M.nlin.(hiv).(strain) = sparse(m - diag(sum(m,1)));
        end
        
        % Adjust for any immune protection
        cols = intersect([s.Lf, s.Ls, s.Rlo, s.Rhi, s.R], s.(hiv));
        m(:,cols) = m(:,cols)*(1-p.imm(ih));
    end
    M.nlin.(age).(strain) = sparse(m - diag(sum(m,1)));
  end
end


% construct basic building block of repeating contact matrix, as per five-age category
tmp = zeros(length(gps.age),i.nstates);
for ia_sus = 1:length(gps.age)
    for ia_inf = 1:length(gps.age)
        cols = intersect(s.infectious,s.(gps.age{ia_inf}));
        tmp(ia_sus,cols) = prm.cm(ia_sus,ia_inf);
    end
end    

% --- Getting force-of-infection
m = zeros(10,i.nstates);            % DS/mdr and 5 age-groups
m(1:5,intersect(s.infectious, s.ds))   = r.beta*tmp(:,intersect(s.infectious, s.ds));
if length(gps.strains)>1
   m(6:10,intersect(s.infectious, s.mdr))  = r.beta*p.rfbeta_mdr*tmp(:,intersect(s.infectious, s.mdr)); %r.beta*p.rfbeta_mdr;
end
if length(gps.hiv)>1
   m(:,[s.h1, s.hart]) = m(:,[s.h1, s.hart])*p.rfbeta_hiv;
end 
m(:,s.Isc) = m(:,s.Isc)*0.75;
M.lambda = sparse(m);

 
% --- Get the mortality rates
m = zeros(i.nstates,5);
%m(:,1)            = r.mort;
m(s.a1,1)         = r.mort(1);
m(s.a2,1)         = r.mort(2);
m(s.a3,1)         = r.mort(3);
m(s.a4,1)         = r.mort(4);
m(s.a5,1)         = r.mort(5);

%m(s.h1,1)         = r.HIV_mort;                                            % HIV+ve, no TB
% HIV-ve TB
inds = intersect([s.h0], s.infectious_wosc);
m(inds,2) = r.mort_TB(1);
% % HIV/TB coinfection
% inds = intersect([s.h1,s.hart], s.infectious_wosc);
m(inds,3) = r.mort_TB(2);
%keyboard

M.mortvec = m;