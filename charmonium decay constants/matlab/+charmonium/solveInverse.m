function out = solveInverse(M2min,M2max,Lambda,phys,varargin)
%SOLVEINVERSE Shared compact-GSVD decay-constant inversion for three channels.
o=charmonium.commonOptions(varargin{:});
lo=double(M2min); hi=double(M2max); lam=double(Lambda);
assert(all(isfinite([lo,hi,lam]))&&lo>0&&hi>lo&&lam>0, ...
    'Require positive M2min<M2max and Lambda');
names={o.Channel};
for name=names.'
    assert(any(strcmp(phys.channel_names,name{1})),'Invalid channel %s',name{1});
end
rg=charmonium.runDec(o.Mu,'Loops',o.RGLoops, ...
    'AlphaMZ',o.RGInputs.alpha_mz,'MZ',o.RGInputs.mz, ...
    'Mb',o.RGInputs.mb,'McBar',o.RGInputs.mbar_c);
mc=rg.mc; alphaS=rg.alpha_s;
sHat=charmonium.spectralGrid(lam,phys.S_H,o.DS); sHat=sHat(:);
t=charmonium.borelGrid(lo,hi); t=t(:); weights=charmonium.trapweights(t);
for name=names.'
    mass=charmonium.massForChannel(phys,name{1});
    assert(lam>phys.S_H&&lam>max(mass.^2),'Lambda must exceed S_H and pole masses squared');
    assert(numel(t)>=numel(sHat)-2+numel(mass),'Borel window has too few free unknowns');
end
if strcmp(phys.kind,'aa')
    [g,pert,d4,d6]=charmonium.aaOpe(lam,t,mc,alphaS,phys.has_nlo,phys.b3_cross);
    ope=struct('total',g,'pert',pert,'D4',d4,'D6',d6);
    boundary=[0;charmonium.aaRhoPert(lam,mc,alphaS,phys.has_nlo)];
else
    [g,loOpe,nloOpe,d4,d6]=charmonium.ppOpe(lam,t,mc,alphaS);
    ope=struct('LO',loOpe,'NLO',nloOpe,'pert',loOpe+nloOpe, ...
        'D4',d4,'D6',d6,'total',g);
    boundary=[0;charmonium.ppRhoPert(lam,mc,alphaS)];
end
Khat=charmonium.hatKernel(t,sHat); H=charmonium.hMassStiffness(sHat);
assert(all(isfinite([g(:);Khat(:);H(:)])), 'Non-finite operator or OPE');
for j=1:numel(names)
    name=names{j}; channelIndex=find(strcmp(phys.channel_names,name),1);
    field=phys.field_names{channelIndex};
    masses=charmonium.massForChannel(phys,name); masses=masses(:); p=numel(masses);
    poles=(exp(-masses.^2./t.')./t.').';
    K=[poles,Khat];
    fac=charmonium.factorsOneWindow(K,g,t,H,p,boundary);
    if strcmp(phys.kind,'aa')
        if isempty(o.Alpha), [alpha,covered]=charmonium.selectionGrid( ...
            fac.sigma_tilde,fac.beta_c,fac.null);
        else, alpha=o.Alpha; covered=false; end
    else
        if isempty(o.Alpha), alpha=10.^linspace(-20,0,2001); else, alpha=o.Alpha; end
        covered=true;
    end
    sel=charmonium.selectRaus(fac.sigma_tilde,fac.beta_c,fac.null,alpha);
    a=sel.alpha; nReg=numel(fac.sigma);
    coef=(fac.sigma.*fac.beta_c)./(fac.sigma.^2+a*fac.mu.^2);
    un=fac.beta_n./fac.sigma_unregularized;
    shifted=fac.X'\[coef;un];
    residues=shifted(1:p);
    cont=[boundary(1);shifted(p+1:end)-fac.q;boundary(2)];
    solution=[residues;cont];
    backsub=K*solution; residual=backsub-g;
    wr=sqrt(sum(weights.*residual.^2));
    penalty=sqrt(max(cont'*H*cont,0));
    positive=residues>0; f=nan(p,1); f(positive)=sqrt(residues(positive));
    if strcmp(phys.kind,'pp')
        f(positive)=2*mc.*sqrt(residues(positive))./(masses(positive).^2);
    end
    diagnostics=sel;
    diagnostics.grid=alpha; diagnostics.on_boundary=any(sel.index==[1,numel(alpha)]);
    diagnostics.small_alpha_peak_covered=covered;
    diagnostics.sigma_min_sq=min(fac.sigma_tilde)^2;
    if strcmp(phys.kind,'pp')
        den=fac.sigma.^2+alpha.*fac.mu.^2;
        coefAll=(fac.sigma.*fac.beta_c)./den;
        derivAll=-(log(10)*alpha.*fac.sigma.*fac.beta_c.*fac.mu.^2)./den.^2;
        unAll=repmat(fac.beta_n./fac.sigma_unregularized,1,numel(alpha));
        shiftedAll=fac.X'\[coefAll;unAll];
        dshifted=fac.X'\[derivAll;zeros(p,numel(alpha))];
        residuesAll=shiftedAll(1:p,:);
        rhoAll=[repmat(boundary(1),1,numel(alpha)); ...
            shiftedAll(p+1:end,:)-fac.q; ...
            repmat(boundary(2),1,numel(alpha))];
        fAll=nan(p,numel(alpha)); positiveAll=residuesAll>0;
        massGrid=masses*ones(1,numel(alpha));
        fAll(positiveAll)=2*mc.*sqrt(residuesAll(positiveAll))./ ...
            massGrid(positiveAll).^2;
        dcont=dshifted(p+1:end,:);
        interiorW=charmonium.trapweights(sHat); interiorW=interiorW(2:end-1);
        qoc=interiorW.'*dcont.^2;
        qmins=charmonium.localMinima(qoc);
        if isempty(qmins), ruleB=-1; else, ruleB=qmins(charmonium.argmin(qoc(qmins))); end
        diagnostics.rule_b_index=ruleB; diagnostics.qoc_argmin_index=charmonium.argmin(qoc);
        diagnostics.qoc=qoc;
    end
    chan=struct('masses',masses,'s',sHat,'s_hat',sHat,'rho',cont, ...
        'hat',cont,'residues',residues,'f',f,'boundary',boundary, ...
        'shifted',shifted, ...
        'alpha_star',a,'alpha_diagnostics',diagnostics,'backsub',backsub, ...
        'backsub_poles',poles*residues,'backsub_continuum',Khat*cont, ...
        'residual',residual,'weighted_residual',wr, ...
        'relative_weighted_residual',wr/max(fac.data_norm,realmin), ...
        'unweighted_residual',norm(residual),'penalty',penalty, ...
        'objective',wr^2+a*penalty^2,'min_rho',min(cont), ...
        'negative_rho_count',sum(cont<0),'factors',fac);
    if strcmp(phys.kind,'pp')
        chan.alpha=alpha; chan.residues_all=residuesAll;
        chan.rho_all=rhoAll; chan.f_all=fAll; chan.qoc=qoc;
    end
    chan.(phys.f_field)=f(1);
    betaAll=[fac.beta_c;fac.beta_n]; y=[coef;un];
    Cy=fac.C*y-betaAll; Sy=fac.S*y;
    stationarity=fac.C'*Cy+a*(fac.S'*Sy);
    transformError=fac.X'*shifted-y;
    directSolutionError=NaN; relativeDirectSolutionError=NaN;
    relativeDirectGradient=NaN;
    if o.CheckSolution
        direct=[fac.A;sqrt(a)*fac.B]\[fac.b0;zeros(size(fac.B,1),1)];
        directGradient=fac.A'*(fac.A*shifted-fac.b0)+ ...
            a*(fac.B'*(fac.B*shifted));
        gradientScale=norm(fac.A'*fac.A,inf)*norm(direct,inf)+ ...
            a*norm(fac.B'*fac.B,inf)*norm(direct,inf);
        directSolutionError=norm(direct-shifted,inf);
        relativeDirectSolutionError=directSolutionError/ ...
            max(norm(direct,inf),realmin);
        relativeDirectGradient=norm(directGradient,inf)/ ...
            max(gradientScale,realmin);
    end
    firstOrder=struct( ...
        'weighted_residual_formula',sqrt(sum(Cy.^2)+fac.null), ...
        'penalty_formula',sqrt(sum(Sy.^2)), ...
        'stationarity_inf',norm(stationarity,inf), ...
        'transform_residual_inf',norm(transformError,inf), ...
        'relative_stationarity_inf',norm(stationarity,inf)/ ...
            max(norm(fac.C'*fac.C*y,inf),realmin), ...
        'relative_transform_residual_inf',norm(transformError,inf)/ ...
            max(norm(y,inf),realmin), ...
        'direct_solution_error_inf',directSolutionError, ...
        'relative_direct_solution_error_inf',relativeDirectSolutionError, ...
        'relative_direct_gradient_inf',relativeDirectGradient);
    chan.first_order=firstOrder;
    out.(field)=chan;
end
out.metadata=struct('schema',phys.schema,'M2_min',lo,'M2_max',hi, ...
    'Lambda',lam,'n_hats',numel(sHat),'ds',o.DS,'dM2',0.01,'mu',o.Mu, ...
    'mc',mc,'alpha_s',alphaS,'G2',0.038,'G3',0.013,'S_H',phys.S_H, ...
    'ope_order',phys.ope_order,'rule',phys.rule,'channels',{names}, ...
    'rg',rg,'gsvd',struct('source','MATLAB built-in gsvd','economy_qr',true, ...
    'convention',"A=U*C*X', B=V*S*X'"));
out.M2=t; out.weights=weights; out.OPE=ope;
end
