function result = test_smoke()
%TEST_SMOKE Self-contained MATLAB GSVD, stationarity, and three-entry check.
here=fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here),'Jpsi_vector','matlab'));
addpath(fullfile(fileparts(here),'eta_c_axial_vector','matlab'));
addpath(fullfile(fileparts(here),'eta_c_pseudoscalar','matlab'));
rng(2718,'twister');
s=[15,15.7,16.9,18]; t=linspace(2,5,17).';
% Two unpenalized poles and four continuum nodes, with a nonzero endpoint.
K=randn(numel(t),6); g=randn(numel(t),1);
H=charmonium.hMassStiffness(s); boundary=[0;0.37];
fac=charmonium.factorsOneWindow(K,g,t,H,2,boundary);
assert(numel(fac.beta_n)==2 && max(abs(fac.mu-1))>0.05);
% Assemble the augmented reference independently of the returned factors.
sw=sqrt(charmonium.trapweights(t));
A=sw.*K(:,[1,2,4,5]); Hii=H(2:3,2:3); Hif=H(2:3,[1,4]);
q=Hii\(Hif*boundary); L=[zeros(2),chol(Hii)];
b=sw.*(g-K(:,[3,6])*boundary)+A(:,3:end)*q;
alpha=[0.01,0.3,2];
result.contract=struct('direct',zeros(size(alpha)),'gradient',zeros(size(alpha)), ...
    'mu_range',[min(fac.mu),max(fac.mu)],'identity', ...
    [fac.identity_relative_error_A,fac.identity_relative_error_B]);
for j=1:numel(alpha)
    coef=(fac.sigma.*fac.beta_c)./(fac.sigma.^2+alpha(j).*fac.mu.^2);
    un=fac.beta_n./fac.sigma_unregularized;
    shifted=fac.X'\[coef;un];
    direct=[A;sqrt(alpha(j))*L]\[b;zeros(size(L,1),1)];
    continuum=[boundary(1);shifted(3:end)-fac.q;boundary(2)];
    weightedResidual=sw.*(K*[shifted(1:2);continuum]-g);
    gradient=A'*weightedResidual+ ...
        alpha(j)*[zeros(2,1);Hii*continuum(2:3)+Hif*boundary];
    scale=max(norm(A'*b)+alpha(j)*norm(L,'fro')^2*norm(shifted),realmin);
    result.contract.direct(j)=norm(shifted-direct)/max(norm(direct),realmin);
    result.contract.gradient(j)=norm(gradient)/scale;
    assert(result.contract.direct(j)<1e-11);
    assert(result.contract.gradient(j)<1e-11);
end
assert(max(fac.identity_relative_error_A,fac.identity_relative_error_B)<1e-11);
solvers={@solve_jpsi_vector,@solve_eta_c_axial_vector,@solve_eta_c_pseudoscalar};
expectedMasses={[3.096900;3.686097],[2.9841;3.51067],[2.9841;3.6378]};
expectedFields={'two_pole';'metadata';'M2';'weights';'OPE'};
result.physical=struct('direct',[],'gradient',[],'factor',[]);
for k=1:numel(solvers)
    out=solvers{k}(5,20,20,'CheckSolution',true);
    assert(isequal(sort(fieldnames(out)),sort(expectedFields)));
    assert(isequal(out.metadata.channels,{'2pole'}) && out.metadata.rg.loops==4);
    ch=out.two_pole; fo=ch.first_order;
    assert(isequal(ch.masses,expectedMasses{k}));
    assert(numel(ch.factors.beta_n)==2);
    assert(all(isnan(ch.f(ch.residues<=0))));
    factorError=max(ch.factors.identity_relative_error_A, ...
        ch.factors.identity_relative_error_B);
    result.physical.direct(end+1)=fo.relative_direct_solution_error_inf;
    result.physical.gradient(end+1)=fo.relative_direct_gradient_inf;
    result.physical.factor(end+1)=factorError;
    assert(fo.relative_direct_solution_error_inf<1e-10);
    assert(fo.relative_direct_gradient_inf<1e-10);
    assert(factorError<1e-10);
    explicit=solvers{k}(5,20,20,'Channel','2pole');
    assert(isequal(sort(fieldnames(explicit)),sort(expectedFields)));
    assert(isequaln(ch.residues,explicit.two_pole.residues));
    assert(ch.alpha_star==explicit.two_pole.alpha_star);
    assert(isnan(explicit.two_pole.first_order.relative_direct_solution_error_inf));
end
options=charmonium.commonOptions();
assert(strcmp(options.Channel,'2pole'));
rejected=false;
try
    charmonium.commonOptions('Channel','invalid');
catch
    rejected=true;
end
assert(rejected,'An unsupported channel must be rejected');
result.max_physical_direct=max(result.physical.direct);
result.max_physical_gradient=max(result.physical.gradient);
result.max_physical_factor=max(result.physical.factor);
end
