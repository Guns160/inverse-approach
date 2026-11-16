function test_pp_bw(scope)
%TEST_PP_BW Affine, physical-width, integration, and public one-window tests.
% test_pp_bw('metadata') checks only the physical input registry/metadata.
here=fileparts(mfilename('fullpath')); addpath(here);
addpath(fullfile(fileparts(here),'eta_c_pseudoscalar+BW','matlab'));
if nargin>0
    assert(strcmp(scope,'metadata'),'Only the metadata subset is supported');
    testMetadata(); return;
end
testMetadata(); testAffine(); testWidths(); testBorel(); testOptions(); testWindows();
disp('PP BW: affine GSVD, widths, adaptive Borel, options and 3S/4S/5S windows PASS');
end

function testMetadata()
channels={'3S','4S','5S'};
usage={'diagnostic: too close to the two delta poles; poor extraction', ...
    'alternative','primary'};
gg=[0.00916,0.00629,0.00437]; photons=[0.00000374,0.00000257,0.00000178];
for j=1:3
    info=pp_bw_state_info(channels{j}); model=pp_bw_model(channels{j});
    assert(strcmp(info.status,'available')&&strcmp(model.status,'available'));
    assert(strcmp(info.usage,usage{j})&&strcmp(model.usage,usage{j}));
    assert(model.rest_components_GeV.gg==gg(j));
    assert(model.rest_components_GeV.gamma_gamma==photons(j));
    assertClose(model.rest_width,gg(j)+photons(j),2e-15,2e-17);
    assert(strcmp(model.rest_width_definition,'constant gg + gamma_gamma'));
    assert(strcmp(model.partial_width_source,'EPJC 75 (2015) 313, Table 3, model b'));
    assert(ismember('https://doi.org/10.1140/epjc/s10052-015-3418-5',model.sources));
    assert(ismember('https://pdg.lbl.gov/2026/',model.sources));
end
reserved=pp_bw_state_info('6S');
assert(strcmp(reserved.status,'reserved_unavailable')&&strcmp(reserved.usage,'reserved'));
assertReject(@()pp_bw_model('6S'),'charmonium:BW6SUnavailable');
disp('  metadata: state usage, public sources, gg+gamma_gamma residual width and 6S reservation PASS');
end

function testAffine()
saved=rng; cleanup=onCleanup(@()rng(saved));
rng(741); n=9; m=35;
P=randn(m,2); G=randn(m,n); k=randn(m,1); w=0.2+rand(m,1);
H=charmonium.hMassStiffness(linspace(14.99,20,n));
g=-3*P(:,1)+0.1*randn(m,1);
sys=pp_bw_assemble(P,G,k,g,w,H,[0.03;0.07],0.52);
assert(isequal(size(sys.E),[n,n-1])&&isequal(size(sys.A),[m,n+1]));
assert(isequal(sys.E(2:end-1,2:end),eye(n-2)));
assert(isequal(sys.E([1,end],1),-[0.03;0.07]));
assert(nnz(all(sys.B==0,1))==2&&norm(sys.B(:,3))>0);
assert(sys.kappa>0&&sys.qr_relative_error<1e-14);
alpha=[1e-8,1e-4,0.1,1]; [family,fac]=pp_bw_solve_family(sys,alpha);
assert(fac.n_regularized==n-1&&fac.n_unregularized==2);
for j=1:numel(alpha)
    % Independent unshifted augmented LS; penalize the FULL affine c=a+E*v.
    L=[zeros(n,2),sys.L_H*sys.E];
    direct=[sqrt(w).*[P,sys.F];sqrt(alpha(j))*L]\ ...
        [sqrt(w).*(g-G*sys.a);-sqrt(alpha(j))*sys.L_H*sys.a];
    directC=sys.a+sys.E*direct(3:end);
    directR=[direct(1:2);direct(3)];
    assertClose(family.residues_all(:,j),directR,2e-11,2e-12);
    assertClose(family.c_all(:,j),directC,2e-11,2e-12);
    x=[family.residues_all(1:2,j); ...
        [family.residues_all(3,j);family.c_all(2:end-1,j)]+sys.q];
    originalPenalty=family.c_all(:,j)'*H*family.c_all(:,j);
    assertClose(originalPenalty,norm(sys.B*x)^2+sys.kappa,2e-12,2e-12);
    assertClose(family.objective_all(j),norm(sqrt(w).*(P*directR(1:2)+ ...
        k*directR(3)+G*directC-g))^2+alpha(j)*originalPenalty,2e-12,2e-12);
end
assert(max(abs(family.endpoint_error_all),[],'all')<1e-14);
assert(max(family.stationarity_relative_all,[],'all')<2e-13);
assert(all(family.residues_all(1,:)<0),'Signed residues must not be clipped');
assert(all(~family.positive_all));
end

function testWidths()
channels={'3S','4S','5S'}; masses=[4.0119,4.2693,4.4508];
rest=[0.00916374,0.00629257,0.00437178]; open=[0.06458,0.02727,0.07515];
counts=[1,4,8];
for j=1:3
    model=pp_bw_model(channels{j});
    assert(model.mass==masses(j)&&model.rest_width==rest(j)&&model.open_width==open(j));
    assert(numel(model.channel_names)==counts(j));
    [shape,gamma,partial]=pp_bw_spectrum(model.mass^2,model);
    assertClose(sum(partial),open(j),3e-15,3e-17);
    assertClose(gamma,open(j)+rest(j),3e-15,3e-17);
    assertClose(shape,1/(pi*model.mass*gamma),3e-15,3e-17);
    shifted=model.raw_partial_widths.*pp_bw_phase_space(model.mass^2,model)./ ...
        pp_bw_phase_space(model.reference_mass^2,model);
    assertClose(partial,open(j)*shifted/sum(shifted),3e-15,3e-17);
    [~,below]=pp_bw_spectrum(1,model); assert(below==rest(j));
    for c=1:counts(j)
        threshold=min(model.channel_thresholds_s{c});
        [~,~,closed]=pp_bw_spectrum(threshold*(1-1e-10),model);
        [~,~,opened]=pp_bw_spectrum(threshold*(1+1e-3),model);
        assert(closed(c)==0&&opened(c)>0,'Wrong threshold in %s/%s',channels{j},model.channel_names{c});
    end
    s=[15;model.mass^2;25]; actual=pp_bw_phase_space(s,model);
    assertClose(actual(:,1),stableAverage(s,model.D,model.Dstar,1,4),2e-14,2e-16);
    if j==2
        expected=[stableAverage(s,model.D,model.Dstar,1,4), ...
            stableReference(s,model.Ds,model.Dsstar,1,4), ...
            stableAverage(s,model.Dstar,model.Dstar,1,4), ...
            stableReference(s,model.Dsstar,model.Dsstar,1,4)];
        assertClose(actual,expected,2e-14,2e-16);
    end
end
model=pp_bw_model('5S'); s=[16.5;model.mass^2;25];
phi=pp_bw_phase_space(s,model);
assertClose(phi(:,3),broadReference(s,model.D,model.D0star,model.D0_width, ...
    model.D,0,model),3e-10,3e-13);
assertClose(phi(:,8),broadReference(s,model.Dstar,model.D1_2420,model.D1_width, ...
    model.Dstar,0,model),3e-10,3e-13);
assertClose(phi(:,4),stableReference(s,model.Ds,model.Ds0star,0,4),2e-14,2e-16);
assertClose(phi(:,5),stableAverage(s,model.D,model.D2star,2,4),2e-14,2e-16);
assertClose(phi(:,6),stableAverage(s,model.Dstar,model.Dstar,1,4),2e-14,2e-16);
assertClose(phi(:,7),stableReference(s,model.Dsstar,model.Dsstar,1,4),2e-14,2e-16);
assert(min(model.channel_thresholds_s{3})<min((model.D+model.D0star).^2));
assert(min(model.channel_thresholds_s{8})<min((model.Dstar+model.D1_2420).^2));
for radius=[0,4]
    dWave=pp_bw_model('5S','DstarD1L',2,'Radius',radius);
    dphi=pp_bw_phase_space(s,dWave);
    assertClose(dphi(:,8),broadReference(s,dWave.Dstar,dWave.D1_2420, ...
        dWave.D1_width,dWave.Dstar,2,dWave),3e-10,3e-13);
end
refined=pp_bw_model('5S','QuadratureNodes',2048);
assertClose(pp_bw_phase_space(s,refined),phi,3e-10,3e-13);
disp('  widths: charge averages, all thresholds, pole sums, S/P/D barriers and GL1024 convergence PASS');
end

function testBorel()
channels={'3S','4S','5S'}; t=[5;10;30];
for j=1:3
    model=pp_bw_model(channels{j}); [k,info]=pp_bw_borel_column(t,20,model);
    reference=zeros(size(t)); points=model.thresholds_s;
    points=unique([points(points>model.S_H&points<20),model.mass^2]);
    for i=1:numel(t)
        reference(i)=quadgk(@(s)exp(-s/t(i))/t(i).*reshape(pp_bw_spectrum(s,model),size(s)), ...
            model.S_H,20,'Waypoints',points,'AbsTol',2e-14,'RelTol',2e-12);
    end
    assertClose(k,reference,3e-10,3e-13);
    assert(~info.finite_area_renormalization&&abs(info.area-1)>1e-4);
end
model=pp_bw_model('5S');
[k,info]=pp_bw_borel_column(t,18,model);
refined=pp_bw_borel_column(t,18,model,'RelTol',2e-13,'AbsTol',2e-15);
assertClose(k,refined,2e-10,2e-13);
assert(~info.pole_inside&&info.area>0&&info.area<1);
disp('  Borel: independent s/angle adaptive integrals, convergence and truncated peak PASS');
end

function testOptions()
reserved=pp_bw_state_info('6S');
assert(reserved.mass==4.5777&&reserved.reference_mass==4.612);
assert(strcmp(reserved.status,'reserved_unavailable'));
assert(contains(reserved.blocker,'4.6472')&&contains(reserved.blocker,'0.35 MeV')&& ...
    contains(reserved.blocker,'DDstar(2S)')&&contains(reserved.blocker,'DstarD1prime'));
assert(~isfield(reserved,'gamma0')&&~isfield(reserved,'pole_partial_widths'));
assertReject(@()solve_eta_c_pseudoscalar_bw(NaN,NaN,NaN,'Channel','6S'), ...
    'charmonium:BW6SUnavailable');
assertReject(@()pp_bw_model('6S'),'charmonium:BW6SUnavailable');
invalid={{'Radius',-1},{'Radius',1i},{'DstarD1L',1},{'QuadratureNodes',15}, ...
    {'QuadratureNodes',16.5},{'Channel','all'},{'Alpha',0},{'Alpha',[]}, ...
    {'Alpha',[1e-3,1e-3]},{'Alpha',[1,0.1]},{'Alpha',1i},{'Alpha',[1,NaN]}};
for j=1:numel(invalid)
    args=invalid{j}; assertReject(@()solve_eta_c_pseudoscalar_bw(5,8,20,args{:}),'');
end
assertReject(@()solve_eta_c_pseudoscalar_bw(5,5.01,20,'DS',1,'Alpha',1e-6),'');
assertReject(@()solve_eta_c_pseudoscalar_bw(5,8,14.99,'Alpha',1e-6),'');
end

function testWindows()
alpha=[1e-10,1e-6,1e-2]; channels={'3S','4S','5S'};
for j=1:3
    out=solve_eta_c_pseudoscalar_bw(5,8,20,'Channel',channels{j},'DS',1,'Alpha',alpha);
    name=['bw',lower(channels{j})]; c=out.(name); s=c.s; n=numel(s);
    assert(numel(fieldnames(out))==5&&strcmp(out.metadata.channel,channels{j}));
    assert(strcmp(out.metadata.schema,'eta_c-pp-bw-matlab-v1'));
    assert(strcmp(out.metadata.alpha_selection,'none')&&~isfield(c,'alpha_star'));
    assert(isequal(size(c.residues_all),[3,3])&&isequal(size(c.c_all),[n,3]));
    assert(isequal(c.hat_all,c.c_all)&&isequal(c.rho_all,c.c_all./s));
    assert(isequal(c.bw_hat_all,c.bw_shape*c.residues_all(3,:)));
    assert(isequal(c.total_hat_all,c.c_all+c.bw_hat_all));
    assert(isequal(c.rho_total_all,c.total_hat_all./s));
    assert(isequal(c.positive_all,all(c.residues_all>0,1)));
    assert(max(abs(c.endpoint_error_all),[],'all')<1e-12);
    assert(max(c.stationarity_relative_all,[],'all')<1e-8);
    assert(all(isnan(c.f_all(c.residues_all<=0))));
    assertClose(c.f_all(3,:),c.f_raw_all(3,:),0,0);
    assertClose(c.anomaly_shift_GeV,[-0.079/c.masses(1)^2;0;0],2e-15,2e-17);
    sys=c.factors.affine;
    assert(nnz(all(sys.B==0,1))==2&&size(sys.B,1)==n-1);
    assertClose(c.backsub_all,sys.P*c.residues_all(1:2,:)+ ...
        c.k_B*c.residues_all(3,:)+sys.G*c.c_all,0,0);
    for a=1:numel(alpha)
        direct=[sys.A;sqrt(alpha(a))*[zeros(n,2),sys.L_H*sys.E]]\ ...
            [sqrt(out.weights).*(out.OPE.total-sys.G*sys.a); ...
            -sqrt(alpha(a))*sys.L_H*sys.a];
        assertClose(c.residues_all(:,a),[direct(1:2);direct(3)],2e-7,2e-8);
        assertClose(c.c_all(:,a),sys.a+sys.E*direct(3:end),2e-7,2e-8);
    end
end
% Scalar alpha, a BW peak outside the integration interval, and Ward options.
base=solve_eta_c_pseudoscalar_bw(5,8,16,'DS',0.5,'Alpha',1e-6);
other=solve_eta_c_pseudoscalar_bw(5,8,16,'DS',0.5,'Alpha',1e-6, ...
    'Anomaly1S',[],'Anomaly2S',0.02);
assert(strcmp(base.metadata.channel,'5S')&&~base.metadata.integration.pole_inside);
assert(isequal(base.bw5s.residues_all,other.bw5s.residues_all));
assert(isequal(base.bw5s.c_all,other.bw5s.c_all));
assert(isequal(other.metadata.anomaly.applied,[false;true;false]));
assertClose(other.bw5s.f_all,other.bw5s.f_raw_all+[0;0.02/3.6378^2;0],0,0);
% The Borel kernel must be independent of the continuum hat step.
fine=solve_eta_c_pseudoscalar_bw(5,8,16,'DS',0.25,'Alpha',1e-6);
assert(isequal(fine.bw5s.k_B,base.bw5s.k_B));
default=solve_eta_c_pseudoscalar_bw(5,6,16,'DS',0.5,'Channel','3S');
assert(isequal(default.bw3s.alpha,10.^linspace(-20,0,2001)));
disp('  windows: 3S/4S/5S backsubstitution, augmented LS, signed f and anomaly options PASS');
end

function value = stableAverage(s,a,b,L,R)
value=(stableReference(s,a(1),b(1),L,R)+stableReference(s,a(2),b(2),L,R))/2;
end

function value = stableReference(s,a,b,L,R)
p=sqrt(max((s-(a+b).^2).*(s-(a-b).^2),0))./(2*sqrt(s));
p=p.*(sqrt(s)>a+b); z2=(p*R).^2;
if L==0, barrier=ones(size(p)); elseif L==1, barrier=1./(1+z2);
else, barrier=1./(9+3*z2+z2.^2); end
value=p.^(2*L+1)./sqrt(s).*barrier;
end

function value = broadReference(s,spectators,poles,width,daughters,L,model)
value=zeros(size(s));
for i=1:numel(s)
    for k=1:2
        lower=[(daughters(3-k)+model.pion(2))^2,(daughters(k)+model.pion(1))^2];
        weights=[2/3,1/3]; upper=(sqrt(s(i))-spectators(k))^2;
        for mode=1:2
            if sqrt(s(i))<=spectators(k)+sqrt(lower(mode)), continue; end
            point=poles(k)^2; point=point(point>lower(mode)&point<upper);
            integralValue=quadgk(@(u)poles(k)*width./ ...
                (pi*((u-poles(k)^2).^2+(poles(k)*width)^2)).* ...
                stableReference(s(i),spectators(k),sqrt(u),L,model.radius), ...
                lower(mode),upper,'Waypoints',point,'RelTol',2e-12,'AbsTol',2e-14);
            value(i)=value(i)+0.5*weights(mode)*integralValue;
        end
    end
end
end

function assertClose(actual,expected,rtol,atol)
assert(isequal(size(actual),size(expected)),'Comparison size mismatch');
assert(isequal(isnan(actual),isnan(expected)),'NaN masks differ');
valid=~isnan(expected);
assert(all(abs(actual(valid)-expected(valid))<=atol+rtol*abs(expected(valid))), ...
    'Numerical comparison failed (max absolute difference %.4g)', ...
    max(abs(actual(valid)-expected(valid)),[],'all'));
end

function assertReject(thunk,identifier)
rejected=false;
try
    thunk();
catch err
    rejected=true;
    if ~isempty(identifier)
        assert(strcmp(err.identifier,identifier),'Unexpected error identifier: %s',err.identifier);
        assert(contains(err.message,'6S')&&contains(err.message,'partial widths')&& ...
            contains(err.message,'unavailable'),'6S error must identify missing physical inputs');
    end
end
assert(rejected,'Expected invalid input to be rejected');
end
