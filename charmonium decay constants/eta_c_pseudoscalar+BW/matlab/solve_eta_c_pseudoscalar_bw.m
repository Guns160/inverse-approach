function out = solve_eta_c_pseudoscalar_bw(M2min,M2max,Lambda,varargin)
%SOLVE_ETA_C_PSEUDOSCALAR_BW Two delta poles plus ONE running-width BW state.
% Channel: '5S' (default), '3S', '4S', or reserved/unavailable '6S'.
% Alpha: positive scalar or increasing vector; default 10.^linspace(-20,0,2001).
% No alpha selection, positivity constraints, spectral clipping, or plotting.
% DS=0.1 GeV^2; Mu=2 GeV; RGLoops=4; Radius=4 GeV^-1;
% DstarD1L=0 (or 2); QuadratureNodes=1024 for broad-daughter convolutions.
% Anomaly1S=-0.079 GeV^3, Anomaly2S=[]; the BW state is density-only.
% c_all/hat_all is the remaining continuum (GeV^2); rho_all=c_all./s.
here=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(here)),'matlab'));
o=options(varargin{:});
% Fail on 6S before RG, grids, quadratures, or any inverse-problem objects.
model=pp_bw_model(o.Channel,'Radius',o.Radius,'DstarD1L',o.DstarD1L, ...
    'QuadratureNodes',o.QuadratureNodes);
validateattributes(M2min,{'numeric'},{'real','finite','positive','scalar'});
validateattributes(M2max,{'numeric'},{'real','finite','scalar','>',M2min});
validateattributes(Lambda,{'numeric'},{'real','finite','scalar','>',model.S_H});
lo=double(M2min); hi=double(M2max); lam=double(Lambda);
deltaMass=[2.9841;3.6378];
assert(lam>max(deltaMass.^2),'Lambda must exceed the two delta-pole masses squared');
% Lambda may truncate the BW peak; it need not exceed model.mass^2.
s=charmonium.spectralGrid(lam,model.S_H,o.DS); s=s(:);
t=charmonium.borelGrid(lo,hi); t=t(:); weights=charmonium.trapweights(t);
assert(numel(t)>=numel(s)+1,'BW model needs at least N+1 Borel data for N+1 free unknowns');
rg=charmonium.runDec(o.Mu,'Loops',o.RGLoops, ...
    'AlphaMZ',o.RGInputs.alpha_mz,'MZ',o.RGInputs.mz, ...
    'Mb',o.RGInputs.mb,'McBar',o.RGInputs.mbar_c);
mc=rg.mc; alphaS=rg.alpha_s;
[g,loOpe,nloOpe,d4,d6]=charmonium.ppOpe(lam,t,mc,alphaS);
G=charmonium.hatKernel(t,s); H=charmonium.hMassStiffness(s);
boundary=[0;charmonium.ppRhoPert(lam,mc,alphaS)];
P=(exp(-deltaMass.^2./t.')./t.').';
[kB,integration]=pp_bw_borel_column(t,lam,model);
[shape,gamma]=pp_bw_spectrum(s,model);
system=pp_bw_assemble(P,G,kB,g,weights,H,shape([1,end]),boundary(2));
if isempty(o.Alpha), alpha=10.^linspace(-20,0,2001); else, alpha=o.Alpha; end
[channel,factors]=pp_bw_solve_family(system,alpha);
masses=[deltaMass;model.mass]; residues=channel.residues_all;
[fDelta,rawDelta,shift,anomaly]=charmonium.ppDecayConstants( ...
    residues(1:2,:),deltaMass,mc,o.Anomaly1S,o.Anomaly2S);
rawBW=nan(1,numel(alpha)); positive=residues(3,:)>0;
rawBW(positive)=2*mc*sqrt(residues(3,positive))/model.mass^2;
channel.s=s; channel.s_hat=s; channel.masses=masses;
channel.f_all=[fDelta;rawBW]; channel.f_raw_all=[rawDelta;rawBW];
channel.anomaly_shift_GeV=[shift;0];
channel.hat_all=channel.c_all; channel.rho_all=channel.c_all./s;
channel.bw_shape=shape; channel.gamma=gamma;
channel.bw_hat_all=shape*residues(3,:);
channel.total_hat_all=channel.c_all+channel.bw_hat_all;
channel.rho_total_all=channel.total_hat_all./s;
channel.k_B=kB; channel.boundary_total=boundary; channel.factors=factors;
channel.validation_pass_all=max(abs(channel.endpoint_error_all),[],1)<= ...
    1e-10*(1+abs(boundary(2))) & max(channel.stationarity_relative_all,[],1)<1e-8;
channel.negative_c_count_all=sum(channel.c_all<0,1);
anomaly.bw_GeV3=[]; anomaly.applied=[anomaly.applied;false];
anomaly.f_definition=[anomaly.f_definition;{'pseudoscalar_density_only'}];
field=['bw',lower(model.channel)]; out.(field)=channel;
out.metadata=struct('schema','eta_c-pp-bw-matlab-v1','M2_min',lo,'M2_max',hi, ...
    'Lambda',lam,'n_hats',numel(s),'ds',o.DS,'dM2',0.01,'mu',o.Mu, ...
    'mc',mc,'alpha_s',alphaS,'G2',0.038,'G3',0.013,'S_H',model.S_H, ...
    'ope_order','PP n=0 LO+NLO(alpha_s)+G2+G3','channel',model.channel, ...
    'alpha_selection','none','rg',rg,'anomaly',anomaly,'bw',model, ...
    'integration',integration,'residue_unit','GeV^4','f_unit','GeV', ...
    'continuum_hat_unit','GeV^2','rho_unit','dimensionless', ...
    'penalty','c^T H c; remaining continuum only', ...
    'endpoints','c(S_H)+R_BW*BW(S_H)=0; c(Lambda)+R_BW*BW(Lambda)=rho_pert(Lambda)', ...
    'gsvd',struct('source','MATLAB built-in gsvd','economy_qr',true, ...
    'convention',"A=Q*U*C*X', B=V*S*X'; x=X'\coeff"));
out.M2=t; out.weights=weights;
out.OPE=struct('LO',loOpe,'NLO',nloOpe,'pert',loOpe+nloOpe,'D4',d4,'D6',d6,'total',g);
end

function o = options(varargin)
p=inputParser; p.PartialMatching=false;
p.addParameter('Channel','5S',@(x)((ischar(x)&&isrow(x))|| ...
    (isstring(x)&&isscalar(x)))&&any(strcmpi(x,{'3S','4S','5S','6S'})));
positive=@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
p.addParameter('DS',0.1,positive); p.addParameter('Mu',2,positive);
p.addParameter('RGLoops',4,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&& ...
    isfinite(x)&&x==fix(x)&&x>=1&&x<=4);
p.addParameter('RGInputs',struct('mbar_c',1.275,'alpha_mz',0.1184, ...
    'mz',91.1876,'mb',4.18),@charmonium.validRGInputs);
validAnomaly=@(x)isnumeric(x)&&isreal(x)&& ...
    (isempty(x)||(isscalar(x)&&isfinite(x)));
p.addParameter('Anomaly1S',-0.079,validAnomaly);
p.addParameter('Anomaly2S',[],validAnomaly);
p.addParameter('Alpha',[],@(x)isnumeric(x)&&isreal(x)&&isvector(x)&&~isempty(x)&& ...
    all(isfinite(x(:)))&&all(x(:)>0)&&all(diff(x(:))>0));
p.addParameter('Radius',4,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0);
p.addParameter('DstarD1L',0,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&any(x==[0,2]));
p.addParameter('QuadratureNodes',1024,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&& ...
    isfinite(x)&&x==fix(x)&&x>=16);
p.parse(varargin{:}); o=p.Results;
o.Channel=upper(char(o.Channel));
if ~isempty(o.Alpha), o.Alpha=double(o.Alpha(:).'); end
end
