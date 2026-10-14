function out = runDec(mu, varargin)
%RUNDEC Shared charmonium MSbar RG: N-loop running, (N-1)-loop matching.
% mu in [McBar,Mb]; charm anchor is ALWAYS nf=4 mc(McBar)=McBar.
% Bottom is its nf=5 MSbar self mass. No independent MassScale is allowed.
if nargin < 1, mu = 2.0; end
positive = @(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
loopValid = @(x)positive(x)&&x==fix(x)&&x>=1&&x<=4;
p = inputParser; p.FunctionName = 'charmonium.runDec';
p.PartialMatching = false; p.KeepUnmatched = false;
p.addParameter('Loops',4,loopValid);
p.addParameter('AlphaMZ',0.1184,@(x)positive(x)&&x<pi);
p.addParameter('MZ',91.1876,positive);
p.addParameter('Mb',4.18,positive);
p.addParameter('McBar',1.275,positive);
p.parse(varargin{:}); o = p.Results;
validateattributes(mu, {'numeric'}, {'real','scalar','finite','positive'});
mu = double(mu); loops = double(o.Loops);
inputs = struct('mbar_c',double(o.McBar),'alpha_mz',double(o.AlphaMZ), ...
    'mz',double(o.MZ),'mb',double(o.Mb));
if ~(inputs.mbar_c < inputs.mb && inputs.mb < inputs.mz)
    error('charmonium:RGInputs','Require McBar < Mb < MZ');
end
if mu < inputs.mbar_c || mu > inputs.mb
    error('charmonium:RGDomain', ...
        'mu must lie in [McBar,Mb]; mass matching outside this nf=4 domain is not implemented');
end
[a5Mb,~] = charmonium.integrateRG(inputs.alpha_mz,1,inputs.mz,inputs.mb,5,loops);
a4Mb = charmonium.matchAlphaDownMS(a5Mb,5,loops-1);
[aAnchor,~] = charmonium.integrateRG(a4Mb,1,inputs.mb,inputs.mbar_c,4,loops);
[alphaS,mc,integration] = charmonium.integrateRG(aAnchor,inputs.mbar_c,inputs.mbar_c,mu,4,loops);
inputSource = 'Caller-supplied RG inputs';
if inputs.mbar_c == 1.275 && inputs.alpha_mz == .1184 && inputs.mz == 91.1876 && inputs.mb == 4.18
    inputSource = ['Project inputs: mcbar(mcbar)=1.275 GeV; alpha_s^(5)(MZ)=0.1184; ' ...
        'MZ=91.1876 GeV; mbbar(mbbar)=4.18 GeV'];
end
out = struct('mu',mu,'mc',mc,'alpha_s',alphaS,'nf',4, ...
    'loops',loops,'matching_loops',loops-1,'mass_method','ode','inputs',inputs, ...
    'input_source',inputSource, ...
    'coefficient_source','Chetyrkin, Kuhn, Steinhauser, hep-ph/0004189v1, Eqs. (2), (7), (20)-(21)', ...
    'scheme','MSbar','a_convention','alpha_s/pi','t_convention','ln(mu^2)', ...
    'mass_anchor',struct('mu',inputs.mbar_c,'mc',inputs.mbar_c,'alpha_s',aAnchor,'nf',4), ...
    'threshold',struct('mu',inputs.mb,'heavy_mass_scheme','MSbar at own scale', ...
        'nf_high',5,'nf_low',4,'alpha_s_high',a5Mb,'alpha_s_low',a4Mb, ...
        'matching_loops',loops-1),'integration',integration);
end
