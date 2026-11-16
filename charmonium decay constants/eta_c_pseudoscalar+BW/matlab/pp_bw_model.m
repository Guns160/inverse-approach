function model = pp_bw_model(channel,varargin)
%PP_BW_MODEL One running-width eta_c BW state, with inclusive partial widths.
% All masses/widths are in GeV; Radius is in GeV^-1. 6S is reserved only.
if nargin<1, channel='5S'; end
registration=pp_bw_state_info(channel); channel=registration.channel;
if strcmp(registration.status,'reserved_unavailable')
    error('charmonium:BW6SUnavailable','%s',registration.blocker);
end
p=inputParser; p.PartialMatching=false;
p.addParameter('Radius',4,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0);
p.addParameter('DstarD1L',0,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&any(x==[0,2]));
p.addParameter('QuadratureNodes',1024,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&& ...
    isfinite(x)&&x==fix(x)&&x>=16);
p.parse(varargin{:}); o=p.Results;
model=struct('channel',channel,'status','available','usage',registration.usage,'S_H',14.99, ...
    'mass',registration.mass,'reference_mass',registration.reference_mass, ...
    'registration',registration,'truncation_note',registration.truncation_note, ...
    'radius',double(o.Radius),'DstarD1L',double(o.DstarD1L), ...
    'quadrature_nodes',double(o.QuadratureNodes), ...
    'D',[1.86484,1.86961],'Dstar',[2.00696,2.01026], ...
    'Ds',1.9683,'Dsstar',2.1121,'D0star',[2.318,2.403], ...
    'D1_2420',[2.4214,2.4232],'D2star',[2.4626,2.4643], ...
    'Ds0star',2.3177,'D0_width',0.2292725269223561, ...
    'D1_width',0.03133891996221452,'pion',[0.1349768,0.13957039]);
switch channel
    case '3S'
        model.open_width=0.06458; model.rest_width=0.00916374;
        model.rest_components_GeV=struct('gg',0.00916,'gamma_gamma',0.00000374);
        model.raw_partial_widths=0.06458;
        model.channel_names={'DDstar'}; model.orbital_L=1;
        model.channel_thresholds_s={(model.D+model.Dstar).^2};
    case '4S'
        model.open_width=0.02727; model.rest_width=0.00629257;
        model.rest_components_GeV=struct('gg',0.00629,'gamma_gamma',0.00000257);
        model.raw_partial_widths=[0.00012,0.00077,0.02610,0.00028];
        model.channel_names={'DDstar','DsDsstar','DstarDstar','DsstarDsstar'};
        model.orbital_L=[1,1,1,1];
        model.channel_thresholds_s={(model.D+model.Dstar).^2, ...
            (model.Ds+model.Dsstar)^2,(2*model.Dstar).^2,(2*model.Dsstar)^2};
    case '5S'
        model.open_width=0.07515; model.rest_width=0.00437178;
        model.rest_components_GeV=struct('gg',0.00437,'gamma_gamma',0.00000178);
        model.raw_partial_widths=[0.00292,0.00021,0.03846,0.00107, ...
            0.00092,0.00366,0.00013,0.02778];
        model.channel_names={'DDstar','DsDsstar','DD0star','DsDs0star', ...
            'DD2star','DstarDstar','DsstarDsstar','DstarD1'};
        model.orbital_L=[1,1,0,0,2,1,1,model.DstarD1L];
        broadD=physicalThresholds(model.D,model.D,model.pion);
        broadV=physicalThresholds(model.Dstar,model.Dstar,model.pion);
        model.channel_thresholds_s={(model.D+model.Dstar).^2, ...
            (model.Ds+model.Dsstar)^2,broadD,(model.Ds+model.Ds0star)^2, ...
            (model.D+model.D2star).^2,(2*model.Dstar).^2, ...
            (2*model.Dsstar)^2,broadV};
end
model.gamma0=model.open_width+model.rest_width;
model.thresholds_s=unique([model.channel_thresholds_s{:}]);
model.phi_pole=pp_bw_phase_space(model.mass^2,model);
model.phi_reference=pp_bw_phase_space(model.reference_mass^2,model);
assert(all(model.phi_pole>0)&&all(model.phi_reference>0),'Nonpositive pole/reference phase space');
shifted=model.raw_partial_widths.*model.phi_pole./model.phi_reference;
model.pole_partial_widths=model.open_width*shifted/sum(shifted);
model.form_factor='B0^2=1; B1^2=1/(1+z^2); B2^2=1/(9+3*z^2+z^4); z=q*R';
model.charge_policy='equal neutral/charged groups; conjugates not counted twice';
model.spectral_policy='physical D(*)pi lower limits; constant daughter widths; no finite-area renormalization';
model.width_policy='REF-to-pole phase-space shift, then normalize inclusive partial widths to OPEN';
model.broad_width_policy='charged D0star/D1 widths equal neutral widths (isospin approximation)';
model.barrier_location='inside daughter spectral convolution';
model.rest_width_definition='constant gg + gamma_gamma';
model.partial_width_source='EPJC 75 (2015) 196, Table 3, model b';
model.sources={'https://doi.org/10.1140/epjc/s10052-015-3418-5', ...
    'https://pdg.lbl.gov/2026/', ...
    'https://pdg.lbl.gov/2023/reviews/rpp2022-rev-resonances.pdf', ...
    'https://arxiv.org/abs/1608.01289'};
model.units=struct('mass','GeV','width','GeV','s','GeV^2', ...
    'radius','GeV^-1','bw_shape','GeV^-2');
end

function thresholds = physicalThresholds(spectators,daughters,pion)
thresholds=zeros(1,4);
for k=1:2
    thresholds(2*k-1)=(spectators(k)+daughters(3-k)+pion(2))^2;
    thresholds(2*k)=(spectators(k)+daughters(k)+pion(1))^2;
end
end
