function [f,fRaw,shift,anomaly] = ppDecayConstants(residues,masses,mc,a1,a2)
%PPDECAYCONSTANTS f=(2*mc*sqrt(R)+a)/m_H^2; [] omits a state's anomaly.
% The omitted state remains density-only. Nonpositive residues remain NaN.
if nargin<4, a1=-0.079; end
if nargin<5, a2=[]; end
valid=@(x)isnumeric(x)&&isreal(x)&& ...
    (isempty(x)||(isscalar(x)&&isfinite(x)));
assert(valid(a1)&&valid(a2),'Anomaly inputs must be finite real scalars or []');
masses=double(masses(:));
assert(numel(masses)==2&&size(residues,1)==2,'Expected two pole rows');
assert(all(isfinite(masses))&&all(masses>0)&&isscalar(mc)&&isfinite(mc)&&mc>0, ...
    'Require positive finite masses');
values=zeros(2,1); supplied=[~isempty(a1);~isempty(a2)];
if supplied(1), values(1)=double(a1); end
if supplied(2), values(2)=double(a2); end
shift=values./masses.^2;
massGrid=repmat(masses,1,size(residues,2));
positive=residues>0; fRaw=nan(size(residues));
fRaw(positive)=2*mc.*sqrt(residues(positive))./massGrid(positive).^2;
f=fRaw+shift;
definitions=repmat({'pseudoscalar_density_only'},2,1);
definitions(supplied)={'axial_with_supplied_anomaly'};
anomaly=struct('eta_c_1s_GeV3',a1,'eta_c_2s_GeV3',a2, ...
    'applied',supplied,'f_definition',{definitions}, ...
    'default_1s_source','Qin et al., Phys. Rev. D 97, 096002 (2018), Eq. (22)', ...
    'convention','h_i=+sqrt(R_i); supplied a_i uses the same state phase');
end
