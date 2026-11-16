function [kB,diagnostics] = pp_bw_borel_column(m2,Lambda,model,varargin)
%PP_BW_BOREL_COLUMN Adaptive Borel integral, independent of the hat grid.
% Theta flattens the pole; physical thresholds are explicit breakpoints.
% The finite BW area is recorded, never used to renormalize the spectrum.
p=inputParser; p.PartialMatching=false;
validTol=@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
p.addParameter('RelTol',2e-12,validTol); p.addParameter('AbsTol',2e-14,validTol);
p.parse(varargin{:}); o=p.Results;
validateattributes(m2,{'numeric'},{'real','finite','positive','vector'});
validateattributes(Lambda,{'numeric'},{'real','finite','scalar','>',model.S_H});
m2=double(m2(:)); width=model.mass*model.gamma0; center=model.mass^2;
points=unique([model.S_H,model.thresholds_s,center,double(Lambda)]);
points=points(points>=model.S_H & points<=Lambda);
theta=atan((points-center)/width);
value=zeros(numel(m2)+1,1); nIntervals=numel(theta)-1;
for j=1:nIntervals
    value=value+integral(@integrand,theta(j),theta(j+1), ...
        'ArrayValued',true,'RelTol',o.RelTol,'AbsTol',o.AbsTol/nIntervals);
end
kB=value(1:end-1);
diagnostics=struct('method','adaptive theta integral with physical-threshold breakpoints', ...
    'relative_tolerance',o.RelTol,'absolute_tolerance',o.AbsTol, ...
    'breakpoints_s',points,'area',value(end), ...
    'pole_inside',model.S_H<center&&center<Lambda,'finite_area_renormalization',false);
assert(all(isfinite(kB))&&all(kB>0),'Non-finite or nonpositive Borel BW column');

    function valueAtTheta = integrand(th)
        s=center+width*tan(th);
        shape=pp_bw_spectrum(s,model);
        density=shape*width/(cos(th)^2);
        valueAtTheta=[exp(-s./m2)./m2;1]*density;
    end
end
