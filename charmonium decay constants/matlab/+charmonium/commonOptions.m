function o = commonOptions(varargin)
%COMMONOPTIONS Shared Name=Value options for the three public MATLAB solvers.
p=inputParser;
p.PartialMatching=false;
p.addParameter('Channel','2pole',@(x)((ischar(x)&&isrow(x))|| ...
    (isstring(x)&&isscalar(x)))&&strcmp(x,'2pole'));
p.addParameter('DS',0.1,@(x)isscalar(x)&&isnumeric(x)&&isfinite(x)&&x>0);
p.addParameter('RGLoops',4,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&& ...
    isfinite(x)&&x==fix(x)&&x>=1&&x<=4);
p.addParameter('Mu',2,@(x)isscalar(x)&&isnumeric(x)&&isfinite(x)&&x>0);
validAnomaly=@(x)isnumeric(x)&&isreal(x)&& ...
    (isempty(x)||(isscalar(x)&&isfinite(x)));
p.addParameter('Anomaly1S',-0.079,validAnomaly);
p.addParameter('Anomaly2S',[],validAnomaly);
p.addParameter('CheckSolution',false,@(x)isscalar(x)&&(islogical(x)||isnumeric(x)));
p.addParameter('RGInputs',struct('mbar_c',1.275,'alpha_mz',0.1184, ...
    'mz',91.1876,'mb',4.18),@charmonium.validRGInputs);
p.addParameter('Alpha',[],@(x)isempty(x)||(isnumeric(x)&&isvector(x)&&numel(x)>=3&& ...
    all(isfinite(x(:)))&&all(x(:)>0)&&all(diff(x(:))>0)));
p.parse(varargin{:}); o=p.Results;
o.AnomalySpecified=~all(ismember({'Anomaly1S','Anomaly2S'},p.UsingDefaults));
o.Channel=char(o.Channel);
if ~isempty(o.Alpha), o.Alpha=double(o.Alpha(:).'); end
end
