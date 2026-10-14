function ok = validRGInputs(x)
%VALIDRGINPUTS Require the complete, coupled MSbar anchor set.
fields=["mbar_c","alpha_mz","mz","mb"];
ok=isstruct(x)&&isscalar(x)&&numfields(x)==numel(fields)&&all(isfield(x,fields));
if ok
    values=[x.mbar_c,x.alpha_mz,x.mz,x.mb];
    ok=all(isfinite(values))&&all(values>0)&&x.mbar_c<x.mb&&x.mb<x.mz;
end
end
