function [lo,nlo] = ppRhoComponents(s,mc,alphaS)
%PPRHOCOMPONENTS LO and NLO P-P perturbative spectral densities.
values=double(s); lo=zeros(size(values)); nlo=lo;
threshold=4*mc^2; above=values>threshold; at=values==threshold;
sa=values(above); v=sqrt(1-threshold./sa);
lo(above)=3*sa.*v/(8*pi^2);
nlo(above)=sa.*alphaS.*charmonium.ppFiniteFunction(v)/(2*pi^3);
nlo(at)=threshold*alphaS/(4*pi);
end
