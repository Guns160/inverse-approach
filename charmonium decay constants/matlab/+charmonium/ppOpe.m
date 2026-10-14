function [total,lo,nlo,g2,g3] = ppOpe(Lambda,m2,mc,alphaS)
%PPOPE P-P n=0 LO+NLO+D4+D6 OPE at supplied RG inputs.
m2=double(m2(:)); threshold=4*mc^2;
v=linspace(0,sqrt(1-threshold/Lambda),16001).';
s=threshold./(1-v.^2); ws=charmonium.trapweights(s);
[loRho,nloRho]=charmonium.ppRhoComponents(s,mc,alphaS);
lo=zeros(size(m2)); nlo=lo;
for first=1:256:numel(m2)
    ids=first:min(first+255,numel(m2)); mm=m2(ids);
    weight=exp(-s.'./mm)./mm;
    lo(ids)=weight*(loRho.*ws); nlo(ids)=weight*(nloRho.*ws);
end
[g2,g3]=charmonium.ppCondensates(m2,mc);
total=lo+nlo+g2+g3;
end
