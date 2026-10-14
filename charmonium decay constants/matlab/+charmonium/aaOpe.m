function [total,pert,g2,g3] = aaOpe(Lambda,m2,mc,alphaS,hasNLO,b3Cross)
%AAOPE Axial/vector n=0 OPE at supplied RG inputs.
m2=double(m2(:)); q=4*mc^2;
v=linspace(0,sqrt(max(1-q/Lambda,0)),16001).';
sg=q./max(1-v.^2,q/Lambda); sg(end)=Lambda;
rw=charmonium.aaRhoPert(sg,mc,alphaS,hasNLO).*charmonium.trapweights(sg);
pert=zeros(size(m2));
for first=1:256:numel(m2)
    ids=first:min(first+255,numel(m2)); mm=m2(ids);
    pert(ids)=(exp(-sg.'./mm)./mm)*rw;
end
x=linspace(1e-8,1-1e-8,4001); xb=1-x; wx=charmonium.trapweights(x.');
a3=-(45/8)*x.*xb.*(x.^3+xb.^3)-(69/72)*x.^2.*xb.^2;
b3=-(3/4)*mc^2*x.*xb.*(x.^4+xb.^4)+ ...
    mc^2*x.^2.*xb.^2.*(-(23/12)*(x.^2+xb.^2)+b3Cross*x.*xb+2);
c3=(2/5)*mc^4*x.*xb.*(x.^5+xb.^5);
g2=zeros(size(m2)); g3=g2;
for first=1:128:numel(m2)
    ids=first:min(first+127,numel(m2)); mm=m2(ids);
    E=exp(-mc^2./(x.*xb.*mm));
    i2=E.*(1./(2*mm.^2)-mc^2*(x.^3+xb.^3)./(2*mm.^3.*x.^2.*xb.^2));
    i3=E.*(a3./(2*mm.^3.*x.^3.*xb.^3)+ ...
        b3./(6*mm.^4.*x.^4.*xb.^4)+c3./(24*mm.^5.*x.^5.*xb.^5));
    g2(ids)=0.038*(i2*wx)/(6*pi);
    g3(ids)=0.013*(i3*wx)/(4*pi)^2;
end
total=pert+g2+g3;
end
