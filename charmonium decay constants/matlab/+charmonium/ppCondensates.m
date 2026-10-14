function [g2,g3] = ppCondensates(m2,mc)
%PPCONDENSATES P-P D4 and D6 Borel contributions at the supplied charm mass.
m2=double(m2(:)); x=linspace(1e-8,1-1e-8,4001); xb=1-x;
q=x.*xb; wx=charmonium.trapweights(x.');
cubic=mc^2*(x.^3+xb.^3); a=mc*x+mc*xb; x3=x.^3+xb.^3;
b4=64-(416/3)*q-45*x3./q;
b6=((45*a.^2./q-(2/3)*mc^2).*x3-6*mc^2*(x.^4+xb.^4)./q- ...
    12*mc^2*x3-(478/9)*mc^2*q+(287/9)*a.^2-(320/9)*q*mc^2);
b8=((6*a.^2./q+28*mc^2).*mc^2.*(x.^4+xb.^4)-4*mc^4.*(x.^4+xb.^4)+ ...
    (16/5)*mc^4*(x.^5+xb.^5)./q+(30*a.^2+88*mc^2*q).*mc^2.*(x.^2+xb.^2)- ...
    (112/3)*a.^4+(64/3)*mc^4*q+128*mc^4*q.^2);
g2=zeros(size(m2)); g3=g2;
for first=1:128:numel(m2)
    ids=first:min(first+127,numel(m2)); mm=m2(ids);
    expo=exp(-mc^2./(q.*mm));
    term4=(-6*mc^2-2*cubic./q)./(mm.^2.*x.^2.*xb.^2);
    term6=a.^2.*cubic./(mm.^3.*x.^4.*xb.^4);
    i2=expo.*(term4+term6);
    i3=expo.*(b4./(mm.^2.*x.^2.*xb.^2)+ ...
        b6./(2*mm.^3.*x.^3.*xb.^3)+b8./(6*mm.^4.*x.^4.*xb.^4)- ...
        (2/15)*a.^2.*mc^4.*(x.^5+xb.^5)./(mm.^5.*x.^6.*xb.^6));
    g2(ids)=0.038*(i2*wx)/(48*pi);
    g3(ids)=0.013*(i3*wx)/(32*(4*pi)^2);
end
end
