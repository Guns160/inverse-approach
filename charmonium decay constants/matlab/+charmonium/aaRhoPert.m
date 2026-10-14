function rho = aaRhoPert(s,mc,alphaS,hasNLO)
%AARHOPERT LO (and optional NLO factor) perturbative spectral density.
s=double(s); rho=zeros(size(s)); mask=s>4*mc^2;
v=sqrt(1-4*mc^2./s(mask));
rho(mask)=3*v.*(1-v.^2/3)/(8*pi^2);
if hasNLO, rho(mask)=rho(mask)*(1+alphaS/pi); end
end
