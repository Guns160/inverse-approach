function [dMD,psiQ] = rausQuantities(sigma,beta,null2,alpha)
%RAUSQUANTITIES d_MD and psi_Q on an alpha grid (not alpha squared).
s=sigma(:); b=beta(:); a=alpha(:).';
den=s.^2+a;
dMD=sqrt(sum(((a.^1.5.*b)./den.^1.5).^2,1)+null2);
psiQ=sqrt(sum(((a.*s.*b)./den.^2).^2,1));
end
