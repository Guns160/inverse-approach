function mq = psiMQTikhonov(psiQ,dMD)
%PSIMQTIKHONOV Monotone-queue psi_MQ for the Raus p0=1 rule.
psiQ=psiQ(:); dMD=dMD(:); n=numel(psiQ); mq=zeros(n,1); mq(n)=psiQ(n);
for j=n-1:-1:1
    mq(j)=max(psiQ(j),(dMD(j)/max(dMD(j+1),realmin))^(2/3)*mq(j+1));
end
end
