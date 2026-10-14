%SELECTIONGRID Dynamic AA alpha grid copied from the frozen Python rule.
function [alpha,covered] = selectionGrid(sigma,beta,null2)
sigMin=max(min(sigma(:)),1e-300);
logMin=min(-50,2*log10(sigMin)-8); logMin=max(logMin,-80); last=[];
for iter=1:10
    cand=10.^linspace(logMin,0,round((0-logMin)/0.01)+1);
    [~,psi]=charmonium.rausQuantities(sigma,beta,null2,cand);
    peaks=charmonium.localMaxima(psi);
    hasPeak=~isempty(peaks)&&peaks(1)>=4&&psi(1)<psi(peaks(1));
    last=cand;
    if hasPeak, covered=true; break; end
    covered=false;
    if logMin<=-80, break; end
    logMin=max(-80,logMin-5);
end
alpha=10.^linspace(log10(last(1)),0,numel(last));
[~,psi]=charmonium.rausQuantities(sigma,beta,null2,alpha);
peaks=charmonium.localMaxima(psi);
covered=~isempty(peaks)&&peaks(1)>=4&&psi(1)<psi(peaks(1));
end
