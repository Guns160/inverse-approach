function sel = selectRaus(sigma,beta,null2,alpha)
%SELECTRAUS Raus 2024 r_H heuristic, using alpha rather than alpha^2.
a=alpha(:); [d,q]=charmonium.rausQuantities(sigma,beta,null2,a);
mq=charmonium.psiMQTikhonov(q,d); psiHR=d./sqrt(a);
inHQ=a>=max(a(1),min(sigma(:)).^2);
hrValues=inf(size(psiHR)); hrValues(inHQ)=psiHR(inHQ); [~,iHR]=min(hrValues);
mins=charmonium.localMinima(q);
if isempty(mins)
    hrMins=charmonium.localMinima(psiHR);
    if ~isempty(hrMins), iH=hrMins(charmonium.argmin(psiHR(hrMins))); fallback='psiHR_local_min';
    else, iH=iHR; fallback='psiHR_global'; end
    iQ=iH; iMQ=iH; source=fallback;
else
    iQ=mins(charmonium.argmin(q(mins))); iMQ=mins(charmonium.argmin(mq(mins)));
    rH=min(1./a(iQ),max(1./a(iMQ),1./a(iHR)));
    [~,iH]=min(abs(a-1./rH));
    if iH==iQ, source='q'; elseif iH==iMQ, source='mq';
    elseif iH==iHR, source='hr'; else, source='nearest'; end
    fallback='';
end
sel=struct('index',iH,'index_q_local',iQ,'index_mq',iMQ,'index_hr',iHR, ...
    'source',source,'fallback',fallback,'n_local_minima',numel(mins), ...
    'status','defined','alpha',a(iH),'alpha_Q',a(iQ),'alpha_MQ',a(iMQ), ...
    'alpha_HR',a(iHR),'n_alpha',numel(a),'d_MD',d,'psi_Q',q,'psi_MQ',mq, ...
    'psi_HR',psiHR);
end
