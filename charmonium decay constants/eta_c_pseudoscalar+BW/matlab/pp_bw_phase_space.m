function phi = pp_bw_phase_space(s,model)
%PP_BW_PHASE_SPACE Charge-averaged phase spaces with Blatt-Weisskopf factors.
% Broad daughters use GL quadrature in u with a sin^2 map and the physical
% D(*)pi lower limit. The barrier is evaluated INSIDE the u integral.
validateattributes(s,{'numeric'},{'real','finite','positive','vector'});
s=double(s(:)); R=model.radius; D=model.D; V=model.Dstar;
if strcmp(model.channel,'6S')
    error('charmonium:BW6SUnavailable','6S running BW partial widths are unavailable');
end
switch model.channel
    case '3S'
        phi=averageStable(s,D,V,1,R);
    case '4S'
        phi=[averageStable(s,D,V,1,R),stable(s,model.Ds,model.Dsstar,1,R), ...
            averageStable(s,V,V,1,R),stable(s,model.Dsstar,model.Dsstar,1,R)];
    case '5S'
        phi=[averageStable(s,D,V,1,R),stable(s,model.Ds,model.Dsstar,1,R), ...
            broad(s,D,model.D0star,model.D0_width,D,0,model), ...
            stable(s,model.Ds,model.Ds0star,0,R),averageStable(s,D,model.D2star,2,R), ...
            averageStable(s,V,V,1,R),stable(s,model.Dsstar,model.Dsstar,1,R), ...
            broad(s,V,model.D1_2420,model.D1_width,V,model.DstarD1L,model)];
    otherwise
        error('charmonium:BWChannel','Unknown BW channel');
end
end

function result = averageStable(s,a,b,L,R)
result=(stable(s,a(1),b(1),L,R)+stable(s,a(2),b(2),L,R))/2;
end

function result = stable(s,a,b,L,R)
q=momentum(s,a,b^2);
result=q.^(2*L+1)./sqrt(s).*barrier(q,L,R);
end

function q = momentum(s,a,u)
b=sqrt(u);
q=sqrt(max((s-(a+b).^2).*(s-(a-b).^2),0))./(2*sqrt(s));
q=q.*(sqrt(s)>=a+b);
end

function result = barrier(q,L,R)
z2=(q*R).^2;
switch L
    case 0, result=ones(size(q));
    case 1, result=1./(1+z2);
    case 2, result=1./(9+3*z2+z2.^2);
    otherwise, error('charmonium:BWAngularMomentum','Only L=0,1,2 are defined');
end
end

function result = broad(s,spectators,poles,width,daughters,L,model)
result=zeros(size(s));
for k=1:2
    charged=(daughters(3-k)+model.pion(2))^2;
    neutral=(daughters(k)+model.pion(1))^2;
    result=result+0.5*((2/3)*convolution(s,spectators(k),poles(k),width,charged,L,model)+ ...
        (1/3)*convolution(s,spectators(k),poles(k),width,neutral,L,model));
end
end

function result = convolution(s,a,b,width,lower,L,model)
result=zeros(size(s)); ids=find(sqrt(s)>a+sqrt(lower));
[x,w]=legendreNodes(model.quadrature_nodes);
theta=(x.'+1)*pi/4; sin2=sin(theta).^2; jacobian=sin(2*theta)*pi/4;
% Bound temporary memory independently of the Borel/spectral grid length.
for first=1:128:numel(ids)
    use=ids(first:min(first+127,numel(ids))); ss=s(use);
    span=(sqrt(ss)-a).^2-lower;
    u=lower+span.*sin2; jac=span.*jacobian;
    rho=b*width./(pi*((u-b*b).^2+(b*width)^2));
    q=momentum(ss,a,u);
    integrand=rho.*q.^(2*L+1)./sqrt(ss).*barrier(q,L,model.radius).*jac;
    result(use)=integrand*w;
end
end

function [x,w] = legendreNodes(n)
% Newton-Legendre roots avoid a dense 1024-by-1024 eigendecomposition.
persistent cachedN cachedX cachedW
if isequal(cachedN,n), x=cachedX; w=cachedW; return; end
x=cos(pi*((1:n).'-0.25)/(n+0.5));
for iteration=1:30
    [pn,pnm1]=legendrePair(x,n);
    derivative=n*(x.*pn-pnm1)./(x.^2-1);
    correction=pn./derivative; x=x-correction;
    if max(abs(correction))<4*eps, break; end
end
assert(iteration<30,'Gauss-Legendre nodes did not converge');
[pn,pnm1]=legendrePair(x,n);
derivative=n*(x.*pn-pnm1)./(x.^2-1);
w=2./((1-x.^2).*derivative.^2);
[x,order]=sort(x); w=w(order);
cachedN=n; cachedX=x; cachedW=w;
end

function [pn,pnm1] = legendrePair(x,n)
pnm1=ones(size(x)); pn=x;
for k=2:n
    next=((2*k-1)*x.*pn-(k-1)*pnm1)/k;
    pnm1=pn; pn=next;
end
end
