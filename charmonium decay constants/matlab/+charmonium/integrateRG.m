function [alphaEnd, massEnd, info] = integrateRG(alphaStart, massStart, muStart, muEnd, nf, loops)
%INTEGRATERG Adaptive fixed-nf ODE, state [alpha_s/pi; ln(mass)].
% Inputs and outputs are alpha_s and mass in GeV; scales are in GeV.
% Returns [alpha_s(end), mass(end)], matching Python evolve_fixed_nf.
% Holds nf fixed; this low-level routine performs no threshold inference.
if nargin < 5, nf = 4; end
if nargin < 6, loops = 4; end
validateattributes(alphaStart, {'numeric'}, {'real','scalar','finite','positive','<',pi});
validateattributes(massStart, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(muStart, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(muEnd, {'numeric'}, {'real','scalar','finite','positive'});
c = charmonium.qcdCoefficients(nf,loops);
alphaStart = double(alphaStart); massStart = double(massStart);
muStart = double(muStart); muEnd = double(muEnd);
loops = double(loops);
rtol = 2e-12; atol = 2e-14;
info = struct('method','ode113','rtol',rtol,'atol',atol,'mass_method','ode');
if muStart == muEnd
    alphaEnd = alphaStart; massEnd = massStart; return
end
t0 = 2*log(muStart); t1 = 2*log(muEnd);
options = odeset('RelTol',rtol,'AbsTol',atol,'Events',@guard);
try
    [t,y,te] = ode113(@rhs,[t0,t1],[alphaStart/pi;log(massStart)],options);
catch cause
    error('charmonium:RGIntegration','QCD RG integration failed: %s',cause.message);
end
if ~isempty(te)
    error('charmonium:RGIntegration','QCD RG reached perturbative guard alpha_s=pi');
end
if isempty(t) || abs(t(end)-t1) > 8*eps(max(1,abs(t1)))
    error('charmonium:RGIntegration','QCD RG integration did not reach target scale');
end
aEnd = y(end,1); logMass = y(end,2);
alphaEnd = pi*aEnd; massEnd = exp(logMass);
if ~all(isfinite([alphaEnd,massEnd])) || alphaEnd <= 0 || alphaEnd >= pi || massEnd <= 0
    error('charmonium:RGIntegration','QCD RG returned an invalid/nonperturbative final state');
end
    function dy = rhs(~,state)
        a = state(1); b = c.beta(end); g = c.gamma(end);
        for j = loops-1:-1:1
            b = c.beta(j)+a*b; g = c.gamma(j)+a*g;
        end
        dy = [-a*a*b; -a*g];
    end
    function [value,isterminal,direction] = guard(~,state)
        value = 1-state(1); isterminal = 1; direction = 0;
    end
end
