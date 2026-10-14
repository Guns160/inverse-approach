function alphaLow = matchAlphaDownMS(alphaHigh, nfHigh, matchingLoops)
%MATCHALPHADOWNMS MSbar matching ONLY at mu=mhbar(mhbar), hence L=0.
% RunDec Eqs. (19)-(21); alphaLow=alphaHigh*zeta_g^2.
% Arguments use alpha_s (NOT a) and actual matching order (NOT running order).
if nargin < 2, nfHigh = 5; end
if nargin < 3, matchingLoops = 3; end
validateattributes(alphaHigh, {'numeric'}, {'real','scalar','finite','positive','<',pi});
validateattributes(nfHigh, {'numeric'}, {'real','scalar','finite','integer','>=',2,'<=',6});
validateattributes(matchingLoops, {'numeric'}, {'real','scalar','finite','integer','>=',0,'<=',3});
a = double(alphaHigh)/pi; nl = double(nfHigh)-1;
d2 = 11/72;
d3 = 564731/124416-82043*1.2020569031595942854/27648-2633*nl/31104;
factor = 1;
if matchingLoops >= 2, factor = factor+d2*a^2; end
if matchingLoops >= 3, factor = factor+d3*a^3; end
alphaLow = double(alphaHigh)*factor;
if ~isfinite(alphaLow) || alphaLow <= 0 || alphaLow >= pi
    error('charmonium:RGMatching','Matching returned invalid/nonperturbative alpha_s');
end
end
