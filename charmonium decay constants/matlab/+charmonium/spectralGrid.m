function s = spectralGrid(Lambda, sh, ds)
%SPECTRALGRID Physical-s hat anchors, retaining the exact Lambda endpoint.
assert(isfinite(Lambda) && Lambda > sh, 'Lambda must be finite and exceed S_H');
assert(isscalar(ds) && isfinite(ds) && ds > 0, 'DS must be a finite positive scalar');
n = floor((Lambda-sh)/ds);
s = sh+ds*(0:n).';
s = s(s < Lambda);
tol = 8*eps(max([1,abs(sh),abs(Lambda)]));
if numel(s) > 1 && Lambda-s(end) <= tol
    s(end) = Lambda;
else
    s(end+1) = Lambda;
end
assert(numel(s) >= 3, 'Need at least three hat nodes');
assert(all(diff(s) > 0), 'Hat grid must be strictly increasing');
end
