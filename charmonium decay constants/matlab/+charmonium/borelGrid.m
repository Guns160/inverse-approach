function grid = borelGrid(lo, hi)
%BORELGRID Fixed 0.01 GeV^2 Borel grid with exact requested endpoints.
step = 0.01;
index = round((lo-1.5)/step);
anchor = 1.5+step*index;
tol = 8*eps(max([1,abs(lo),abs(hi)]));
n = floor((hi-lo)/step);
if abs(anchor-lo) <= tol
    grid = 1.5+step*(index+(0:n).');
else
    grid = lo+step*(0:n).';
end
grid(1) = lo;
grid = grid(grid < hi);
if numel(grid) > 1 && hi-grid(end) <= tol
    grid(end) = hi;
else
    grid(end+1) = hi;
end
assert(all(diff(grid) > 0), 'Borel grid must be strictly increasing');
end
