function w = trapweights(x)
%TRAPWEIGHTS Trapezoidal quadrature weights on an increasing grid.
x = double(x(:));
w = [diff(x(1:2))/2; (diff(x(1:end-1))+diff(x(2:end)))/2; diff(x(end-1:end))/2];
end
