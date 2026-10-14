function K = hatKernel(m2, s)
%HATKERNEL Analytic Borel transform of piecewise-linear physical-s hats.
m2 = double(m2(:)); s = double(s(:));
K = zeros(numel(m2),numel(s));
for j = 1:numel(s)-1
    h = s(j+1)-s(j); z = h./m2;
    small = abs(z) < 1e-4;
    asc = zeros(size(z)); des = asc;
    t = z(small);
    asc(small) = 0.5-t/3+t.^2/8-t.^3/30+t.^4/144;
    des(small) = 0.5-t/6+t.^2/24-t.^3/120+t.^4/720;
    t = z(~small); e = exp(-t);
    asc(~small) = (1-(1+t).*e)./t.^2;
    des(~small) = (t-1+e)./t.^2;
    fac = exp(-s(j)./m2).*z;
    K(:,j) = K(:,j)+fac.*des;
    K(:,j+1) = K(:,j+1)+fac.*asc;
end
end
