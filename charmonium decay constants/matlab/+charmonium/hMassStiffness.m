function H = hMassStiffness(s)
%HMASSSTIFFNESS H1 penalty H=M+M1 on the physical-s hat grid.
s = double(s(:)); H = zeros(numel(s));
for j = 1:numel(s)-1
    h = s(j+1)-s(j);
    H(j:j+1,j:j+1) = H(j:j+1,j:j+1)+ ...
        [h/3+1/h,h/6-1/h;h/6-1/h,h/3+1/h];
end
end
