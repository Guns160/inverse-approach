function f = compactGsvd(A,B)
%COMPACTGSVD Economy-QR followed by MATLAB's built-in gsvd.
% Convention: R=U*C*X', L=V*S*X'; solution coefficients live in y=X'*x.
A=double(A); B=double(B); [m,n]=size(A); p=size(B,1);
assert(size(B,2)==n && m>=n && p<n, 'Compact M3 GSVD requires m>=n>p');
[Q,R]=qr(A,0);
[U,V,X,C,S]=gsvd(R,B,'econ');
cn=sqrt(sum(C.^2,1)); sn=sqrt(sum(S.^2,1));
tol=max(1e-12,10*eps(max([1,norm(cn,inf),norm(sn,inf)])));
regularized=sn>tol; unregularized=sn<=tol & cn>tol;
assert(nnz(regularized)==p && nnz(unregularized)==n-p, ...
    'Unexpected MATLAB gsvd rank structure');
order=[find(regularized),find(unregularized)];
U=U(:,order); X=X(:,order); C=C(order,order); S=S(:,order);
cn=cn(order); sn=sn(order);
identityA=norm(R-U*C*X','fro')/max(norm(R,'fro'),realmin);
identityB=norm(B-V*S*X','fro')/max(norm(B,'fro'),realmin);
assert(identityA<1e-12 && identityB<1e-12, 'MATLAB gsvd factor identities failed');
assert(all(abs(cn(1:p).^2+sn(1:p).^2-1)<1e-10), ...
    'Regularized gsvd columns are not normalized');
assert(all(abs(cn(end-(n-p)+1:end)-1)<1e-10), ...
    'Unregularized gsvd columns are not data modes');
cc=cn.'; mu=sn.';
assert(all(isfinite([cc;mu])), 'gsvd returned non-finite factors');
assert(all(cc(1:p)>0) && all(mu(1:p)>realmin), 'gsvd returned a zero compact mode');
f=struct('Q',Q,'R',R,'U',U,'X',X,'C',C,'S',S,'sigma',cc,'mu',mu, ...
    'sigma_unregularized',cc(end-(n-p)+1:end),'beta',U'*0, ...
    'source','MATLAB gsvd after economy QR','library','MATLAB R2026a built-in gsvd', ...
    'column_order',order,'identity_relative_error_A',identityA, ...
    'identity_relative_error_B',identityB,'A',Q*R,'B',B);
end
