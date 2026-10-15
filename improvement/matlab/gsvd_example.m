%% gsvd_example —— 用一次 compact GSVD 取代逐 alpha 正规方程求解的示例（MATLAB 版）
%
% 对应仓库正式实现（同一套算术）：
%   +charmonium/compactGsvd.m      economy QR + built-in gsvd 的紧凑 GSVD
%   +charmonium/factorsOneWindow.m 端点消去 + Cholesky 罚项装配
%   +charmonium/solveInverse.m     单 alpha 与全 alpha 网格的滤波解
%
% 与本目录 ../python/gsvd_example.py 构造同一个合成反问题，输出可逐项对照。
% 运行：matlab -batch "run('<本文件绝对路径>')"，或在 MATLAB 中 cd 到本目录后运行。
% 只需基础 MATLAB；第 5 节的交叉验证需要 Hansen 的 Regularization Tools
% （cgsvd / fil_fac；未安装时自动跳过）。

%% ---------- 0) 被替换的原始解法（师兄 Pion_LCDA_test.m，解方程部分原文引用） ----------
% for i=1:n+mm
%     for j=1:n+mm
%        A2(i,j)=trapz(y,Aphi_New(:,i).*Aphi_New(:,j));
%     end
%     D2(i,1)=trapz(y,Aphi_New(:,i).*g_delta);
% end
% count_alphak=1000;
% alpha_range=logspace(log10(alpha_min),log10(alpha_max),count_alphak);
% for k=1:count_alphak
%     alpha2(k)=alpha_range(count_alphak-k+1);
%    T=(A2+alpha2(k)*H2);
%    Final= T\ [D2, e1, en];      % 每个 alpha 都重新做一次完整分解
%    ...
% end
% 更早版本（显式逆矩阵，精度更差）：
% FA(:,k)=inv(A2+alpha*H2)*D2;
% FA5(:,k)=inv(A2+alpha2(k)*H2)*(D2-Lag);
%
% 本示例把两端谱端点预先消去（正式程序 factorsOneWindow 的做法），因此基线
% 只需解 D2 一个右端项；师兄原代码同时解 [D2, e1, en]，用两个单位右端项处理
% 边界条件，每个 alpha 的代价只会更高。

%% ---------- 1) 合成反问题（结构同 factorsOneWindow.m） ----------
M2min=5; M2max=20; dM2=0.05;            % Borel 网格 t（数据点）
S_H=4; Lam=20; ds=0.1;                  % 谱 s 帽基网格
masses=[3.0969;3.6861];                 % 示例极点：J/psi、psi(2S)
kPoles=numel(masses);

t=(M2min:dM2:M2max).';  m2=numel(t);
w=trapwUniform(t);
s=(S_H:ds:Lam).';       nHats=numel(s);
H=massPlusStiffness(s);                % 帽基质量+刚度（师兄文件中的 M+M1）
boundary=[0;0];                        % 两端谱值固定为 0（示例取 0）

Poles=(exp(-masses.^2./t.')./t.').';    % 极点列，与 solveInverse.m 相同
Khat=exp(-(t*s.')./60)*ds;             % 光滑示意核（正式程序用解析 hatKernel）
K=[Poles,Khat];                        % m2 x (kPoles+nHats)

ni=nHats-2;
Hii=H(2:end-1,2:end-1); Hif=H(2:end-1,[1,end]);
Rb=chol(Hii,'upper');                  % 罚项行 = 内部刚度块的 Cholesky 因子
% 非零边界时的移位（本例 boundary=0，故 q=0、b0 不含移位项）：
% q=Hii\(Hif*boundary);
% b0=sqrt(w).*(g-K(:,fixed)*boundary)+A(:,kPoles+1:end)*q;
B=[zeros(ni,kPoles),Rb];               % 罚项行：极点不受罚
fixed=[kPoles+1,size(K,2)];            % 固定的两端帽列
free=[1:kPoles,kPoles+2:size(K,2)-1];  % 极点 + 内部帽
A=sqrt(w).*K(:,free);                  % 最小二乘矩阵（保持高矩阵，不平方）

xTrue=[1.2;0.45;0.8*exp(-(s(2:end-1)-12).^2/12)];
b0=A*xTrue;                            % 无噪声数据（重点在解法对比）
n=size(A,2); p=size(B,1);
assert(m2>=n && n>p, 'compact M3 GSVD 需要 m>=n>p');
assert(numel(free)==n && ni==p, '自由未知量结构不符');

%% ---------- 2) 基线：师兄的逐 alpha 正规方程解法 ----------
alpha=logspace(-16,-2,1001); nAlpha=numel(alpha);
A2=A'*A; D2=A'*b0;                     % 正规方程（条件数平方）
H2=zeros(n); H2(kPoles+1:end,kPoles+1:end)=Hii;  % 与 ||B*x||^2 完全相同的罚项
Xold=zeros(n,nAlpha);
% 旧法在小 alpha 处对接近奇异的正规方程逐点求解，MATLAB 自身也会告警
% （这本身就是精度问题的表现）；此处暂时屏蔽以免刷屏，也使计时公平。
warnState=warning;
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singular');
warning('off','MATLAB:illConditionedMatrix');
Xold(:,1)=(A2+alpha(1)*H2)\D2;         % 预热，使计时处于稳态
tic;
for j=1:nAlpha
    T=A2+alpha(j)*H2;                  % 每个 alpha 重新组装并完整分解
    Xold(:,j)=T\D2;
end
tOld=toc;
warning(warnState);                    % 恢复告警状态

%% ---------- 3) 改进：一次 compact GSVD + 全 alpha 网格滤波 ----------
[Qw,Rw]=qr(A,0); [Uw,Vw,Xw,Cw,Sw]=gsvd(Rw,B,'econ'); %#ok<*ASGLU>  % 预热
tic;
[Qa,Ra]=qr(A,0);                       % economy QR（同 compactGsvd.m）
[U,V,Xg,C,S]=gsvd(Ra,B,'econ');        % built-in gsvd
cn=sqrt(sum(C.^2,1)); sn=sqrt(sum(S.^2,1));
tol=max(1e-12,10*eps(max([1,norm(cn,inf),norm(sn,inf)])));
reg=sn>tol; unr=sn<=tol & cn>tol;
assert(nnz(reg)==p && nnz(unr)==n-p,'gsvd 秩结构与预期不符');
order=[find(reg),find(unr)];
U=U(:,order); Xg=Xg(:,order); C=C(order,order); S=S(:,order);
cn=cn(order); sn=sn(order);
sigma=cn(1:p).'; mu=sn(1:p).';         % 正则化模的广义奇异值
sigmaUn=cn(end-(n-p)+1:end).';         % 不受罚模（应全为 1）
assert(all(abs(sigmaUn-1)<1e-10),'不受罚模的 sigma 不为 1');
idA=norm(Ra-U*C*Xg','fro')/max(norm(Ra,'fro'),realmin);
idB=norm(B-V*S*Xg','fro')/max(norm(B,'fro'),realmin);
assert(idA<1e-12 && idB<1e-12,'gsvd 因子恒等式未通过');
tGsvd=toc;

betaAll=U'*(Qa'*b0);                   % 解系数空间右端项（U 对应 Ra，需先 Qa'*b0）
betaC=betaAll(1:p); betaN=betaAll(p+1:end);
recon=Qa*(U*betaAll);
null2=(b0-recon)'*(b0-recon);          % 与 alpha 无关的正交残差

tic;
den=sigma.^2+mu.^2.*alpha;            % p x nAlpha（隐式扩展）
coef=(sigma.*betaC)./den;             % Tikhonov 广义滤波系数
un=betaN./sigmaUn;
Y=[coef; repmat(un,1,nAlpha)];
Xnew=Xg'\Y;                           % 一次分解，整条网格同时回代
resAn=sqrt(sum((C*Y-betaAll).^2,1)+null2);  % 解析 ||A*x-b0||_w
penAn=sqrt(sum((S*Y).^2,1));          % 解析 ||B*x||
tFilter=toc;

%% ---------- 4) 验证与对比 ----------
idx=round(linspace(1,nAlpha,5));
fprintf('问题规模: m2=%d 数据点, n=%d 未知量(%d 极点 + %d 内部帽), p=%d 罚项行, alpha 网格 %d 点\n', ...
    m2,n,kPoles,ni,p,nAlpha);
fprintf('cond(A)=%.2e; 正规方程把条件数平方到≈%.1e，浮点下 A''A 已分辨不出小奇异方向（实测 cond(A''A)=%.2e）\n', ...
    cond(A),cond(A)^2,cond(A2));
fprintf('耗时: 旧法 %.1f ms; 新法 = GSVD %.1f ms + 全网格滤波 %.1f ms; 加速比 %.1f 倍\n', ...
    1e3*tOld,1e3*tGsvd,1e3*tFilter,tOld/(tGsvd+tFilter));
fprintf('\n  alpha        旧法相对误差   新法相对误差（相对增广最小二乘直接解）\n');
devR=0; devP=0;
for jj=idx
    a=alpha(jj);
    ref=[A;sqrt(a)*B]\[b0;zeros(p,1)]; % 直接解参考（同 solveInverse.m 的 CheckSolution）
    scaleRef=max(norm(ref),realmin);
    eOld=norm(Xold(:,jj)-ref)/scaleRef;
    eNew=norm(Xnew(:,jj)-ref)/scaleRef;
    rDirect=norm(A*Xnew(:,jj)-b0); pDirect=norm(B*Xnew(:,jj));
    devR=max(devR,abs(rDirect-resAn(jj))/max(rDirect,realmin));
    devP=max(devP,abs(pDirect-penAn(jj))/max(pDirect,realmin));
    fprintf('%8.1e   %10.3e   %10.3e\n',a,eOld,eNew);
end
relGrid=vecnorm(Xnew-Xold,2,1)./max(vecnorm(Xnew,2,1),realmin);
[relMax,iMax]=max(relGrid);
fprintf('解析范数自检: 残差范数最大相对偏差 %.2e, 罚项范数最大相对偏差 %.2e\n',devR,devP);
fprintf('全网格新旧解最大相对差 %.2e（出现在 alpha=%.1e；小 alpha 处旧法失真）\n',relMax,alpha(iMax));

%% ---------- 5) 可选交叉验证：Hansen Regularization Tools（cgsvd + fil_fac） ----------
if exist('cgsvd','file')==2 && exist('fil_fac','file')==2
    [Uh,sm,Xh]=cgsvd(A,B);             % 文献标准紧凑 GSVD（直接对高矩阵 A 分解）
    sh=sm(:,1); mh=sm(:,2);
    stOurs=sort(sigma./mu);
    relSv=norm(sort(sh./mh)-stOurs,inf)/max(norm(stOurs,inf),realmin);
    % 该包 ''Tikh'' 滤波用 lambda^2 约定：f=sh^2/(sh^2+lambda^2*mh^2)；
    % 本示例与仓库实现的 alpha 对应 lambda^2。
    bh=Uh'*b0; maxRelX=0;
    for jj=idx
        f=fil_fac(sm,sqrt(alpha(jj)),'Tikh');
        yh=[f.*(bh(1:p)./sh); bh(p+1:n)];
        xh=Xh*yh;
        maxRelX=max(maxRelX,norm(xh-Xnew(:,jj))/max(norm(Xnew(:,jj)),realmin));
    end
    fprintf('cgsvd 交叉验证: 广义奇异值相对差 %.2e, 解最大相对差 %.2e\n',relSv,maxRelX);
else
    disp('cgsvd 交叉验证: 未检测到 Regularization Tools，跳过');
end

%% ---------- 局部函数 ----------
function w=trapwUniform(t)
%TRAPWUNIFORM 均匀网格梯形权（同 charmonium.trapweights 的均匀情形）。
dt=t(2)-t(1); w=dt*ones(size(t)); w(1)=dt/2; w(end)=dt/2;
end

function H=massPlusStiffness(s)
%MASSPLUSSTIFFNESS 帽基质量 Gram（M）+ 一阶差分刚度（M1），H=M+M1。
% 构造与师兄 Pion_LCDA_test.m 中的 M、M1、H1 相同（h 为帽基步长）。
h=s(2)-s(1); nH=numel(s);
M=zeros(nH); M1=zeros(nH);
M(1,1)=4; M(end,end)=4; M(1,2)=2; M(end,end-1)=2; M(2,1)=2; M(end-1,end)=2;
M1(1,1)=1; M1(end,end)=1; M1(1,2)=-1; M1(end,end-1)=-1; M1(2,1)=-1; M1(end-1,end)=-1;
for j=2:nH-1
    M(j,j)=8; M(j,j-1)=2; M(j,j+1)=2;
    M1(j,j)=2; M1(j,j-1)=-1; M1(j,j+1)=-1;
end
H=M*(h/12)+M1*(1/h);
end
