# improvement：解方程环节的改进示例（compact GSVD 滤波 vs 逐 alpha 正规方程）

这个文件夹只包含**解方程环节**的最小可运行示例，用来说明本仓库对原始解法的改进。两个示例（Python + MATLAB）构造**同一个**合成反问题，输出可逐项对照。

## 原方法与改进方法

| | 原方法 | 改进方法（本仓库三套程序共用） |
|---|---|---|
| 出处 | 师兄 `Pion_LCDA_test.m`（解方程部分原文引用在 `matlab/gsvd_example.m` 头部注释） | [`../charmonium decay constants/`](../charmonium%20decay%20constants/) 的公共求解核心 |
| 装配 | 正规方程 `A2=A'WA`、`D2=A'Wb`、罚项 `H2` | 保持高矩阵 `A=sqrt(w)*K_free`，罚项行 `B=[0, chol(H_ii)]`，谱端点消去 |
| 求解 | 每个 alpha 重新组装并完整分解 `(A2+alpha*H2)`（更早版本用 `inv(A2+alpha*H2)*D2`） | **一次** economy-QR + compact GSVD，之后任意 alpha 的解都是标量滤波 |
| 整条 alpha 网格 | 逐点求解，nAlpha 次完整分解（L-curve / 准最优准则需要对整条网格求解） | 一次分解 + 向量化滤波，`X'\[coef; un]` 一次回代整条网格 |
| 诊断量 | 每点重新回代算残差 / 罚项 | 残差范数 `norm(C*y-beta)^2+null`、罚项范数 `norm(S*y)` 有解析式 |

滤波公式（与 `solveInverse.m`、`reconstruction.py` 相同）：

$$
y_i(\alpha)=\frac{\sigma_i\,\beta_i}{\sigma_i^2+\alpha\,\mu_i^2}\ \ (i\le p),
\qquad
y_i=\frac{\beta_i}{\sigma_i}\ \ (i>p,\ \sigma_i=1),
\qquad
x_\alpha=X\,y,
$$

其中 `(sigma_i, mu_i)` 是矩阵对 `(A, B)` 的广义奇异值，`beta=U'*(Q'*b)`。旧法求解的正规方程 `(A'A+alpha*B'B)x=A'b` 与它对应**同一个** Tikhonov 泛函；两个示例中只有解法不同，对比公平。

对应正式实现：

- Python：`../charmonium decay constants/inverse_common/gsvd.py`（LAPACKE `dggsvd3`）、`factors.py`（装配）、`reconstruction.py`（全 alpha 网格滤波）
- MATLAB：`../charmonium decay constants/matlab/+charmonium/compactGsvd.m`、`factorsOneWindow.m`、`solveInverse.m`

## 为什么快

- 原方法选 alpha 需要对整条网格逐点求解：师兄代码 1000 个点就是 1000 次 O(n^3) 完整分解（原代码每次还同时解 `[D2, e1, en]` 三个右端项）。
- 改进方法只做一次 GSVD；之后每个 alpha 的解是 O(p) 标量滤波，整条网格一次算完。示例规模（n=161、1001 个 alpha）实测滤波部分仅 1.8–3.7 ms，对应旧法逐点求解 103–140 ms。

## 为什么准

- 正规方程把条件数平方：示例中 `cond(A)≈1e+18`，理论上 `cond(A'A)≈cond(A)^2≈1e+36`，`A'A` 的最小奇异方向在浮点下已被舍入淹没。实测小 alpha 处旧法相对误差 30–40%，新法保持在 1e-9 量级（见下表）。
- 不使用显式逆矩阵；GSVD 滤波解与增广最小二乘直接解 `[A; sqrt(alpha)*B]\[b; 0]` 一致到机器精度量级（这正是正式程序 `solveInverse.m` 的 `CheckSolution` 校验）。
- 残差 / 罚项范数有解析式；极小 alpha 处比直接回算更稳（示例最小 alpha 处直接回算有 4e-7–3e-6 量级的抵消误差，解析式无此问题）。

## 实测（2026-10-15）

环境：Python 3.9.6 + NumPy 2.0.2（GSVD 调 SciPy 自带 OpenBLAS 的 LP64 `LAPACKE_dggsvd3`）；MATLAB R2026a（built-in `gsvd`）。问题规模：301 个 Borel 数据点，161 个未知量（2 个不受罚极点 + 159 个内部帽基），159 条罚项行，alpha 网格 1001 点（logspace 1e-16…1e-2）。误差为相对增广最小二乘直接解的相对误差：

| alpha | 旧法（Python / MATLAB） | 新法（Python / MATLAB） |
|---|---|---|
| 1e-16 | 3.9e-01 / 3.0e-01 | 8.3e-09 / 2.3e-09 |
| 3.2e-13 | 1.9e-04 / 1.0e-04 | 2.0e-10 / 5.8e-11 |
| 1e-09 | 6.3e-08 / 2.5e-08 | 1.3e-12 / 3.1e-12 |
| 3.2e-06 | 5.8e-10 / 2.0e-10 | 3.1e-13 / 4.7e-13 |
| 1e-02 | 2.1e-12 / 3.9e-12 | 2.3e-14 / 9.3e-15 |

耗时（1001 个 alpha 的整条网格）：

| | 旧法（逐点求解） | 新法（一次 GSVD + 全网格滤波） | 加速比 |
|---|---|---|---|
| Python（LAPACKE dggsvd3） | 103.3 ms | 56.1 ms + 1.8 ms | 1.8 倍 |
| MATLAB（built-in gsvd） | 140.2 ms | 5.0 ms + 3.7 ms | 16.2 倍 |

说明：新法的分解是一次性成本，alpha 网格越长、问题规模越大，摊销越明显（旧法代价 ∝ nAlpha*n^3，新法滤波 ∝ n^2*nAlpha + 一次 GSVD）。示例为演示规模；正式程序默认 alpha 网格 2001 点，历史生产谱网格（ds=0.01）更大。

其它自检：

- 罚项范数解析式 vs 直接回算：Python 6.0e-15、MATLAB 7.9e-16。
- MATLAB 端用 Hansen Regularization Tools 的 `cgsvd`+`fil_fac` 独立复核：广义奇异值（sigma/mu）相对差 0.0，解最大相对差 3.8e-10。
- 两语言构造同一问题（网格值仅差最后几位舍入），各输出项一致到舍入 / 量级水平。

## 运行

Python（需 NumPy；GSVD 需要 LP64 `LAPACKE_dggsvd3`，自动探测系统 / NumPy / SciPy 的 OpenBLAS，也可用环境变量 `JPSI_LAPACKE_LIBRARY` 指定。SciPy 没有暴露 gsvd 接口，示例与正式实现一样直接用 ctypes 调 LAPACKE）：

```sh
python3 python/gsvd_example.py
```

MATLAB（只需基础 MATLAB；末段 cgsvd 交叉验证需要 Regularization Tools，未安装则自动跳过）：

```matlab
matlab -batch "run('matlab/gsvd_example.m')"   % 或在 MATLAB 中 cd 到 matlab/ 后运行 gsvd_example
```

## 与 Regularization Tools 的关系

紧凑 GSVD 的排序约定沿用 Hansen《Regularization Tools》的 `cgsvd`（正式实现 `gsvd.py` 的说明也注明了这一点）。约定差异：该包 `'Tikh'` 滤波因子写作 `sigma^2/(sigma^2+lambda^2*mu^2)`，本仓库的 alpha 对应 `lambda^2`。

## 边界

示例用光滑示意核（`exp(-s*t/60)`）构造合成数据，只演示解方程环节；正式物理管道（解析 hatKernel、非零边界 q 移位、Raus / 准最优 alpha 选择、三通道物理与 RG 约定）见 [`../charmonium decay constants/`](../charmonium%20decay%20constants/)。
