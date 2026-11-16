# 验证记录：2026-10-14

**中文** | [English](VALIDATION.md)

这些检查验证实现及选定的数值契约。它们不确立新的物理提取窗口，也不替代现有误差预算分析。

## 环境

- Python 3.9.6；NumPy 2.0.2；SciPy 1.13.1。
- Python GSVD：来自现有 SciPy OpenBLAS 库的 LP64 `LAPACKE_dggsvd3`。
- MATLAB R2026a Update 4，版本 26.1.0.3312084，内置 `gsvd`。
- macOS 主机。本次运行未在其他操作系统上验证相同的 Python 代码和库发现机制。

## 分开的变更

1. 原始源码/布局：仓库提交 `06fb9a2`。
2. 共享 Python 基础设施：仓库提交 `bfe62ea`。保留旧输入常数。在 `(M2_min,M2_max,Lambda)=(5,20,20)`、`ds=0.1` 下，三个双极点模型的 M2、s、OPE、连续谱、留数、f（包括 NaN）和回代结果均与原始数组完全一致。此阶段通过了十项测试。
3. 可选跑动阶数及广义 P-P 滤波器：仓库提交 `849f178`。这一改动有意改变数值结果。未声称复现原始冻结输入下的物理输出。
4. MATLAB 实现和跨语言测试记录在后续提交中。

## RG 参考值

四圈质量和耦合直接积分截断的 beta 和 gamma 函数得到。独立的 55 位逆 beta 求积及 gamma/beta 积分给出：

| 量 | 独立参考值 |
|---|---:|
| alpha_s^(5)(4.18 GeV) | 0.226141845816242854 |
| alpha_s^(4)(4.18 GeV) | 0.226374296730276341 |
| alpha_s^(4)(1.275 GeV) | 0.391619277160514745 |
| alpha_s^(4)(2 GeV) | 0.304398404402813591 |
| mc^(4)(2 GeV), GeV | 1.093276847184784040 |

按所示无量纲/GeV 单位，Python 和 MATLAB ODE 结果与这些参考值的偏差小于 4e-13。额外位数用于数值检查，并非物理输入精度。测试还恢复先前的三圈数值，验证一圈解析解、往返演化、匹配方向、移动的质量/标度锚点、无效输入，以及独立 JSON CLI 执行。

## Python 检查

将警告视为错误时，四十一项 Python 测试通过。可选 MATLAB 测试由独立开关控制。测试覆盖物理回代、精确的请求网格端点、非零边界位移、不受惩罚的极点、标准/广义奇异值转换及 QOC 有限差分。一个具有非平凡 GSVD mu 值的良态人工问题与独立组装的增广最小二乘解进行了比较；这可检出原始 P-P 中缺失 mu 平方的错误。

三个 `run.py` 入口还从另一工作目录运行，采用四圈 RG、JSON 输出和无需 pickle 的 NPZ 输出。三者使用相同的导出质量/耦合。

额外的默认网格运行采用双极点模型和四圈输入，完成后得到有限的谱与留数：

| 程序包 | M2_min | M2_max | Lambda | ds | 帽函数节点数 |
|---|---:|---:|---:|---:|---:|
| Jpsi_vector | 18.5 | 80 | 68 | 0.1 | 542 |
| eta_c_axial_vector | 18.5 | 80 | 68 | 0.1 | 532 |
| eta_c_pseudoscalar | 15.5 | 80 | 67.5 | 0.1 | 527 |

这些是执行检查，而非完整扫描或不确定度的重新估计。

## MATLAB 与 Python 一致性

在 `(5,20,20)` 和 `ds=0.1` 下，对全部三个双极点模型分别进行了三圈和四圈比较。测试比较 RG 元数据、Borel 网格、每个 OPE 分量、alpha 选择、极点留数、谱、f（包括其 NaN 掩码）及残差。alpha 比较要求 `abs(log10(alpha_matlab/alpha_python)) <= 1e-8`。

独立 MATLAB 冒烟测试检查非平凡惩罚、非零连续谱边界，以及恒等式 `A=U*C*X'`、`B=V*S*X'`，并直接检查一阶条件。最后三个双极点用例使用默认四圈输入。直接解、缩放梯度和因子重建检查均通过各自的 1e-10 容差。玩具系统也恰有两个不受惩罚的极点模式。额外的直接解检查通过 `CheckSolution=true` 显式启用。

在程序目录中重复 MATLAB 比较：

```sh
RUN_MATLAB_TESTS=1 MATLAB_PARITY_LOOPS=4 python -B -W error -m unittest discover -s tests -p test_matlab_parity.py -v
RUN_MATLAB_TESTS=1 MATLAB_PARITY_LOOPS=3 python -B -W error -m unittest discover -s tests -p test_matlab_parity.py -v
```

如果 MATLAB 不在 PATH 中，也不位于已测试的 macOS 路径，请设置 `MATLAB_EXECUTABLE`。测试会创建自己的临时结果，并以当前 Python 实现作为跨语言参考。
