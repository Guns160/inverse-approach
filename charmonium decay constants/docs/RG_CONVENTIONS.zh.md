# 共享四圈 QCD RG 约定

**中文** | [English](RG_CONVENTIONS.md)

## 来源与范围

来源：K. G. Chetyrkin、J. H. Kuhn 和 M. Steinhauser，**RunDec：用于强耦合和
夸克质量跑动与退耦的 Mathematica 程序包**，hep-ph/0004189v1（19 April 2000）。
已通过提取文本阅读所提供的 PDF；还在渲染后的 PDF 第 6 和 11 页
（印刷页码第 5 和 10 页）核对了式 (7) 和式 (20)-(21)。

代码在 MSbar 重夸克自身标度实现式 (1)、(2)、(6)、(7)、(19) 和 (20)/(21)。
所有质量和标度均以 GeV 为单位。
这升级的是**输入的跑动**；不改变 OPE 的微扰阶数、
凝聚量演化、谱模型或反问题求解器。

## Python 求解器集成 API

```python
from inverse_common.rg import RGInputs, derive, evolve_fixed_nf

r = derive()  # mu=2 GeV, loops=4
mc, alpha_s = r["mc"], r["alpha_s"]
r3 = derive(loops=3)  # compare with frozen three-loop physics
r_varied = derive(inputs=RGInputs(mbar_c=1.305))
# Low-level mathematical fixed-nf evolution, no automatic thresholds:
alpha_end, mass_end = evolve_fixed_nf(alpha_s, mc, 2.0, 3.0, nf=4, loops=4)
```

- `RGInputs(mbar_c=1.275, alpha_mz=0.1184, mz=91.1876, mb=4.18)` 是
  冻结且可哈希的数据类。输入会经过验证并转换为浮点数。
- `derive(mu=2.0, loops=4, inputs=None) -> dict`：字段包括 `mu`、
  `mc`、`alpha_s`、`nf=4`、`loops`、`matching_loops`、`inputs`、
  `input_source`、`coefficient_source`、`mass_method="ode"`、`mass_anchor`、`threshold` 及
  积分容差。每次调用返回新的可变元数据。
- `beta_coefficients(nf)` 和 `gamma_coefficients(nf)` 返回不可变
  元组，采用下述归一化，索引从 0 到 3。
- `decouple_alpha_down(alpha_s, nf=5, matching_loops=3)` **仅在**
  `mu=mhbar(mhbar)` 处返回低味数耦合。
- `evolve_fixed_nf(alpha_s, mass, start, end, nf=4, loops=4)` 返回
  `(alpha_s(end), mass(end))`。它固定显式给出的 nf，按所请求的 beta
  和 gamma 阶数精确积分，支持两个演化方向。

`derive` 支持 `mbar_c <= mu <= mb`。在 `mu=mb` 处，它返回阈值
nf=4 **低味数侧**的值。区间外的目标会触发错误。
未实现跨越底夸克阈值的轻夸克质量匹配。粲夸克始终保持
活跃：在粲夸克锚点处不进行 nf=4 到 nf=3 的退耦。底层
固定 nf API 是数学运算，不推断阈值。

输入链为：

1. 将 `alpha_s^(5)(MZ)=0.1184` 从 `MZ=91.1876` 演化至 `mb=4.18`，nf=5。
2. 在 `mu=mbbar(mbbar)=4.18` 处将耦合从 nf=5 匹配至 nf=4。
3. 将 nf=4 耦合从 mb 演化至**移动的**锚点 `mu=mbar_c`。
4. 施加 `mc^(4)(mbar_c)=mbar_c`，随后将耦合和质量演化至目标 mu。

步骤 1 和 3 在固定 nf 例程中使用值为 1 的占位质量；该质量随后
丢弃。物理粲夸克质量不跨越底夸克阈值。改变
`mbar_c` 会同时改变其边界质量和边界标度。改变 `mb`
会同时改变 MSbar 底夸克自身质量和阈值；这并非在固定底夸克质量下
任意改变匹配标度。

## 独立跑动程序（Python 和 MATLAB）

一个独立 RG 程序支持 1、2、3 和 4 阶。调用者选择
阶数，无需在每个反问题求解器中重复 RG 系数。
默认采用四圈。在项目目录中运行：

```sh
python -m inverse_common.rg --mu 2 --loops 4
python -m inverse_common.rg --mu 2 --loops 3
python inverse_common/rg.py --mu 2 --loops 4
python -m inverse_common.rg --mu 2 --loops 4 --mcbar 1.305 --alpha-mz 0.119 --mb 4.2 --mz 91.2
```

两个入口均向 stdout 打印一个标准 JSON 对象，包含
物理结果、输入来源、圈数/匹配阶数，以及
`mass_method="ode"`。它们不写入结果文件。锚点默认值与
`RGInputs` 一致；`--mu` 默认为 2 GeV，`--loops` 默认为 4。
使用 `--help` 列出选项。无效参数或定义域会向 stderr 报告清晰的
错误，并以状态码 2 退出。积分失败时以状态码 1 退出。
发生错误时不输出不完整的 JSON 结果或 traceback。

MATLAB 已提供独立跑动入口：

```matlab
r4 = charmonium.runDec(2, "Loops", 4);
r3 = charmonium.runDec(2, "Loops", 3);
r_default = charmonium.runDec(2); % four loops by default
```

反问题程序的调用者选择跑动阶数；参数变化时
还可提供可选锚点。

系统 Python 验证通过了全部 19 个 RG/CLI unittest 方法。CLI 检查
对两个入口、1-4 圈阶数、默认值、自定义锚点、JSON 输出、
空临时输出目录、帮助、无效输入（状态码 2）和积分失败（状态码 1）
运行真实子进程。

## MATLAB 所用的精确系数

令 `a=alpha_s/pi`、`t=ln(mu^2)`、`n=nf`。

```text
da/dt       = -beta0*a^2 - beta1*a^3 - beta2*a^4 - beta3*a^5
d ln(mc)/dt = -gamma0*a - gamma1*a^2 - gamma2*a^3 - gamma3*a^4

beta0 = (11 - 2*n/3)/4
beta1 = (102 - 38*n/3)/16
beta2 = (2857/2 - 5033*n/18 + 325*n^2/54)/64
beta3 = (149753/6 + 3564*zeta3
         + (-1078361/162 - 6508*zeta3/27)*n
         + (50065/162 + 6472*zeta3/81)*n^2
         + 1093*n^3/729)/256

gamma0 = 1
gamma1 = (202/3 - 20*n/9)/16
gamma2 = (1249 + (-2216/27 - 160*zeta3/3)*n - 140*n^2/81)/64
gamma3 = (4603055/162 + 135680*zeta3/27 - 8800*zeta5
          + (-91723/27 - 34192*zeta3/9 + 880*zeta4 + 18400*zeta5/9)*n
          + (5242/243 + 800*zeta3/9 - 160*zeta4/3)*n^2
          + (-332/243 + 64*zeta3/27)*n^3)/256

zeta3 = 1.2020569031595942854
zeta4 = pi^4/90
zeta5 = 1.0369277551433699263
```

这些系数不额外乘以 2：自变量已是
`ln(mu^2)`。若在 `ln(mu)` 上积分，则将两个
右端均乘以 2。积分 `ln(mc)` 以保持质量为正。
在 N 圈时，将两个系数向量截断至索引 `0..N-1`。

| nf | beta0 | beta1 | beta2 | beta3 |
|---|---:|---:|---:|---:|
| 4 | 2.08333333333333333 | 3.20833333333333333 | 6.34924768518518519 | 31.38744695230514049 |
| 5 | 1.91666666666666667 | 2.41666666666666667 | 2.82667824074074074 | 18.85217315933944027 |

| nf | gamma0 | gamma1 | gamma2 | gamma3 |
|---|---:|---:|---:|---:|
| 4 | 1 | 3.65277777777777778 | 9.94704026107295732 | 27.30293363421854576 |
| 5 | 1 | 3.51388888888888889 | 7.41986321214366579 | 11.03432124206805957 |

## 精确匹配约定

RunDec 式 (19) 为 `alpha_low = zeta_g^2 * alpha_high`。在
MSbar 自身质量处，对数 `L=ln(mu^2/mhbar(mu)^2)=0`。式 (20) 和
(21) 在此处一致：

```text
a_high = alpha_high/pi
nl = nf_high - 1
d1 = 0
d2 = 11/72
d3 = 564731/124416 - (82043/27648)*zeta3 - (2633/31104)*nl
alpha_low = alpha_high * (1 + d2*a_high^2 + d3*a_high^3)
```

在底夸克阈值处 `nl=4`、`d3=0.633450719548517290326579933996`。匹配多项式
使用**高味数侧**的耦合。这并非极点质量
匹配公式式 (22)。N 圈跑动采用 N-1 圈匹配：

| loops | matching_loops | L=0 时 zeta_g 平方中的项 |
|---|---|---|
| 1 | 0 | 1 |
| 2 | 1 | 1 (d1=0) |
| 3 | 2 | 1 + d2*a_high^2 |
| 4 | 3 | 1 + d2*a_high^2 + d3*a_high^3 |

## 数值结果与验证

默认输入，`mu=2 GeV`；以下为 RG 输出，并非拟合得到的衰变常数：

| loops | matching_loops | alpha_s^(4)(2) | mc^(4)(2) GeV |
|---|---|---:|---:|
| 1 | 0 | 0.26976783425397743 | 1.171916344348932 |
| 2 | 1 | 0.2992873428686697 | 1.115375821681162 |
| 3 | 2 | 0.30294509881312925 | 1.0990536546335794 |
| 4 | 3 | 0.3043984044027081 | 1.0932768471848942 |

生产积分采用 SciPy `solve_ivp(method="DOP853")`，
`rtol=2e-12`、`atol=2e-14`。它直接积分有限阶 beta/gamma
ODE，而非使用 RunDec 式 (10) 中截断的质量比展开。
基于该展开的实现可能因形式上的高阶项而产生差异，
并不属于相同的数值约定。
结果显式记录 `mass_method="ode"`。2026-10-14 提供的独立物理
审计通过 55 位求积确认了 ODE 参考值。对相同四圈输入，截断的 `c(a)` 比值给出
`mc(2)=1.092612428731451 GeV`，低 `0.664418453333 MeV`。该值
不是 ODE 回归目标。

一项独立的 50 位小数 mpmath 计算采用精确有理数
系数和 `mp.zeta`，改在耦合空间中求解。定义
`B(a)=sum(beta_i*a^(i+2))` 和 `G(a)=sum(gamma_i*a^(i+1))`，求解

```text
integral[a_start -> a_end] da/B(a) + 2*ln(mu_end/mu_start) = 0
m_end = m_start * exp(integral[a_start -> a_end] G(a)/B(a) da)
```

参考结果为：

```text
alpha_s^(5)(4.18)   = 0.226141845816242853648784456086645162252737161
alpha_s^(4)(4.18)   = 0.226374296730276340625973984089232024321228921
alpha_s^(4)(1.275)  = 0.391619277160514744523610869382628899904314963
alpha_s^(4)(2)      = 0.304398404402813591292861525129163528655322391
mc^(4)(2)          = 1.09327684718478403978880450988362917080470269 GeV
```

所列各量的 DOP853/参考值差异均低于 2e-13。
这些参考值来自独立求积，而非已安装的 RunDec
二进制程序。临时安装 `rundec==0.7` 验证包失败，
因为其发行包缺少 `CRunDec.3.1.h`；因此没有声称完成 RunDec 库交叉核对。
独立的精确系数和数值方法提供了四圈参考值。

测试覆盖冻结的三圈数值（绝对容差 2e-11）、
独立四圈数值、nf=4/5 系数、解析一圈
跑动、1-4 圈往返与复合演化、阈值修正/阶数、
mb 的 nf=4 侧、移动粲夸克锚点 1.245/1.275/1.305、改变圈数时的实际数值
变化、不可变输入/缓存隔离、无效输入
阶数、不支持的目标范围，以及非微扰积分保护。

```text
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 -m unittest discover -s tests -p test_rg.py -v
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 -m unittest discover -s tests -v
```

验证环境：系统 `/usr/bin/python3`（Python 3.9.6）、
NumPy 2.0.2、SciPy 1.13.1 和标准库 `unittest`；独立
参考计算使用 mpmath 1.4.1。NumPy/SciPy 是运行时依赖。
测试套件不要求安装 pytest 或 mpmath；其数值
参考值已冻结，实时求积则使用运行时 SciPy。
此次变更未向项目安装任何依赖。

## 缓存与失败行为

内部有界 LRU 缓存存储不可变元组：128 组边界
配置 `(RGInputs, loops)` 和 512 组固定 nf 演化。缓存键
包括每个数值输入和圈数。不存在可变的缓存
结果字典。拒绝非有限/非正的标度、质量和耦合；
拒绝 bool 以及非整数的味数/圈数。拒绝耦合 `alpha_s>=pi`，
作为宽松的微扰保护；这并不保证靠近该边界时的精度。
积分到达 `a=1` 时停止并抛出错误。
求解器失败和无效的最终状态会引发明确错误。

## MATLAB RG 实现与实测验证

单独授权的 MATLAB RG 修复使用以下包文件：
`matlab/+charmonium/runDec.m`、`qcdCoefficients.m`、`integrateRG.m` 和
`matchAlphaDownMS.m`。错误的 `matchAlphaDownOS.m` 已删除。

```matlab
r = charmonium.runDec(2, 'Loops', 4, 'AlphaMZ', .1184, ...
    'MZ', 91.1876, 'Mb', 4.18, 'McBar', 1.275);
[alpha_end, mass_end] = charmonium.integrateRG(.3, 1.1, 2, 3, 4, 4);
c = charmonium.qcdCoefficients(4, 4);
alpha_low = charmonium.matchAlphaDownMS(alpha_high, 5, 3);
```

`runDec` 返回与 Python 相同的物理字段名和约定，
包括 `mass_method='ode'`、`inputs`、`mass_anchor` 和 `threshold`。
目标范围为 `[McBar,Mb]`。没有独立的 `MassScale`
选项。`integrateRG` 接受/返回 alpha_s，将内部状态转换
为 `[alpha_s/pi; ln(mass)]`，并在 `2*ln(mu)` 上使用自适应
`ode113` 积分，`RelTol=2e-12`、`AbsTol=2e-14`。

2026-10-14 使用 MATLAB R2026a 实测：

| loops | alpha_s^(4)(2) | mc^(4)(2) GeV |
|---|---:|---:|
| 3 | 0.30294509881300485 | 1.099053654633755 |
| 4 | 0.30439840440251487 | 1.0932768471851821 |

四圈底夸克阈值：高味数侧 `0.22614184581613855`，低味数侧
`0.22637429673017168`。每项四圈参考值偏差均低于
4e-13。临时 MATLAB 验证通过了绝对容差为 2e-11 的三圈/四圈数值
参考值、nf=4/5 gamma3、1-4 圈往返/复合演化、
移动粲夸克锚点、mb 的低味数侧、零对数匹配
阶数、无效标度/阶数，以及对 `MassScale` 的拒绝测试。

系统 Python `unittest` 发现并通过了全部 13 个 RG 测试方法
（先前的参数化用例保留为子测试）。较早的一次
完整套件发现过程出现了两个不相关的失败：`test_contract` 中的 ndarray JSON
序列化，以及导入 `test_matlab_parity` 时的语法错误。
这些文件属于另一位实现代理；此阶段未声称该次完整套件运行通过。
未操作 Git 状态和生产数据集。
