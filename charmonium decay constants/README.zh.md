# 粲偶素衰变常数

**中文** | [English](README.md)

三个道的程序包使用公共数值基础设施和共享的 MSbar QCD 跑动处方。输入和输出以 GeV 及 GeV 的幂为单位。

| 程序包 | 极点模型名称 | 第二个极点 | 留数转换 |
|---|---|---|---|
| `Jpsi_vector` | `2pole` | psi(2S) | `f = sqrt(R)` |
| `eta_c_axial_vector` | `2pole` | chi_c1 | `f = sqrt(R)` |
| `eta_c_pseudoscalar` | `2pole` | eta_c(2S) | 提供 a_i 时为 `f = (2*mc(mu)*sqrt(R)+a_i)/m_i^2` |

非正的极点留数产生 `NaN` 衰变常数。极点留数和重建的连续谱均未施加正性约束。

## Python

需要 Python 3.9+、NumPy 和 SciPy。通过 `requirements.txt` 安装 Python 依赖。GSVD 后端还需要一个导出 LP64 C 接口 `LAPACKE_dggsvd3` 的库。它会检测受支持的系统、NumPy 和 SciPy OpenBLAS 库；也可通过 `JPSI_LAPACKE_LIBRARY` 显式指定库。`ETA_C_LAPACKE_LIBRARY` 保留为兼容性别名。没有正规方程或 Gram 矩阵回退路径。

在本目录中运行：

```sh
python -m Jpsi_vector 5 20 20 --ds 0.1 --mu 2 --rg-loops 4
python -m eta_c_axial_vector 5 20 20 --ds 0.1 --mu 2 --rg-loops 4
python -m eta_c_pseudoscalar 5 20 20 --ds 0.1 --mu 2 --rg-loops 4
```

位置参数为 `M2_min`、`M2_max` 和 `Lambda`，单位均为 GeV²。每个程序计算其双极点模型，标识为 `2pole`。这些小窗口用于演示执行过程，并非推荐的物理提取窗口。

```python
from Jpsi_vector import solve
from inverse_common.rg import RGInputs

result = solve(5, 20, 20, ds=0.1,
               mu=2.0, rg_loops=4, rg_inputs=RGInputs())
channel = result["channels"]["2pole"]
print(result["metadata"]["rg"])
print(channel["f"], channel["alpha_star"])
```

API 返回数据而不写入文件。CLI 运行会创建包含 `result.npz` 和 `summary.json` 的唯一目录；`--output` 用于选择输出根目录。NPZ 值无需 pickle。JSON 将非有限值替换为 `null`。

## 赝标量轴反常转换

PP 求解器现在将轴 Ward 转换纳入主要结果：
`f_i = f_raw_i + a_i/m_i^2`，其中 `f_raw_i = 2*mc(mu)*sqrt(R_i)/m_i^2`。
这里 `m_i` 是强子质量，而非 Borel 参数。默认基态
矩阵元为 `a_1=-0.079 GeV^3`，取自 Qin 等人，Phys. Rev. D 97, 096002
(2018)，式 (22)。其位移约为 -8.872 MeV。不假定 2S 矩阵元：
默认情况下，第二项仍只含密度贡献，这一状态记录在
`metadata.anomaly.applied` 和 `metadata.anomaly.f_definition` 中。

```python
from eta_c_pseudoscalar import solve
out = solve(5, 20, 20)  # a_1=-0.079; a_2 omitted
raw = solve(5, 20, 20, anomaly_1s=None)  # omit both anomaly inputs
custom = solve(5, 20, 20, anomaly_1s=-0.079, anomaly_2s=0.02)
```

上述自定义 2S 值仅用于演示 API，并非本项目推荐的物理输入。
矩阵元以 GeV^3 为单位提供，并采用与 `h_i=+sqrt(R_i)` 相同的态相位约定。
`None` 表示缺失，而零表示显式提供零矩阵元。

```sh
python -m eta_c_pseudoscalar 5 20 20 --anomaly-1s none
python -m eta_c_pseudoscalar 5 20 20 --anomaly-1s -0.079 --anomaly-2s 0.02
```

PP 输出模式 v3 使用 `f`、`f_all` 和 `f_eta_c` 表示已计入所提供
反常项的值。`f_raw` 和 `f_raw_all` 保留密度贡献；
`anomaly_shift_GeV` 给出每个态的加性位移。汇总 JSON 同时包含
`f_GeV` 和 `f_raw_GeV`。非正留数仍为未定义。反常转换
不改变 OPE、留数、谱和 alpha 选择。
读取 `f` 的绘图程序不得再次添加同一修正；
自行执行修正的旧后处理程序应改为读取 `f_raw`。

## 跑动标度

跑动由独立程序 `inverse_common/rg.py` 实现，无需运行反问题即可使用：

```sh
python -m inverse_common.rg --mu 2 --loops 4
python -m inverse_common.rg --mu 2 --loops 3
```

输出为 JSON；省略 `--loops` 时选择四圈。三个求解器通过 `solve(..., rg_loops=4)` 或 `--rg-loops 4` 调用同一个模块。

默认输入为 `alpha_s^(5)(MZ)=0.1184`、`MZ=91.1876`、`mbbar(mbbar)=4.18`、`mcbar(mcbar)=1.275` 和 `mu=2`。三个求解器均使用四圈 beta 函数、质量反常维数函数及三圈底夸克阈值匹配。支持的物理目标范围为 `mcbar <= mu <= mbbar`，粲夸克保持活跃（`nf=4`）。

同一组导出的质量和耦合用于 OPE、微扰边界和赝标量留数转换。改变 `RGInputs.mbar_c` 会同时改变锚点质量及其标度。调用不会修改全局物理参数。

`rg_loops=1,2,3,4` 选择一致的跑动阶数，匹配阶数为 `rg_loops-1`。历史 eta_c 领头对数质量处方与 `rg_loops=1` 并不等价；原始版本保留在 Git 历史中。

公式、匹配方向和数值参考见 [RG 约定](docs/RG_CONVENTIONS.zh.md)。四圈 RG 不会提高 OPE 截断阶数：VV 保留现有 NLO 因子，A-A 保留 LO 加凝聚项，P-P 保留现有 LO/NLO 加凝聚项。

## 数值约定

- Borel 步长为 0.01 GeV²；精确包含所请求的端点。
- 默认谱帽函数步长为 0.1 GeV²。P-P 历史生产网格使用 0.01 GeV²；这是不同的离散化，对应规模显著更大的 GSVD。
- 连续谱起点为 `S_H=13.9105153024`（VV）或 `14.99`（A-A 和 P-P）。这些是各道的模型输入。
- 目标函数为 `||sqrt(W)*(K*z-g)||² + alpha*rho_hat'*H*rho_hat`，其中 `H` 为质量加刚度矩阵。Alpha 不取平方。极点留数不受惩罚。
- 连续谱的两个端点分别固定为零和 Lambda 处的微扰谱密度。保留端点交叉项。
- VV 和 A-A 使用各自的自适应 alpha 网格；P-P 使用历史 2001 点 `10^[-20,0]` 网格。各道均使用 Raus Q/MQ/HR 规则。
- P-P 的全 alpha 重建使用广义滤波器 `sigma/(sigma²+alpha*mu_gsvd²)`。其选择器使用 `sigma/mu_gsvd`。这修正了原始独立 P-P 滤波器约定。

共享基础设施位于 `inverse_common/`。三个程序包目录保留各自的物理核和极点定义。复制或安装代码时，请将共享目录与这些程序包放在一起。

## MATLAB

三个 MATLAB 版本使用共享的 `matlab/+charmonium` 包和 MATLAB 内置的 `gsvd`。无需单独的 GSVD 工具箱。在本目录中运行：

```matlab
addpath("Jpsi_vector/matlab", "eta_c_axial_vector/matlab", ...
        "eta_c_pseudoscalar/matlab", "matlab");
rg = charmonium.runDec(2, "Loops", 4);
j = solve_jpsi_vector(5, 20, 20, "RGLoops", 4);
a = solve_eta_c_axial_vector(5, 20, 20, "RGLoops", 4);
p = solve_eta_c_pseudoscalar(5, 20, 20, "RGLoops", 4);
```

省略 `Loops` 或 `RGLoops` 时选择四圈。将其设为 1、2 或 3 可改变阶数。由于字段名不能以数字开头，三个 MATLAB 程序均在 `two_pole` 中返回双极点结果。例如，`p.two_pole.f` 存放两个赝标量衰变常数。公共标度信息位于 `p.metadata.rg`。

选项和验证说明见 [MATLAB 使用说明](matlab/README.zh.md)。`CheckSolution=true` 启用额外的增广最小二乘比较；正常计算时关闭。

## 验证

```sh
python -B -W error -m unittest discover -s tests -v
```

测试覆盖 RG 系数和独立求积参考值、参数传播、非零边界项、驻定性、留数约定、QOC 导数和序列化。MATLAB 一致性测试为需显式启用的可选测试，相关说明随 MATLAB 实现提供。

已测试的环境、测试窗口及本次验证的范围见[验证记录](docs/VALIDATION.zh.md)。
