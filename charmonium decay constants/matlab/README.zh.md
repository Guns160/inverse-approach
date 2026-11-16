# MATLAB 实现

**中文** | [English](README.md)

共享数值核心是本目录中的 `+charmonium` 包。它实现物理 s 变量上的帽函数基、H1 惩罚、各道的 OPE 公式、经济型 QR 紧致 GSVD、重建，以及 Raus 2024 选择器。需要 MATLAB 内置的 `gsvd`；无需工具箱或外部 GSVD 包。

调用求解器之前，将一个公共入口目录加入路径。每个入口也会加入此共享 `matlab` 目录：

```matlab
addpath("Jpsi_vector/matlab");
out = solve_jpsi_vector(5,20,20,"RGLoops",3);
```

无需运行反问题即可调用独立 RG 程序：

```matlab
addpath("matlab");
rg = charmonium.runDec(2,"Loops",3);   % alpha_s and mc at mu=2 GeV
```

三个公共求解器默认采用四圈，通过 `RGLoops` 选择阶数：

```matlab
addpath("Jpsi_vector/matlab");
addpath("eta_c_axial_vector/matlab");
addpath("eta_c_pseudoscalar/matlab");
j = solve_jpsi_vector(5,20,20,"RGLoops",3);
a = solve_eta_c_axial_vector(5,20,20,"RGLoops",3);
p = solve_eta_c_pseudoscalar(5,20,20,"RGLoops",3);
```

三个求解器均只运行 `2pole`，并在 `out.two_pole` 中返回结果。`Channel="2pole"` 是该默认值的可选显式写法。元数据、Borel 点、权重和 OPE 分量仍保存在 `out.metadata`、`out.M2`、`out.weights` 和 `out.OPE` 中。两个极点质量仍由各道分别定义：

| 求解器 | 极点质量（GeV） |
| --- | --- |
| `solve_jpsi_vector` | 3.096900, 3.686097 |
| `solve_eta_c_axial_vector` | 2.9841, 3.51067 |
| `solve_eta_c_pseudoscalar` | 2.9841, 3.6378 |

公共选项为 `Channel`、`DS`、`RGLoops`、`Mu`、`RGInputs`、`Alpha` 和 `CheckSolution`。`RGLoops` 接受 1 至 4，默认为 4；主求解器按该阶数调用独立的 `charmonium.runDec` 程序。其他默认值为 `DS=0.1`、`Mu=2` 和 `CheckSolution=false`。设置 `CheckSolution=true` 会增加独立的增广最小二乘回归和直接梯度检查。`RGInputs` 是完整的耦合结构体，字段为 `mbar_c`、`alpha_mz`、`mz` 和 `mb`；粲夸克质量及其锚点标度均使用 `mbar_c`。

## 赝标量轴反常

仅 PP 入口接受 `Anomaly1S` 和 `Anomaly2S`（GeV^3）。默认值为
`Anomaly1S=-0.079` 和 `Anomaly2S=[]`。空输入使该态保持
仅含密度贡献；零则显式提供零矩阵元。

```matlab
p = solve_eta_c_pseudoscalar(5,20,20); % ground-state anomaly included
raw = solve_eta_c_pseudoscalar(5,20,20,"Anomaly1S",[]);
```

`two_pole.f`、`f_all` 和 `f_eta_c` 包含所提供的 `a_i/m_i^2` 项。
`f_raw` 和 `f_raw_all` 保留 `2*mc*sqrt(R_i)/m_i^2`，
`anomaly_shift_GeV` 记录位移。`metadata.anomaly` 记录输入，
以及哪些态仍仅含密度贡献。绘制 `f` 时不要再次添加修正。
OPE 和谱反演不依赖反常输入。

加入本目录后，运行 `test_pp_anomaly` 以验证 Ward 恒等式、
可选输入，以及谱反演保持不变。Python/MATLAB 一致性
检查也会运行此测试。

使用以下命令运行自包含数值冒烟测试：

```matlab
addpath("matlab");
result = test_smoke;
```

它将具有两个不受惩罚极点和非零端点的广义滤波器与直接增广最小二乘解进行比较。随后检查全部三个物理入口的默认四圈输入、独立驻定性、GSVD 因子恒等式，以及唯一的 `2pole` 结果契约。`tests/test_matlab_parity.py` 中的 Python/MATLAB 跨语言测试由于会启动 MATLAB，需由 `RUN_MATLAB_TESTS=1` 启用。设置 `MATLAB_PARITY_LOOPS=3` 或 `MATLAB_PARITY_LOOPS=4` 以选择验证阶数。
