# 反问题方法

**中文** | [English](README.md)

用于从 Borel 求和规则提取粲偶素极点留数和衰变常数的科研代码。

程序位于 [`charmonium decay constants/`](charmonium%20decay%20constants/)：

| 程序 | 流 | 显式极点模型 |
|---|---|---|
| `Jpsi_vector` | 矢量流 | J/psi + psi(2S) |
| `eta_c_axial_vector` | 轴矢量流 | eta_c + chi_c1 |
| `eta_c_pseudoscalar` | 赝标量密度 | eta_c + eta_c(2S) |

轴矢量道的第二个极点是 chi_c1，而非 eta_c(2S)。程序保留了这三个道各自的 OPE 和留数约定。
三个程序均使用模型标识符 `2pole`。

数值输出取决于跑动耦合处方、谱网格、Borel 窗口、极点模型和正则化规则。求解成功本身并不确立物理稳定性窗口。

## 版本历史

- `ed20389`：最初上传的 Python 源码。
- `06fb9a2`：按粒子/流组织的目录布局，以及修正后的入口导入。
- 后续提交分别处理公共核心重构、跑动标度更新和 MATLAB 实现。

源码包不包含历史结果和参数扫描。

## 求解器改进示例

[`improvement/`](improvement/) 包含一个最小化、自包含的 Python + MATLAB 求解器改进示例：以一次经济型 QR + 紧致 GSVD 替代原先逐 alpha 求解正规方程的过程，并通过向量化广义滤波器生成整个 alpha 网格。MATLAB 示例引用了其替代的原始逐 alpha 求解循环。

## 不确定度分析笔记

[`uncertainty-analysis/`](uncertainty-analysis/) 记录了此类反问题的不确定度分析方法。误差按来源划分为方法外部与方法内部，每类误差各有一篇笔记：

| 笔记 | 误差类型 | 分类 |
|---|---|---|
| [`input-sampling.md`](uncertainty-analysis/input-sampling.md) | 物理输入的高斯抽样 | 输入诱导 |
| [`ope-truncation.md`](uncertainty-analysis/ope-truncation.md) | OPE 截断 | 输入诱导 |
| [`mu-scale.md`](uncertainty-analysis/mu-scale.md) | 重整化标度（固定，不作为误差） | 输入诱导 |
| [`platform-parameters.md`](uncertainty-analysis/platform-parameters.md) | 平台参数 | 方法内禀 |
| [`alpha-regularization.md`](uncertainty-analysis/alpha-regularization.md) | 正则化近似 | 方法内禀 |
| [`reference-comparison.md`](uncertainty-analysis/reference-comparison.md) | 与参考研究的输入列表差异 | — |

每篇笔记给出操作流程、检查项和报告格式。由于算子、网格和输出量各不相同，实现保留在各个项目中。[`README.md`](uncertainty-analysis/README.md) 为索引。

这些笔记还记录了为各误差源指定的分布：具有测量中心值及不确定度的输入采用高斯分布，平台参数在平台内采用均匀分布，重整化标度、OPE 截断和正则化近似则完全不抽样。无论逐源抽样还是联合抽样，这一分配都保持不变。
