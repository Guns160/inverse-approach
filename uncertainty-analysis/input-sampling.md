# 输入抽样误差

输入误差来自方法外部。它不依赖反问题求解流程本身，换一个反演程序，这部分误差依然存在。

本文件只讲**高斯撒点**这一类。其余三类误差各有单独文件：

| 误差种类 | 文件 |
|---|---|
| 平台参数 | [platform-parameters.md](platform-parameters.md) |
| 输入高斯撒点 | 本文件 |
| OPE 截断 | [ope-truncation.md](ope-truncation.md) |
| 重整化标度 $\mu$ | [mu-scale.md](mu-scale.md) |

> **本文只记录方法。** 文中的区间和代码里的数值都是示例或当前取值，不代表最终选择。输入的中心值和 $1\sigma$ 是已确认的，可以直接使用。

## 一、撒点清单

本课题进入高斯撒点的输入共 5 项。

| 输入 | 代码位置 | 中心值 | $1\sigma$ | 出处 |
|---|---|---|---|---|
| $G_2$ | `physics.G2` | $0.038\ \mathrm{GeV}^4$ | 0.011 | 1801.06832 |
| $G_3$ | `physics.G3` | $0.013\ \mathrm{GeV}^6$ | 0.007 | 1801.06832 |
| $\bar m_c(\bar m_c)$ | `RGInputs.mbar_c` | 1.275 GeV | 0.020 GeV | 项目输入表 |
| $M_{1S}$ | `MASSES` | 见下表 | 见下表 | PDG |
| $M_{2S}$ | `MASSES` | 见下表 | 见下表 | PDG |

极点质量按通道不同：

| 通道 | $M_{1S}$ | $\sigma$ | $M_{2S}$ | $\sigma$ |
|---|---|---|---|---|
| `Jpsi_vector` | $J/\psi$ 3.096900 | $6\times10^{-6}$ GeV | $\psi(2S)$ 3.686097 | $6\times10^{-5}$ GeV |
| `eta_c_axial_vector` | $\eta_c$ 2.9841 | $3\times10^{-4}$ GeV | $\chi_{c1}$ 3.51067 | $5\times10^{-5}$ GeV |
| `eta_c_pseudoscalar` | $\eta_c$ 2.9841 | $3\times10^{-4}$ GeV | $\eta_c(2S)$ 3.6378 | $6\times10^{-4}$ GeV |

$\bar m_c$ 的变化通过 RG 传播，同时改变 $m_c(\mu)$ 和微扰谱，不只是 OPE 里的显式质量因子。前一版预算曾漏掉后一部分，导致输入误差偏小。重算时必须确认所有用到 $m_c$ 的位置都重新取值。

不进入本文件的输入：

| 输入 | 处理 | 文件 |
|---|---|---|
| $\alpha_s(2\ \mathrm{GeV})=0.3029$ | 固定，不赋误差 | 本文件，见下文 |
| $M_Z$、$\bar m_b$ | 仅作 RG 内部端点 | 本文件 |
| $\mu$ | 固定或诊断扫描 | [mu-scale.md](mu-scale.md) |
| $a_1=-0.079$ GeV³ | 暂不赋误差 | 本文件 |
| OPE 截断 | 差分 | [ope-truncation.md](ope-truncation.md) |

$\alpha_s$ 已由 RG 跑到计算能标，项目约定不计独立误差。$a_1$ 需要先确定来源才能撒点。

## 二、分布

### 2.1 单参数

给定中心值 $\theta_0$ 和误差 $\sigma$，抽

$$
\theta = \theta_0 + \sigma z, \qquad z \sim \mathcal{N}(0,1)
$$

误差按 $1\sigma$ 解释。这是“测量值服从正态分布”的标准含义。

### 2.2 非对称误差

上下误差不同时，用分裂正态：先抽 $z$，再按符号选误差。

$$
\theta = \theta_0 + z \times
\begin{cases}
\sigma_{-}, & z < 0\\
\sigma_{+}, & z \ge 0
\end{cases}
$$

两侧各占一半概率质量。中心值因此是中位数，不是均值。报告时两套中心值都要写清楚。

本课题目前 5 项输入的上下误差都对称，用不到这一条。保留它是为了以后加入非对称输入。

### 2.3 相关性

多个输入相关时，用相关矩阵 $C$。先做特征分解 $C=V\sqrt{D}$，再抽

$$
\boldsymbol\theta = \boldsymbol\theta_0 + V\sqrt{D}\,\boldsymbol z,
\qquad \boldsymbol z\sim\mathcal{N}(\mathbf 0,I)
$$

$C$ 必须对称、对角元为 1、半正定。相关性从哪来必须有依据，不能凭感觉设。

**本课题暂按独立处理。** 报告中必须写明“假设独立”。

## 三、两套撒点模式

同一批输入可以按两种方式组织。两者都是高斯撒点，区别在于一次动几个输入。

**分布类型与模式无关。** 单参数和联合只是组织方式，5 项输入都是高斯，不会因为换了模式就变成别的分布。平台参数则始终是均匀，见 [platform-parameters.md](platform-parameters.md)。

### 3.1 模式一：单参数撒点（逐项）

一次只动一个输入，其余固定在中心值。每个输入各得一个输出分布，最后各项按平方和合成。

```text
PROCEDURE per_source(inputs, N):
    FOR each input i:
        values ← []
        FOR k IN 1..N:
            θ ← inputs with only θ_i resampled
            y ← SOLVE(θ)
            IF finite(y): values.append(y)
        Δ_i ← std(values)
    Δ_total ← sqrt(Σ_i Δ_i²)
    RETURN {Δ_i}, Δ_total
```

特点：

- 能出分项表，看得出哪一项主导。
- 隐含假设各输入独立、响应近似线性。
- 计算量随输入个数线性增加。

### 3.2 模式二：联合撒点

所有输入同时抽，每个样本得到一组完整的输入，跑一次反演。

```text
PROCEDURE joint(inputs, N):
    values ← []
    FOR k IN 1..N:
        θ ← resample every input at once
        y ← SOLVE(θ)
        IF finite(y): values.append(y)
    Δ_total ← std(values)
    RETURN Δ_total
```

特点：

- 只给一个总的不确定度，给不出分项。
- 自动包含输入之间的相关性和响应的非线性。
- 计算量与输入个数无关，只由样本数决定。

### 3.3 两套模式的关系

| 情形 | 单参数结果 | 联合结果 |
|---|---|---|
| 独立、线性 | $\sqrt{\sum\Delta_i^2}$ | 与左边一致 |
| 有相关性 | 平方和偏大或偏小 | 直接给出联合分布 |
| 非线性明显 | 只是一阶近似 | 更接近真实分布 |

所以两者不是互相替代，而是用途不同。

### 3.4 怎么选

- 要出误差预算表，逐项列来源，用模式一。
- 要报一个总的不确定度，或输入之间有相关性，用模式二。
- 两种都做：模式一出分项表，模式二核对总和。若两者差得远，说明线性或独立假设不成立，报告中必须说明。

本课题的 $\bar m_c$ 响应在 1S 和 2S 上方向不同、量级差一倍以上，属于非线性明显的情形，建议两种都做。

## 四、算法

```text
PROCEDURE input_error(theta0, sigma_minus, sigma_plus, n_samples, seed):
    rng.seed(seed)
    values ← empty list
    FOR i IN 1..n_samples:
        z ← N(0,1)
        theta ← theta0 + z * (sigma_minus if z < 0 else sigma_plus)
        IF theta violates physical constraint:
            record and resample
        y ← SOLVE(theta)
        IF y is finite:
            values.append(y)
    RETURN mean(values), std(values), quantiles(values, [0.025,0.1587,0.5,0.8413,0.975])
```

要点：

- `SOLVE` 每次调用内部自动重选 $\alpha$。撒点代码不干预 $\alpha$。
- 约束检查只在生成输入时做。进入反演后失败的样本不补抽。
- 标准差用样本标准差（分母 $n-1$），不是标准误。
- 有效样本数和被拒样本数都要记录。

## 五、统计量

### 5.1 输出统计

对每个输出量分别统计：

```text
mean     ← arithmetic mean of finite values
sd       ← sample standard deviation (n-1)
q_lo     ← 15.87th percentile
q_med    ← median
q_hi     ← 84.13rd percentile
err_minus ← q_med - q_lo
err_plus  ← q_hi - q_med
```

报样本标准差，不报标准误。这里要的是“输入变化时结果能变多少”，不是“均值的精度”。样本数增加不会让这个数变小。

### 5.2 检查

每个参数跑完后检查四项：

1. **有效样本数。** 明显少于 `n_samples` 时，报出的统计量只代表子集。
2. **分布形状。** 直方图是否单峰、是否严重偏斜。偏斜明显时分位数比标准差更可靠。
3. **非物理输出。** 极点残量为负导致的 `NaN` 数量。这个数量本身随参数变化，要记录。
4. **$\alpha$ 是否跑出平台。** 若某些样本的 $\alpha$ 选到端点，这些样本的可靠性下降。

## 六、报告格式

每个通道、每个输出量分别一张表。

| 输入 | 中心值 | $1\sigma$ | $\Delta y_1$ | $\Delta y_2$ | 有效样本 |
|---|---|---|---|---|---|
| $G_2$ | $0.038\ \mathrm{GeV}^4$ | 0.011 | — | — | — |
| $G_3$ | $0.013\ \mathrm{GeV}^6$ | 0.007 | — | — | — |
| $\bar m_c(\bar m_c)$ | 1.275 GeV | 0.020 GeV | — | — | — |
| $M_{1S}$ | 见第一节 | 见第一节 | — | — | — |
| $M_{2S}$ | 见第一节 | 见第一节 | — | — | — |

表下方必须写明：撒点模式、样本数、随机种子、分布类型、相关性假设、约束拒绝率、无效样本数、输出单位，以及分项是否按平方和合成。

## 七、注意事项

**1. 有效样本数必须报。**

极点残量为负的样本被置为 `NaN` 后排除。排除后的统计量是条件样本的统计量。

**2. 标准差不是标准误。**

这里要的是输入传播后的分散程度。不要除以 $\sqrt{n}$。

**3. 非对称误差要报两个数。**

$\Delta_-$ 和 $\Delta_+$ 分别相对中位数报告。

**4. 相关性没有依据就不要引入。**

假设独立是可以接受的，但必须写明。设了相关矩阵就必须说明相关系数的来源。

**5. 平方和合成的前提是独立且近似线性。**

单参数结果按平方和合成前，先检查响应是否近似线性。若某个参数使输出跨越符号变化或阈值，平方和不成立。

**6. $S_H$ 与 $\bar m_c$ 的联动待确认。**

当前 `S_H` 在三个通道里都是硬编码常数，不随 `mbar_c` 变化。若 $S_H$ 的物理定义是 $4m_c^2$，那么改变 $\bar m_c$ 时它也应当移动。这一点需要先确认，再决定是否把 $S_H$ 当作独立输入。
