# 平台参数误差

本文件说明粲偶素三通道反问题中，平台参数的不确定度如何估计。方法与参考论文一致，参数名称按本仓库的 Borel 表述改写。

> **本文只记录方法。** 文中的具体数值和区间都是示例，不代表本项目的实际选择。平台边界和扫描范围要由具体计算确定，确定之后再填进报告。

## 一、参数清单

本问题的平台参数有三个，单位都是 GeV²。

| 参数 | 代码名 | 含义 |
|---|---|---|
| Borel 下界 | `M2_min` | Borel 窗口下端 |
| Borel 上界 | `M2_max` | Borel 窗口上端 |
| 连续谱分界 | `Lambda` | 非微扰区与 OPE 区的分界，同时是 $s$ 网格上端 |

与参考论文的对应关系：

| 参考论文 | 本工作 | 说明 |
|---|---|---|
| $q^2_{\rm up}$ | `M2_max` | 输入区间上端 |
| $q^2_{\rm low}$ | `M2_min` | 输入区间下端 |
| $\Lambda$ | `Lambda` | 对应 |
| $\sigma_0$（减法点） | 无 | Borel 变换消掉了减法，本问题不设 $\sigma_0$ |

参考论文的四个平台参数，在本问题中只剩三个。

## 二、平台的含义

平台不是数值技巧，是物理条件的交集。每个参数都有上下两侧的约束，平台是两侧都满足的区间。

### 2.1 `M2_min`

$M^2$ 越小，Borel 权重 $e^{-s/M^2}$ 对低 $s$ 越敏感，高激发态和连续谱的贡献被压得越低。但 $M^2$ 太小，算符乘积展开的高阶项相对变大，OPE 本身不可靠。

所以下界由**OPE 收敛性**决定。

### 2.2 `M2_max`

$M^2$ 越大，OPE 收敛越好，但 Borel 权重对大 $s$ 的压制变弱，高激发态和连续谱重新进入，基态不再主导。

所以上界由**基态主导性**决定。

两端合起来给出一条窗口 $[M^2_{\min}, M^2_{\max}]$。窗口内部的每个 $M^2$ 都同时满足两个条件，窗口端点在平台内移动时结果不应明显变化。

### 2.3 `Lambda`

$\Lambda$ 是非微扰区与 OPE 区的分界。$s<\Lambda$ 的谱由反问题求解，$s>\Lambda$ 直接取 OPE。

- $\Lambda$ 太小：可求解的区间窄，信息量不足，解对边界条件敏感。
- $\Lambda$ 太大：OPE 在接近 $\Lambda$ 处已经不成立，边界条件本身有偏。

所以 $\Lambda$ 的平台由两端同时决定。

## 三、误差定义

对每个平台参数 $p$，在它的平台 $\mathcal{P}$ 内扫描，记录输出 $y$ 的最大值和最小值，取

$$
\Delta_{\rm plat}(p) \;=\; \frac{\max_{p\in\mathcal{P}} y(p) \;-\; \min_{p\in\mathcal{P}} y(p)}{2}
$$

平台参数按**均匀**处理，不用高斯。原因：它们没有可作中心的测量值，只有合法区间。具体做法有两种，都是均匀：在平台内取均匀网格逐点扫，或在平台内均匀随机抽样。

平台的分布类型不随单参数或联合模式改变。

$y$ 可以是基态衰变常数、激发态衰变常数，或谱泛函。每一类输出分别算一行。

## 四、三条约定

**第一，取的是输出的极差，不是参数区间的宽度。**

$\Lambda$ 的平台可能是 $[20,30]$ GeV²（示例数值），区间宽 10 GeV²；但 $f_{J/\psi}$ 在其中可能只变 2 MeV。取的是后者的极差再除二，不是 $10/2$。

**第二，每次只动一个参数，其余固定在基准值。**

扫描 `Lambda` 时，`M2_min` 和 `M2_max` 保持基准值不变。扫描 `M2_min` 时，`M2_max` 和 `Lambda` 保持基准值不变。

这样做的代价是忽略参数之间的相关性。二维窗口的两个端点实际上可以同时移动，单独扫描只覆盖矩形区域的两条边，不覆盖角点。若需要更严格的结果，见第七节。

**第三，平台边界必须有依据，且不能事后调整。**

边界由第二节的物理条件决定。不能为了把误差做小而收窄区间，也不能为了显得稳健而放宽区间。选定边界所依据的判据要写进报告。

## 五、步骤

整体流程：

```text
PROCEDURE platform_error(p, P_range, base, n_scan):
    y ← empty list
    FOR v IN linspace(P_range, n_scan):
        args ← base with p ← v          # 其余两个参数保持基准值
        (f1, f2, alpha, flags) ← SOLVE(args)
        IF alpha on boundary OR f1 not finite:
            record (v, reason); continue
        y.append((v, f1, f2))
    IF count(y) < 2:
        RETURN undefined
    Δ(p) ← (max(y) − min(y)) / 2
    RETURN Δ(p), count(y)
```

每个参数调用一次，得到三个 $\Delta$。下面五步展开这一流程。

### 5.1 第一步：确定基准点和平台

先在三通道各自的求解程序里找到基准窗口。基准点取平台中心，不是平台端点。

```python
base = dict(M2_min=M2_MIN_CENTER, M2_max=M2_MAX_CENTER, Lambda=LAMBDA_CENTER)
```

设置好三个参数各自的扫描区间。区间必须落在平台内，并且要有物理依据。

### 5.2 第二步：单参数扫描

每个扫描点独立调用一次 `solve()`。$\alpha$ 由程序内部自动重选，扫描代码不需要处理。

```python
import numpy as np
from Jpsi_vector import solve

def scan_parameter(name, values, base):
    """扫描一个平台参数，其余两个固定在 base 中。"""
    rows = []
    for v in values:
        args = dict(base)
        args[name] = v
        res = solve(
            args["M2_min"], args["M2_max"], args["Lambda"],
            channel="2pole", mu=2.0, rg_loops=4,
        )
        ch = res["channels"]["2pole"]
        rows.append(dict(
            value=v,
            f1=float(ch["f"][0]),
            f2=float(ch["f"][1]),
            alpha=float(ch["alpha_star"]),
            on_boundary=bool(ch["alpha_diagnostics"]["on_boundary"]),
            rel_res=float(ch["relative_weighted_residual"]),
            neg_rho=int(ch["negative_rho_count"]),
        ))
    return rows
```

调用示例：

```python
rows = scan_parameter("Lambda", np.arange(LAMBDA_LO, LAMBDA_HI + 1e-9, LAMBDA_STEP), base)
```

`eta_c_axial_vector` 和 `eta_c_pseudoscalar` 的接口相同，只是导入路径不同。赝标量通道注意区分 `f` 与 `f_raw`：`f` 已经包含轴反常修正，`f_raw` 不含。同一张误差表里只能用同一种定义。

### 5.3 第三步：计算极差

```python
def plateau_error(rows, key):
    values = np.array([r[key] for r in rows], float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return np.nan, values.size
    return 0.5 * (values.max() - values.min()), values.size
```

`np.isfinite` 过滤掉极点残量为负导致的 `NaN`。**过滤后的样本数必须记录**，见第九节。

### 5.4 第四步：检查

每个扫描点都要检查下面三项。任何一项不合格，该点不能进入极差统计。

```python
ok = [
    r for r in rows
    if not r["on_boundary"]
    and np.isfinite(r["f1"])
    and r["rel_res"] < RESIDUAL_TOL
]
```

`RESIDUAL_TOL` 取基准点残差的若干倍，例如 3 倍。

### 5.5 第五步：记录

一次扫描至少记录下面这些量，写进结果表：

| 字段 | 含义 |
|---|---|
| 参数名与基准值 | 扫描的是谁，其余两个参数取什么值 |
| 扫描区间与点数 | 平台边界及其依据 |
| 每个扫描点的 $\alpha_{\rm opt}$ | 用于判断是否落在稳定窗口内 |
| `on_boundary` | $\alpha$ 是否落在候选网格端点 |
| `relative_weighted_residual` | 数据拟合程度 |
| 有效样本数 | 过滤后剩几个点 |
| $\Delta_{\rm plat}$ | 最终误差 |

## 六、检查项

### 6.1 $\alpha$ 是否落在端点

`solve()` 返回的 `alpha_diagnostics.on_boundary` 为真时，说明该点的 $\alpha$ 选到了候选序列的两端。这个点不可用。

跨点比较时这一点尤其重要：不同 $\Lambda$ 对应不同的求解区间，$\alpha$ 的稳定窗口也会移动。如果某些点跑出窗口，说明这次扫描超出了方法的适用范围。

### 6.2 残差是否稳定

`relative_weighted_residual` 应在整个扫描区间内量级一致。若某点突然上升一个量级，说明该点的 OPE 输入与解的匹配出了问题。

### 6.3 解是否非物理

`negative_rho_count` 统计重建谱中为负的网格点数。少量负值可接受，随参数变化应保持同一量级。若某个参数取值使负值数量骤增，说明该点不在平台内。

### 6.4 是否单调

平台内输出应接近平坦。若输出随参数单调变化且斜率明显，说明扫描区间没有完全落在平台内，或者平台判定有误。

画一张 $y$ 对 $p$ 的图，比只看数字更容易发现这个问题。参考论文的 Fig. 3 就是这个用途。

### 6.5 极点残量是否为负

残量为负时 `f` 被置为 `NaN`，不是零，也不是绝对值。这会影响极差统计，必须记录被排除的样本数。

## 七、二维窗口的处理

`M2_min` 和 `M2_max` 共同定义一个二维窗口，两个端点可以同时移动。第五节的单参数扫描只覆盖矩形的两条边，不覆盖角点。

三种处理方式，按严格程度排列：

**方式一：单参数扫描（参考论文做法）。**

每条边单独扫，两个 $\Delta$ 最后按平方和合成。计算量最小，缺点是忽略角点，且隐含假设两个端点独立。

**方式二：二维网格。**

在矩形内取网格，每个点都算一次。

```python
for m_lo in np.arange(lo_plateau[0], lo_plateau[1] + 1e-9, 1.0):
    for m_hi in np.arange(hi_plateau[0], hi_plateau[1] + 1e-9, 1.0):
        if m_hi <= m_lo:
            continue
        ...
```

最终误差取二维极差的一半：

$$
\Delta_{\rm window} = \frac{\max y - \min y}{2}
$$

这个量已经包含两个端点的联合效应，**不要再和单参数结果相加**，否则重复计算。

**方式三：平台内随机抽样。**

在矩形内均匀抽若干个点，取极差。点少时极差会低估真实跨度，需要足够多的样本。

三种方式的结论应互相接近。若方式二明显大于方式一，说明角点效应显著，报告中必须说明采用了哪一种。

## 八、报告格式

每个通道单独一张表。

| 参数 | 平台 | 依据 | $\Delta f_1$ [MeV] | $\Delta f_2$ [MeV] | 有效点数 |
|---|---|---|---|---|---|
| `M2_min` | $[\ldots,\ldots]$ GeV² | OPE 收敛 | — | — | — |
| `M2_max` | $[\ldots,\ldots]$ GeV² | 基态主导 | — | — | — |
| `Lambda` | $[\ldots,\ldots]$ GeV² | OPE 匹配 | — | — | — |

空白处填实际算出的数。表下方必须写明：

- 基准点 $(M^2_{\min}, M^2_{\max}, \Lambda)$。
- $\alpha$ 在每个点重选，说明选择准则。
- 被排除的点及其原因。
- 采用第七节的哪一种窗口处理方式。

三通道（`Jpsi_vector`、`eta_c_axial_vector`、`eta_c_pseudoscalar`）分别出表，不能合并成一行。三者的极点模型和残量转换不同，误差不能互相替代。

## 九、注意事项

**1. 参数区间宽度不是误差。** 见第四节第一条。这是最容易犯的错误。

**2. 平台边界不能事后调。** 先在物理判据下定边界，再算误差。反过来做会让误差失去意义。

**3. 过滤 `NaN` 会引入选择偏差。** 残量为负的点被排除后，剩下的是条件样本。有效点数和原始点数都要报。

**4. $\alpha$ 不在本文件的范围内。** 本文件只处理平台参数，$\alpha$ 的误差另议。

**5. 两个 $\Delta$ 不能重复计入。** 若用二维网格直接得到联合跨度，就不要再把单参数结果加进去。

**6. 扫描步长影响极差。** 极差是 $\max-\min$，步长越粗越容易漏掉真实极值。步长要足够密，并且报告中要写明。
