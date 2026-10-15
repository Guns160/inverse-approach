# 分系综程序接口

Python 参考程序位于 `scripts/uncertainty/`，只需 Python 3.9+ 与 NumPy 1.24+。不需要安装包；`scripts/run_uncertainty.py` 可按绝对路径运行。模型由 `evaluate(parameters, **options)` 提供，可以封装现有 Python、MATLAB 或命令行求解器。

## 模型接口

```python
def evaluate(parameters, **fixed_options):
    result = existing_solver(parameters, **fixed_options)
    return {"observable_name": float(result)}
```

返回值可以是一个数、一维向量或若干具名标量。每次调用的输出名和维数必须一致。模型用 `ValueError`、`RuntimeError`、`FloatingPointError` 或线性代数异常报告物理/数值无效点；适配器编程错误会中止运行，不会混入物理失败率。

## 配置结构

顶层字段：`model`、`central`、`ensembles`、`provenance`、`notes`。`central` 包含所有独立数值输入、固定参数和约束使用的网格量。`ensembles` 只允许 `statistical` 与 `systematic` 两个子对象。

每个系综字段：

| 字段 | 含义 |
|---|---|
| `variables` | 本系综随机变化的参数及其分布 |
| `dependence` | 本系综内部联合抽样结构 |
| `constraints` | 可选物理约束；违反的候选输入被拒绝并重新抽样 |
| `n_samples` | 接受的输入样本数 |
| `seed` | 随机种子 |
| `max_draw_attempts` | 候选输入最大尝试数，默认 200000 |

支持分布：

```json
{"kind": "normal", "sigma": 0.2}
{"kind": "uniform", "low": 2.0, "high": 4.0}
{"kind": "loguniform", "low": 1e-12, "high": 1e-10}
```

正态均值、均匀中心与对数均匀范围由 `central` 决定。`group`、`unit`、`source` 等只是元数据，不改变抽样。

### 联合分布


独立声明：

```json
{"kind": "independent"}
```

相关正态块加独立剩余参数：

```json
{
  "kind": "normal_covariance",
  "names": ["input_a", "input_b"],
  "matrix": [[0.04, 0.01], [0.01, 0.01]],
  "independent_remaining": true
}
```

`matrix` 是实际输入单位下的协方差，不是相关系数；对角必须等于对应 `sigma²`，矩阵必须对称半正定。块外剩余参数必须显式声明独立，否则拒绝运行。更复杂的联合分布可由上层生成成对样本后调用底层 `generate_draws(..., supplied_draws=...)`；提供的行必须已经满足约束，程序不打乱、不重采样。

### 约束

支持四种声明式约束：

```json
[
  {"kind": "positive", "names": ["input_a", "input_b"]},
  {"kind": "less_than_margin", "left": "lower_scale", "right": "upper_scale", "margin_name": "grid_step"},
  {"kind": "greater_than_margin", "left": "upper_scale", "right": "threshold", "margin_name": "grid_step"},
  {"kind": "min_difference", "lower": "Y1", "upper": "Y2", "minimum": 0.0}
]
```

`less_than_margin` 要求 `left < right - margin_name`；`greater_than_margin` 要求 `left > right + margin_name`。这些约束用于生成物理有效输入。拒绝率记录在 `design.input_rejection_fraction`。求解阶段的失败不走这个拒绝流程，而是保留在 `samples.csv`。

## 命令行

```sh
python scripts/run_uncertainty.py --config examples/toy_config.json \
  --ensemble statistical --output results

python scripts/run_uncertainty.py --config examples/toy_config.json \
  --ensemble systematic --samples 1000 --seed 15002 --output results
```

`--samples` 和 `--seed` 只覆盖当前运行，不修改配置文件。每次运行创建独立目录，包含：

- `settings.json`：原始配置。
- `samples.csv`：每个接受输入点、配对输出、有效标志和失败原因。
- `distribution.csv`：输出直方图分箱。
- `summary.json`：中心值、模型哈希、抽样设计、拒绝率、分布统计、协方差和系综专属报告。

统计系综的 `uncertainty_reporting.kind` 是 `statistical_distribution`，每个输出含 `std`、`median`、`q16`、`q84` 及相对中心值的 68% 偏差。系统系综是 `systematic_envelope`，每个输出含包含中心值的最小/最大值和 `envelope_minus` / `envelope_plus`。没有 `total` 字段。

## 测试

```sh
python -B -W error -m unittest discover -s uncertainty-analysis/tests -v
```

测试覆盖：一个系综运行时另一类参数逐行固定、相关正态块、独立均匀参数、约束拒绝采样、种子可复现、模型失败不补抽、输出 schema 错误与适配器错误不误记为物理失败、CLI 两个系综分别运行、旧的单/联合总误差配置被拒绝。
