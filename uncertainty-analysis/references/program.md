# 通用程序接口

附带工具是 Skill 的一个参考实现；物理模型不属于抽样器。Python 3.9+ 和 NumPy 1.24+ 即可。脚本可直接按绝对路径运行，不需要安装 Skill 或 Python 包。

## 最小模型接口

```python
def evaluate(parameters, **fixed_options):
    # 调用已有求解器；必要时重算派生参数。
    # 参数名、单位和返回值含义由本问题定义。
    result = existing_solver(parameters, **fixed_options)
    return {"observable": float(result)}
```

支持一个数、一维数值向量或若干具名标量。每次调用的输出名称和维数须一致。通过 `ValueError`、`RuntimeError`、`FloatingPointError` 或 NumPy 线性代数异常报告无法求值的输入；其他编程异常会中止，以免把代码错误混为物理失效点。非有限输出按无效行记录。

模型必须真正消费被采样参数。如果原代码把输入写死在全局常数中，Agent 应先建立局部显式输入接口并验证；不要只修改输出标签。适配器可以调用 MATLAB、现有 CLI 或远程求解服务，但执行权限和任务范围仍按当前用户请求处理。

## 配置

读取 [examples/toy_config.json](../examples/toy_config.json) 和 [examples/toy_model.py](../examples/toy_model.py) 可运行一个三个输出的代数示例。示例中的相关系数和数值没有物理含义。

| 字段 | 含义 |
|---|---|
| `model.file` 或 `model.module` | 显式选择可信 Python 适配器；二者只填一个。相对文件路径以配置所在目录为准 |
| `model.function` | 函数名，默认 `evaluate` |
| `model.options` | 传给函数的固定关键字选项，例如求解精度、模型选择；不参与撒点 |
| `central` | 全部独立数值输入的中心值，包括本次保持固定的参数 |
| `variables` | 本次纳入误差分析的参数分布；名称必须在 `central` 中 |
| `dependence` | 联合输入分布的相关结构；`joint` 不默认假设独立 |
| `provenance` / `notes` | 分布来源、平台依据、输入单位、输出单位和其他分析约定 |

`variables` 支持这些边缘分布：

```json
{"kind": "normal", "sigma": 0.2, "group": "statistical"}
{"kind": "uniform", "low": 1.0, "high": 3.0, "group": "systematic"}
{"kind": "loguniform", "low": 1e-12, "high": 1e-11, "group": "systematic"}
```

正态均值取 `central` 对应值；均匀和对数均匀用显式上下界。`group` 仅为来源标签，不触发额外采样或误差合成。需要截断分布、经验后验、条件分布、离散选择或曲线平台时，由 Agent 生成联合样本表；不通过裁剪正态样本伪造分布。

### 联合相关结构

明确独立假设：

```json
{"kind": "independent"}
```

正态输入的协方差：

```json
{
  "kind": "normal_covariance",
  "names": ["input_a", "input_b"],
  "matrix": [[0.04, 0.01], [0.01, 0.01]]
}
```

矩阵顺序严格对应 `names`；矩阵是实际输入单位下的协方差，不是相关系数。要求对称、半正定，对角元与各参数 `sigma²` 一致。若存在不在该正态块内的其他抽样变量，需设置 `independent_remaining: true` 明确其与该块及彼此独立；否则拒绝运行。

对于更一般的联合分布，传入已生成的成对样本 CSV；没有内置边缘分布描述时，可将对应变量写成 `{"kind": "provided"}`。列名与 `variables` 参数名一致，行数与 `--samples` 一致，每行是一组联合输入。`central` 中其他输入保持固定。程序不独立打乱列、不重采样该表；`single` 模式只取指定参数的一列。CSV 的生成方法和关联结构由 Agent 记录在配置中。

使用样本表时，`dependence` 留空、设为 `null`，或明确设为 `{"kind": "provided_paired_rows"}`。不要同时声明另一个生成分布；此时实际联合关系由输入表确定。汇总里的协方差矩阵均带对应的名称顺序，固定列可能有零方差。

## 运行

```sh
python scripts/run_uncertainty.py --config examples/toy_config.json \
  --mode single --parameter Q2 --samples 10000 --seed 42 --output results

python scripts/run_uncertainty.py --config examples/toy_config.json \
  --mode joint --samples 10000 --seed 42 --output results

python scripts/run_uncertainty.py --config /path/to/problem.json \
  --mode joint --samples 10000 --draws /path/to/paired_inputs.csv --output results
```

以上相对路径以 Skill 目录为当前目录。也可将 `scripts/` 加入 Python 模块搜索路径后使用：

```python
from uncertainty.sampling import run_sampling

result = run_sampling(
    evaluate, central, variables,
    mode="single", parameter="arbitrary_input_name", n=1000, seed=42,
)
```

程序不自动决定样本数，也不把单参数实验自动加总为联合结果。要查看多项单来源影响，分别运行若干次 `single`；`joint` 单独运行一次完整配置。

## 输出与复核

每次运行创建唯一输出目录：

- `settings.json`：原始配置。
- `samples.csv`：所有输入样本、配对输出、有效标志和失败原因。
- `distribution.csv`：用于绘图的直方图分箱，不替代原始样本。
- `summary.json`：中心点结果、模型文件哈希、种子/采样设计、标准差、分位区间、输出协方差及有效点统计。

结果的单位由模型/配置决定，抽样器不换算单位。标准差和中央 68.268949% 区间同时保留；解释原则见 [interpretation.md](interpretation.md)。

## 检查

```sh
python -B -W error -m unittest discover -s tests -v
```

测试包括单参数保持其他输入不变、联合相关结构、任意输出维数、确定性种子、无效点记录及完整 CLI 输出。它们验证抽样工具；具体问题仍需模型自身的物理和数值验证。
