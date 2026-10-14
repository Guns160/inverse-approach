"""Configurable workflow around a user/agent-provided evaluate(parameters).

The core has no particle names, condensate definitions, kinematic variables,
platform-finding rules or fixed number of observables. Model files are ordinary
Python programs explicitly selected by the caller, not mathematical strings
passed to eval().
"""
from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

from .sampling import normalize_output, run_sampling

CONFIG_KEYS = {"model", "central", "variables", "dependence",
               "provenance", "notes"}


def _json_safe(value):
    if isinstance(value, np.ndarray) and value.ndim == 0:
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def load_evaluator(model, *, base_directory="."):
    """Load the explicitly selected local model and bind its fixed options."""
    if not isinstance(model, dict):
        raise ValueError("model must identify a Python file or importable module and a function")
    allowed = {"file", "module", "function", "options"}
    if set(model) - allowed:
        raise ValueError(f"unknown model keys: {sorted(set(model) - allowed)}")
    if ("file" in model) == ("module" in model):
        raise ValueError("model must specify exactly one of file or module")
    function = model.get("function", "evaluate")
    if not isinstance(function, str) or not function.isidentifier():
        raise ValueError("model.function must be a function name")
    options = model.get("options", {})
    if not isinstance(options, dict):
        raise ValueError("model.options must be an object of fixed keyword arguments")
    provenance = {"function": function, "options": options}
    if "file" in model:
        path = Path(model["file"]).expanduser()
        if not path.is_absolute():
            path = Path(base_directory) / path
        path = path.resolve()
        if not path.is_file() or path.suffix != ".py":
            raise ValueError("model.file must point to an existing Python source file")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        name = "_uncertainty_model_" + digest[:20]
        spec = importlib.util.spec_from_file_location(name, path)
        if spec is None or spec.loader is None:
            raise ValueError(f"cannot load model source {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        previous_path = list(sys.path)
        try:
            sys.path.insert(0, str(path.parent))
            spec.loader.exec_module(module)
        finally:
            sys.path[:] = previous_path
        provenance.update(file=str(path), sha256=digest)
    else:
        module = importlib.import_module(model["module"])
        provenance["module"] = model["module"]
        source = getattr(module, "__file__", None)
        if source and Path(source).is_file():
            provenance["sha256"] = hashlib.sha256(Path(source).read_bytes()).hexdigest()
    callback = getattr(module, function, None)
    if not callable(callback):
        raise ValueError(f"{function} is not callable in the selected model")

    def evaluate(parameters):
        previous_path = list(sys.path)
        try:
            if "file" in provenance:
                sys.path.insert(0, str(Path(provenance["file"]).parent))
            return callback(dict(parameters), **copy.deepcopy(options))
        finally:
            sys.path[:] = previous_path

    return evaluate, provenance


def analyze(config, *, mode, n, parameter=None, seed=0, supplied_draws=None,
            evaluator=None, base_directory="."):
    """Evaluate one source's marginal draws or all configured joint draws."""
    if not isinstance(config, dict) or set(config) - CONFIG_KEYS:
        raise ValueError("configuration must contain model, central, variables, dependence and optional provenance/notes only")
    central = config.get("central")
    variables = config.get("variables")
    if evaluator is None:
        evaluator, model_info = load_evaluator(config.get("model"), base_directory=base_directory)
    else:
        if not callable(evaluator):
            raise ValueError("evaluator must be callable")
        model_info = {"source": "caller_supplied_callable"}
    # A bad central point is reported before sampling, not silently replaced.
    if not isinstance(central, dict):
        raise ValueError("central must be an object of independent numeric inputs")
    baseline = evaluator(dict(central))
    output_names, baseline_array = normalize_output(baseline)
    if not np.isfinite(baseline_array).all():
        raise ValueError("the model must return finite scalar/vector/named-scalar results at the central point")
    sampled = run_sampling(
        evaluator, central, variables, n=n, mode=mode, parameter=parameter,
        seed=seed, dependence=config.get("dependence"), supplied_draws=supplied_draws,
        output_names=output_names,
    )
    result = dict(sampled)
    result.update({
        "schema": "general-monte-carlo-v1",
        "central_parameters": dict(central), "central_result": _json_safe(baseline),
        "model": model_info, "configuration": config,
        "interpretation": ("one_input_marginal_random_sampling_others_fixed"
                           if mode == "single" else "joint_random_sampling_of_all_configured_inputs"),
        "software": {"python": sys.version.split()[0], "numpy": np.__version__},
    })
    return result


def save_analysis(result, output_root):
    """Write raw paired draws, model results, distribution and configuration."""
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = Path(tempfile.mkdtemp(prefix=f"mc_{result['mode']}_{stamp}_", dir=root))
    inputs = list(result["central_parameters"])
    outputs = list(result["output_names"])
    with (dest / "samples.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["draw_index", "valid", "reason",
                         *["input:" + name for name in inputs],
                         *["output:" + name for name in outputs]])
        for row in result["samples"]:
            values = row.get("values") or [None] * len(outputs)
            writer.writerow([row["index"], row["valid"], row.get("reason", ""),
                             *[row["inputs"][name] for name in inputs], *values])
    with (dest / "distribution.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["observable", "bin_left", "bin_right", "count"])
        valid = [row["values"] for row in result["samples"] if row["valid"]]
        if valid:
            values = np.asarray(valid, float)
            for j, name in enumerate(outputs):
                counts, edges = np.histogram(values[:, j], bins=min(50, max(1, int(np.sqrt(len(values))))))
                writer.writerows((name, left, right, int(count))
                                 for count, left, right in zip(counts, edges[:-1], edges[1:]))
    (dest / "settings.json").write_text(json.dumps(_json_safe(result["configuration"]),
        ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    summary = {key: value for key, value in result.items() if key != "samples"}
    summary.update(result_directory=str(dest),
                   samples_sha256=hashlib.sha256((dest / "samples.csv").read_bytes()).hexdigest())
    summary = _json_safe(summary)
    (dest / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="General single-parameter and joint Monte Carlo uncertainty workflow")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--mode", choices=("single", "joint"), required=True)
    parser.add_argument("--parameter", help="one input name, required for single mode")
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--draws", type=Path, help="CSV with pre-generated paired joint input draws")
    parser.add_argument("--output", type=Path, default=Path("results") / "uncertainty")
    args = parser.parse_args(argv)
    try:
        config = json.loads(args.config.read_text())
        supplied = None
        draw_source = None
        if args.draws:
            with args.draws.open(newline="", encoding="utf-8-sig") as stream:
                reader = csv.DictReader(stream)
                if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
                    raise ValueError("draw CSV must have distinct named columns")
                supplied = [{key: float(value) for key, value in row.items()} for row in reader]
            draw_source = {"path": str(args.draws.resolve()),
                           "sha256": hashlib.sha256(args.draws.read_bytes()).hexdigest()}
        result = analyze(config, mode=args.mode, n=args.samples, parameter=args.parameter,
                         seed=args.seed, supplied_draws=supplied, base_directory=args.config.parent)
        result["supplied_draw_source"] = draw_source
        summary = save_analysis(result, args.output)
    except (ValueError, TypeError, KeyError, RuntimeError, OSError, ImportError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(summary, ensure_ascii=False, allow_nan=False))
    return 0 if summary["summary"]["n_valid"] >= 2 else 1
