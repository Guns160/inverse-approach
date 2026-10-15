"""Configurable separated statistical/systematic Monte Carlo workflow.

The model is supplied by the caller. This workflow never combines the two
ensembles into a total uncertainty.
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

from .sampling import ENSEMBLES, normalize_output, run_sampling

CONFIG_KEYS = {"model", "central", "ensembles", "provenance", "notes"}
ENSEMBLE_KEYS = {"variables", "dependence", "constraints", "n_samples", "seed", "max_draw_attempts"}


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
    """Load an explicitly selected local model and bind its fixed options."""
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


def _ensemble(config, name):
    if not isinstance(config, dict) or set(config) - CONFIG_KEYS:
        raise ValueError("configuration may contain model, central, ensembles, provenance and notes only")
    groups = config.get("ensembles")
    if not isinstance(groups, dict) or set(groups) - set(ENSEMBLES):
        raise ValueError("ensembles must contain only statistical and/or systematic groups")
    if name not in ENSEMBLES or name not in groups:
        raise ValueError(f"configuration does not define the {name} ensemble")
    result = groups[name]
    if not isinstance(result, dict) or set(result) - ENSEMBLE_KEYS:
        raise ValueError(f"{name} ensemble contains unknown fields")
    for required in ("variables", "dependence", "n_samples", "seed"):
        if required not in result:
            raise ValueError(f"{name} ensemble is missing {required}")
    return result


def _reporting(ensemble, sampled, central_values):
    """Attach the reporting convention used by the corresponding ensemble."""
    summary = sampled["summary"]
    valid = [row["values"] for row in sampled["samples"] if row["valid"]]
    kind = "statistical_distribution" if ensemble == "statistical" else "systematic_envelope"
    if not valid:
        return {"kind": kind, "outputs": {name: None for name in sampled["output_names"]}}
    values = np.asarray(valid, float)
    reporting = {}
    for j, name in enumerate(sampled["output_names"]):
        output = summary["outputs"][name]
        if ensemble == "statistical":
            reporting[name] = {
                "std": output["std"], "median": output["median"],
                "q16": output["q16"], "q84": output["q84"],
                "q68_minus_from_central": float(central_values[j] - output["q16"]),
                "q68_plus_from_central": float(output["q84"] - central_values[j]),
            }
        else:
            low = min(float(central_values[j]), output["minimum"])
            high = max(float(central_values[j]), output["maximum"])
            reporting[name] = {
                "minimum_including_central": low,
                "maximum_including_central": high,
                "envelope_minus": float(central_values[j] - low),
                "envelope_plus": float(high - central_values[j]),
            }
    return {"kind": kind, "outputs": reporting}


def analyze(config, *, ensemble, n=None, seed=None, evaluator=None, base_directory="."):
    """Run exactly one statistical or systematic ensemble."""
    if ensemble not in ENSEMBLES:
        raise ValueError("ensemble must be 'statistical' or 'systematic'")
    selected = _ensemble(config, ensemble)
    central = config.get("central")
    if not isinstance(central, dict):
        raise ValueError("central must be an object of numeric model inputs")
    if evaluator is None:
        evaluator, model_info = load_evaluator(config.get("model"), base_directory=base_directory)
    elif callable(evaluator):
        model_info = {"source": "caller_supplied_callable"}
    else:
        raise ValueError("evaluator must be callable")
    baseline = evaluator(dict(central))
    output_names, baseline_array = normalize_output(baseline)
    if not np.isfinite(baseline_array).all():
        raise ValueError("the model must return finite results at the central point")
    n_samples = selected["n_samples"] if n is None else n
    selected_seed = selected["seed"] if seed is None else seed
    sampled = run_sampling(
        evaluator, central, selected["variables"], n=n_samples, ensemble=ensemble,
        seed=selected_seed, dependence=selected["dependence"],
        constraints=selected.get("constraints"),
        max_draw_attempts=selected.get("max_draw_attempts", 200000),
        output_names=output_names,
    )
    result = dict(sampled)
    result.update({
        "schema": "separated-ensemble-monte-carlo-v1",
        "central_parameters": dict(central),
        "central_result": _json_safe(baseline),
        "model": model_info,
        "configuration": config,
        "interpretation": ("statistical_inputs_random_with_systematic_parameters_fixed"
                           if ensemble == "statistical" else
                           "systematic_parameters_random_with_statistical_inputs_fixed"),
        "combination_policy": "statistical_and_systematic_are_separate; no_total_uncertainty_is_computed",
        "uncertainty_reporting": _reporting(ensemble, sampled, baseline_array.tolist()),
        "software": {"python": sys.version.split()[0], "numpy": np.__version__},
    })
    return result


def save_analysis(result, output_root):
    """Write every accepted draw, every model result and one ensemble summary."""
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = Path(tempfile.mkdtemp(prefix=f"mc_{result['ensemble']}_{stamp}_", dir=root))
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
    parser = argparse.ArgumentParser(description="Separated statistical/systematic Monte Carlo workflow")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--ensemble", choices=ENSEMBLES, required=True)
    parser.add_argument("--samples", type=int, help="override the configured n_samples")
    parser.add_argument("--seed", type=int, help="override the configured seed")
    parser.add_argument("--output", type=Path, default=Path("results") / "uncertainty")
    args = parser.parse_args(argv)
    try:
        config = json.loads(args.config.read_text())
        result = analyze(config, ensemble=args.ensemble, n=args.samples,
                         seed=args.seed, base_directory=args.config.parent)
        summary = save_analysis(result, args.output)
    except (ValueError, TypeError, KeyError, RuntimeError, OSError, ImportError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(summary, ensure_ascii=False, allow_nan=False))
    return 0 if summary["summary"]["n_valid"] >= 2 else 1
