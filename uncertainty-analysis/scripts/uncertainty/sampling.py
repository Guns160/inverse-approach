from __future__ import annotations

from collections.abc import Mapping
import math

import numpy as np

COVERAGE = math.erf(1 / math.sqrt(2))
ENSEMBLES = ("statistical", "systematic")


class OutputSchemaError(ValueError):
    """A model changed the meanings or dimension of its outputs."""


def _number(value, label):
    if isinstance(value, (bool, np.bool_, str, bytes)) or not np.isscalar(value):
        raise ValueError(f"{label} must be a finite real scalar")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a finite real scalar") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be a finite real scalar")
    return result


def _names(names, label):
    if not names or any(not isinstance(x, str) or not x for x in names) or len(set(names)) != len(names):
        raise ValueError(f"{label} must contain distinct nonempty names")
    return list(names)


def _prepare_constraints(constraints, center):
    if constraints is None:
        return []
    if not isinstance(constraints, list) or not constraints:
        raise ValueError("constraints must be a nonempty list when supplied")
    prepared = []
    for i, item in enumerate(constraints):
        if not isinstance(item, Mapping):
            raise ValueError(f"constraint {i} must be an object")
        kind = item.get("kind")
        if kind == "positive":
            names = _names(item.get("names"), f"constraint {i} names")
            if set(names) - set(center):
                raise ValueError(f"constraint {i} names are not central parameters")
            prepared.append((kind, names))
        elif kind in ("less_than_margin", "greater_than_margin"):
            left, right, margin_name = item.get("left"), item.get("right"), item.get("margin_name")
            if any(not isinstance(x, str) or not x for x in (left, right, margin_name)):
                raise ValueError(f"constraint {i} requires string left, right, and margin_name")
            if set((left, right, margin_name)) - set(center):
                raise ValueError(f"constraint {i} refers to a parameter absent from central")
            if _number(center[margin_name], margin_name) < 0:
                raise ValueError("constraint margin must be nonnegative")
            prepared.append((kind, left, right, margin_name))
        elif kind == "min_difference":
            lower, upper = item.get("lower"), item.get("upper")
            minimum = _number(item.get("minimum"), f"constraint {i}.minimum")
            if any(not isinstance(x, str) or not x for x in (lower, upper)):
                raise ValueError(f"constraint {i} requires string lower and upper")
            if set((lower, upper)) - set(center):
                raise ValueError(f"constraint {i} refers to a parameter absent from central")
            prepared.append((kind, lower, upper, minimum))
        else:
            raise ValueError(f"unsupported constraint kind: {kind!r}")
    return prepared


def _constraint_ok(row, constraints):
    for constraint in constraints:
        kind = constraint[0]
        if kind == "positive" and any(row[name] <= 0 for name in constraint[1]):
            return False
        if kind == "less_than_margin":
            _, left, right, margin_name = constraint
            if not row[left] < row[right] - row[margin_name]:
                return False
        if kind == "greater_than_margin":
            _, left, right, margin_name = constraint
            if not row[left] > row[right] + row[margin_name]:
                return False
        if kind == "min_difference":
            _, lower, upper, minimum = constraint
            if not row[upper] - row[lower] > minimum:
                return False
    return True


def _prepare(central, variables, n, ensemble, seed, supplied):
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("n must be an integer >= 2")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    if ensemble not in ENSEMBLES:
        raise ValueError("ensemble must be 'statistical' or 'systematic'")
    if not isinstance(central, Mapping) or not central:
        raise ValueError("central must be a nonempty parameter mapping")
    _names(list(central), "central parameter names")
    center = {name: _number(value, name) for name, value in central.items()}
    if not isinstance(variables, Mapping) or not variables:
        raise ValueError("variables must define at least one sampled input")
    _names(list(variables), "variable names")
    if set(variables) - set(center):
        raise ValueError("all variable names must exist in central")
    specs = {}
    metadata = {"group", "unit", "source", "provenance", "notes"}
    for name, item in variables.items():
        if not isinstance(item, Mapping):
            raise ValueError(f"variable {name} must have a distribution object")
        spec = dict(item)
        kind = spec.get("kind")
        if kind == "normal":
            allowed = {"kind", "sigma"} | metadata
            spec["sigma"] = _number(spec.get("sigma"), name + ".sigma")
            if spec["sigma"] < 0:
                raise ValueError("normal sigma must be nonnegative")
        elif kind in ("uniform", "loguniform"):
            allowed = {"kind", "low", "high"} | metadata
            spec["low"] = _number(spec.get("low"), name + ".low")
            spec["high"] = _number(spec.get("high"), name + ".high")
            if spec["low"] > spec["high"]:
                raise ValueError("distribution low must not exceed high")
            if kind == "loguniform" and spec["low"] <= 0:
                raise ValueError("loguniform bounds must be positive")
            if not spec["low"] <= center[name] <= spec["high"]:
                raise ValueError(f"central {name} is outside its distribution bounds")
        elif kind == "provided" and supplied:
            allowed = {"kind"} | metadata
        else:
            raise ValueError(f"unsupported distribution for {name}: {kind!r}")
        if set(spec) - allowed:
            raise ValueError(f"unknown distribution fields for {name}: {sorted(set(spec) - allowed)}")
        specs[name] = spec
    return center, specs


def _marginal(rng, center, spec, n):
    kind = spec["kind"]
    if kind == "normal":
        return rng.normal(center, spec["sigma"], n)
    if kind == "uniform":
        return rng.uniform(spec["low"], spec["high"], n)
    if kind == "loguniform":
        return np.exp(rng.uniform(math.log(spec["low"]), math.log(spec["high"]), n))
    raise ValueError(f"cannot generate marginal distribution {kind!r}")


def _normal_block(center, specs, dependence, rng, n):
    names = _names(dependence.get("names", []), "covariance names")
    if any(name not in specs or specs[name]["kind"] != "normal" for name in names):
        raise ValueError("covariance names must identify configured normal variables")
    remaining = [name for name in specs if name not in names]
    if remaining and dependence.get("independent_remaining") is not True:
        raise ValueError("declare independent_remaining=true for variables outside the normal block")
    matrix = np.asarray(dependence.get("matrix"), dtype=float)
    if matrix.shape != (len(names), len(names)) or not np.isfinite(matrix).all():
        raise ValueError("covariance matrix has invalid shape or nonfinite entries")
    if not np.allclose(matrix, matrix.T, rtol=1e-10, atol=0):
        raise ValueError("covariance matrix must be symmetric")
    sigmas = np.array([specs[name]["sigma"] for name in names])
    if not np.allclose(np.diag(matrix), sigmas**2, rtol=1e-10, atol=0):
        raise ValueError("covariance diagonal must equal the declared marginal sigma squared")
    nonzero = sigmas > 0
    if np.any(matrix[~nonzero] != 0) or np.any(matrix[:, ~nonzero] != 0):
        raise ValueError("zero-variance inputs must have zero covariances")
    draws = np.tile([center[name] for name in names], (n, 1)).astype(float)
    if nonzero.any():
        scale = sigmas[nonzero]
        corr = matrix[np.ix_(nonzero, nonzero)] / np.outer(scale, scale)
        eigenvalues, eigenvectors = np.linalg.eigh((corr + corr.T) / 2)
        if float(eigenvalues.min()) < -1e-10:
            raise ValueError("covariance matrix must be positive semidefinite")
        factor = eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))[None, :]
        correlated = np.dot(rng.standard_normal((n, len(scale))), factor.T)
        draws[:, nonzero] += correlated * scale[None, :]
    return names, draws, remaining


def _candidate_rows(center, specs, dependence, rng, n):
    if not isinstance(dependence, Mapping):
        raise ValueError("each ensemble requires an explicit dependence model")
    names = list(specs)
    columns = {}
    kind = dependence.get("kind")
    allowed = {"kind", "source", "notes"}
    if kind == "normal_covariance":
        allowed |= {"names", "matrix", "independent_remaining"}
    if set(dependence) - allowed:
        raise ValueError(f"unknown dependence fields: {sorted(set(dependence) - allowed)}")
    if kind == "independent":
        remaining = names
    elif kind == "normal_covariance":
        block_names, block, remaining = _normal_block(center, specs, dependence, rng, n)
        columns.update({name: block[:, j] for j, name in enumerate(block_names)})
    else:
        raise ValueError(f"unsupported ensemble dependence: {kind!r}")
    for name in remaining:
        columns[name] = _marginal(rng, center[name], specs[name], n)
    return [dict(center, **{name: float(columns[name][i]) for name in names}) for i in range(n)]


def generate_draws(central, variables, *, n, ensemble, seed=0,
                   dependence=None, constraints=None, max_draw_attempts=200000,
                   supplied_draws=None):
    """Generate accepted input rows for one separated ensemble."""
    center, specs = _prepare(central, variables, n, ensemble, seed, supplied_draws is not None)
    prepared_constraints = _prepare_constraints(constraints, center)
    n, seed = int(n), int(seed)
    if isinstance(max_draw_attempts, (bool, np.bool_)) or not isinstance(max_draw_attempts, (int, np.integer)):
        raise ValueError("max_draw_attempts must be an integer")
    max_draw_attempts = int(max_draw_attempts)
    if max_draw_attempts < n:
        raise ValueError("max_draw_attempts must be at least n")
    names = list(specs)
    design = {"ensemble": ensemble, "n": n, "seed": seed,
              "variable_names": names, "active_parameters": names,
              "distributions": specs, "dependence": dependence,
              "constraints": constraints,
              "fixed_parameters": [name for name in center if name not in names],
              "resampling_for_physical_inputs": bool(prepared_constraints)}
    rows, attempted = [], 0
    if supplied_draws is not None:
        if dependence is not None and dependence.get("kind") != "provided_paired_rows":
            raise ValueError("provided draws require dependence=null or kind=provided_paired_rows")
        if len(supplied_draws) != n:
            raise ValueError("supplied draw count must equal n")
        for i, row in enumerate(supplied_draws):
            if not isinstance(row, Mapping) or set(row) != set(specs):
                raise ValueError("each supplied row must contain exactly the configured variable names")
            candidate = dict(center, **{name: _number(row[name], f"draw {i}: {name}") for name in names})
            if prepared_constraints and not _constraint_ok(candidate, prepared_constraints):
                raise ValueError("provided rows must already satisfy the declared constraints")
            rows.append(candidate)
        design.update(sampling_source="provided_paired_rows", seed_used=False,
                      dependence_used="provided_paired_rows")
    else:
        rng = np.random.default_rng(seed)
        while len(rows) < n:
            batch = min(max(50, 2 * (n - len(rows))), max_draw_attempts - attempted)
            if batch <= 0:
                break
            candidates = _candidate_rows(center, specs, dependence, rng, batch)
            for candidate in candidates:
                if len(rows) == n:
                    break
                if not prepared_constraints or _constraint_ok(candidate, prepared_constraints):
                    rows.append(candidate)
                # Only candidates actually considered before filling the requested
                # sample count belong to the attempted-input count.
                attempted += 1
        design.update(sampling_source="numpy_default_rng", seed_used=True,
                      dependence_used=dependence.get("kind"))
        if len(rows) < n:
            raise ValueError("unable to draw enough physical samples before max_draw_attempts")
    design.update(n_candidate_inputs_attempted=attempted, n_accepted_inputs=len(rows),
                  input_rejection_fraction=(attempted - len(rows)) / attempted if attempted else 0.0)
    matrix = np.array([[row[name] for name in names] for row in rows])
    design["sample_input_covariance_names"] = names
    design["sample_input_covariance"] = np.atleast_2d(np.cov(matrix, rowvar=False, ddof=1)).tolist()
    return {"draws": rows, "design": design}


def normalize_output(value, output_names=None):
    """Return stable names and a real vector; nonfinite entries remain visible."""
    if isinstance(value, Mapping):
        names = _names(list(value), "observable names")
        if output_names is not None:
            if set(names) != set(output_names):
                raise OutputSchemaError("model output names changed between draws")
            names = list(output_names)
        raw = list(value[name] for name in names)
        if any(not np.isscalar(item) for item in raw):
            raise OutputSchemaError("each named model output must be a scalar")
    else:
        raw = value
        try:
            shape = np.asarray(raw).shape
        except (TypeError, ValueError) as exc:
            raise OutputSchemaError("model output must be a scalar or a one-dimensional numeric vector") from exc
        if len(shape) > 1:
            raise OutputSchemaError("model output must be a scalar or a one-dimensional numeric vector")
        size = shape[0] if shape else 1
        names = list(output_names) if output_names is not None else (["output"] if not shape else [f"output_{i}" for i in range(size)])
    try:
        if np.iscomplexobj(raw):
            raise OutputSchemaError("model output must be real")
        array = np.asarray(raw, dtype=float).reshape(-1)
    except (TypeError, ValueError) as exc:
        raise OutputSchemaError("model output must contain real numeric scalars") from exc
    if array.size == 0 or array.size != len(names):
        raise OutputSchemaError("model output dimension changed or is empty")
    return names, array


def summarize(samples, output_names):
    valid = [row["values"] for row in samples if row["valid"]]
    n, good = len(samples), len(valid)
    invalid = n - good
    summary = {"n_requested": n, "n_valid": good, "n_invalid": invalid,
               "invalid_fraction": invalid / n, "conditional_on_valid": invalid > 0,
               "status": "insufficient_valid_draws" if good < 2 else ("conditional_on_valid" if invalid else "ok"),
               "coverage_probability": COVERAGE, "quantile_method": "linear",
               "outputs": {}, "covariance": None, "correlation": None,
               "output_covariance_names": list(output_names), "failure_counts": {}}
    for row in samples:
        if not row["valid"]:
            reason = row["reason"]
            summary["failure_counts"][reason] = summary["failure_counts"].get(reason, 0) + 1
    values = np.asarray(valid, float).reshape(good, len(output_names)) if good else np.empty((0, len(output_names)))
    points = [0.025, (1-COVERAGE)/2, 0.16, 0.5, 0.84, (1+COVERAGE)/2, 0.975]
    for j, name in enumerate(output_names):
        if not good:
            summary["outputs"][name] = {key: None for key in
                ("mean", "std", "median", "q025", "q16", "q84", "q975",
                 "central_68_interval", "minimum", "maximum")}
            continue
        q025, qlow, q16, median, q84, qhigh, q975 = np.quantile(values[:, j], points, method="linear")
        summary["outputs"][name] = {
            "mean": float(values[:, j].mean()), "std": float(values[:, j].std(ddof=1)),
            "median": float(median), "q025": float(q025), "q16": float(q16),
            "q84": float(q84), "q975": float(q975),
            "central_68_interval": [float(qlow), float(qhigh)],
            "minimum": float(values[:, j].min()), "maximum": float(values[:, j].max()),
        }
    if good > 1:
        covariance = np.atleast_2d(np.cov(values, rowvar=False, ddof=1))
        scale = np.sqrt(np.maximum(0, np.diag(covariance)))
        denom = np.outer(scale, scale)
        summary["covariance"] = covariance.tolist()
        summary["correlation"] = [[float(np.clip(covariance[i, j]/denom[i, j], -1, 1))
            if denom[i, j] > 0 else None for j in range(len(scale))] for i in range(len(scale))]
    return summary


def run_sampling(evaluate, central, variables, *, n, ensemble, seed=0,
                 dependence=None, constraints=None, max_draw_attempts=200000,
                 supplied_draws=None, output_names=None):
    """Evaluate accepted draws once; solver failures remain in the record."""
    if not callable(evaluate):
        raise ValueError("evaluate must be callable")
    names = _names(output_names, "output_names") if output_names is not None else None
    generated = generate_draws(central, variables, n=n, ensemble=ensemble, seed=seed,
                               dependence=dependence, constraints=constraints,
                               max_draw_attempts=max_draw_attempts, supplied_draws=supplied_draws)
    records = []
    for index, point in enumerate(generated["draws"]):
        row = {"index": index, "inputs": dict(point), "values": None, "valid": False, "reason": ""}
        try:
            value = evaluate(dict(point))
        except OutputSchemaError:
            raise
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
            row["reason"] = type(exc).__name__ + ": " + str(exc)
        else:
            names, vector = normalize_output(value, names)
            row["values"] = [float(x) if np.isfinite(x) else None for x in vector]
            row["valid"] = bool(np.isfinite(vector).all())
            if not row["valid"]:
                row["reason"] = "nonfinite_model_output"
        records.append(row)
    names = names or []
    for row in records:
        if row["values"] is None:
            row["values"] = [None] * len(names)
    return {"ensemble": ensemble, "output_names": names,
            "design": generated["design"], "samples": records,
            "summary": summarize(records, names)}
