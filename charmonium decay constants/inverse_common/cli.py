"""Shared pickle-free persistence and command-line behaviour."""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile

import numpy as np


def flatten(value, prefix=""):
    out = {}
    for key, item in value.items():
        name = f"{prefix}__{key}" if prefix else key
        if isinstance(item, dict):
            out.update(flatten(item, name))
        elif item is not None:
            out[name] = np.asarray(item)
    return out


def json_finite(value):
    if isinstance(value, dict):
        return {key: json_finite(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_finite(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_result(result, output_root, *, prefix, include_peak_covered=False):
    """Create a unique result directory; never silently overwrite files."""
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = Path(tempfile.mkdtemp(prefix=f"{prefix}_{stamp}_", dir=root))
    with (dest / "result.npz").open("xb") as stream:
        np.savez_compressed(stream, **flatten(result))
    summary = {
        "result_directory": str(dest),
        "metadata": result["metadata"],
        "n_borel": len(result["M2"]),
        "channels": {},
    }
    for channel, item in result["channels"].items():
        d = item["alpha_diagnostics"]
        channel_summary = dict(
            alpha=item["alpha_star"], residues=item["residues"].tolist(),
            f_GeV=item["f"].tolist(), weighted_residual=item["weighted_residual"],
            relative_weighted_residual=item["relative_weighted_residual"],
            source=d["source"], fallback=d["fallback"], on_boundary=d["on_boundary"])
        if include_peak_covered:
            channel_summary["small_alpha_peak_covered"] = d["small_alpha_peak_covered"]
        summary["channels"][channel] = channel_summary
    summary = json_finite(summary)
    with (dest / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return summary


def add_common_arguments(parser, *, channels, ds_default):
    parser.add_argument("M2_min", type=float)
    parser.add_argument("M2_max", type=float)
    parser.add_argument("Lambda", type=float)
    parser.add_argument("--channel", choices=channels, default=None)
    parser.add_argument("--ds", type=float, default=ds_default,
                        help="spectral hat step in GeV^2; Borel step is channel-fixed")
    parser.add_argument("--mu", type=float, default=2.0,
                        help="RG scale in GeV")
    parser.add_argument("--rg-loops", type=int, choices=(1, 2, 3, 4), default=4,
                        help="QCD RG loops used for mc and alpha_s")
    parser.add_argument("--output", type=Path)


def run_and_print(parser, solve, args, *, prefix, include_peak_covered=False):
    default_output = Path(__file__).resolve().parents[1] / "results"
    try:
        result = solve(args.M2_min, args.M2_max, args.Lambda,
                       channel=args.channel, ds=args.ds, mu=args.mu,
                       rg_loops=args.rg_loops)
        summary = save_result(result, args.output or default_output, prefix=prefix,
                              include_peak_covered=include_peak_covered)
    except (ValueError, RuntimeError, FloatingPointError, OSError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(summary, ensure_ascii=False, allow_nan=False))
