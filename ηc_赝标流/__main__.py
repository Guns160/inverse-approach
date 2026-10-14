"""Three-position-argument CLI; unique result directory, pickle-free NPZ."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import tempfile
import numpy as np
from .solver import solve


def _flatten(value, prefix=""):
    out = {}
    for key, item in value.items():
        name = f"{prefix}__{key}" if prefix else key
        if isinstance(item, dict):
            out.update(_flatten(item, name))
        elif item is not None:
            out[name] = np.asarray(item)
    return out


def _json_finite(value):
    if isinstance(value, dict):
        return {k: _json_finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_finite(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_result(result, output_root):
    """Create a unique result subdirectory; never silently overwrite files."""
    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = Path(tempfile.mkdtemp(prefix=f"eta_c_pp_{stamp}_", dir=root))
    with (dest / "result.npz").open("xb") as stream:
        np.savez_compressed(stream, **_flatten(result))
    summary = dict(result_directory=str(dest), metadata=result["metadata"], n_borel=len(result["M2"]), channels={})
    for channel, r in result["channels"].items():
        d = r["alpha_diagnostics"]
        summary["channels"][channel] = dict(alpha=r["alpha_star"], residues=r["residues"].tolist(),
            f_GeV=r["f"].tolist(), weighted_residual=r["weighted_residual"],
            relative_weighted_residual=r["relative_weighted_residual"], source=d["source"],
            fallback=d["fallback"], on_boundary=d["on_boundary"])
    summary = _json_finite(summary)
    with (dest / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description="eta_c pseudoscalar-current (P-P) inversion, mu=2 GeV (GeV^2 inputs)")
    parser.add_argument("M2_min", type=float)
    parser.add_argument("M2_max", type=float)
    parser.add_argument("Lambda", type=float)
    parser.add_argument("--channel", choices=("1pole", "2pole"), default=None)
    parser.add_argument("--ds", type=float, default=0.1,
                        help="spectral hat step; 0.01 reproduces the accepted production chain")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    try:
        result = solve(args.M2_min, args.M2_max, args.Lambda, channel=args.channel, ds=args.ds)
        summary = save_result(result, args.output)
    except (ValueError, RuntimeError, FloatingPointError, OSError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(summary, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
