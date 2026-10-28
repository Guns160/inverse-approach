"""Shared CLI entry: unique result directory and pickle-free NPZ."""
from __future__ import annotations
import argparse
from functools import partial
from .solver import solve, DS_DEFAULT, DEFAULT_ANOMALY_1S
from inverse_common.cli import (
    add_common_arguments, flatten as _flatten, json_finite as _json_finite,
    run_and_print, save_result as _save_result,
)


def save_result(result, output_root):
    return _save_result(result, output_root, prefix="eta_c_pp")


def _optional_anomaly(value):
    if value.lower() == "none":
        return None
    return float(value)


def main():
    parser = argparse.ArgumentParser(
        description="eta_c pseudoscalar-current (P-P) inversion, mu=2 GeV (GeV^2 inputs)")
    add_common_arguments(parser, channels=("2pole",), ds_default=DS_DEFAULT)
    parser.add_argument("--anomaly-1s", type=_optional_anomaly, default=DEFAULT_ANOMALY_1S,
                        help="eta_c(1S) anomaly matrix element in GeV^3 (default -0.079); none omits it")
    parser.add_argument("--anomaly-2s", type=_optional_anomaly, default=None,
                        help="eta_c(2S) anomaly matrix element in GeV^3 (default none)")
    args = parser.parse_args()
    solve_with_anomaly = partial(solve, anomaly_1s=args.anomaly_1s,
                                anomaly_2s=args.anomaly_2s)
    run_and_print(parser, solve_with_anomaly, args, prefix="eta_c_pp")


if __name__ == "__main__":
    main()
