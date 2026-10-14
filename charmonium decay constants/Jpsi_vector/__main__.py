"""Shared CLI entry: unique result directory and pickle-free NPZ."""
from __future__ import annotations
import argparse
from .solver import solve
from inverse_common.cli import (
    add_common_arguments, flatten as _flatten, json_finite as _json_finite,
    run_and_print, save_result as _save_result,
)


def save_result(result, output_root):
    return _save_result(result, output_root, prefix="jpsi", include_peak_covered=True)


def main():
    parser = argparse.ArgumentParser(description="J/psi n=0 spectral inversion (GeV^2 inputs)")
    add_common_arguments(parser, channels=("1delta", "2delta"), ds_default=0.1)
    run_and_print(parser, solve, parser.parse_args(), prefix="jpsi", include_peak_covered=True)


if __name__ == "__main__":
    main()
