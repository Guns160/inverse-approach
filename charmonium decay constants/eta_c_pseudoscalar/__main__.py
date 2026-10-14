"""Shared CLI entry: unique result directory and pickle-free NPZ."""
from __future__ import annotations
import argparse
from .solver import solve, DS_DEFAULT
from inverse_common.cli import (
    add_common_arguments, flatten as _flatten, json_finite as _json_finite,
    run_and_print, save_result as _save_result,
)


def save_result(result, output_root):
    return _save_result(result, output_root, prefix="eta_c_pp")


def main():
    parser = argparse.ArgumentParser(
        description="eta_c pseudoscalar-current (P-P) inversion, mu=2 GeV (GeV^2 inputs)")
    add_common_arguments(parser, channels=("1pole", "2pole"), ds_default=DS_DEFAULT)
    run_and_print(parser, solve, parser.parse_args(), prefix="eta_c_pp")


if __name__ == "__main__":
    main()
