"""Shared CLI persistence for the PP + one-BW alpha family."""
from __future__ import annotations

import argparse
from functools import partial

from eta_c_pseudoscalar.__main__ import _optional_anomaly
from inverse_common.cli import add_common_arguments, run_and_print
from .solver import DEFAULT_ANOMALY_1S, DS_DEFAULT, solve


def main():
    parser = argparse.ArgumentParser(
        description="eta_c P-P + one running-width BW: 5S primary; 4S alternative; "
                    "3S diagnostic (nearby poles); 6S unavailable (missing partial widths)")
    add_common_arguments(parser, channels=("3S", "4S", "5S", "6S"), ds_default=DS_DEFAULT)
    parser.add_argument("--alpha", type=float, nargs="+", default=None,
                        help="positive scalar or increasing values; default PP alpha grid, no selection")
    parser.add_argument("--radius", type=float, default=4.0, help="Blatt-Weisskopf radius in GeV^-1")
    parser.add_argument("--dstar-d1-l", type=int, choices=(0, 2), default=0)
    parser.add_argument("--quadrature-nodes", type=int, default=1024)
    parser.add_argument("--anomaly-1s", type=_optional_anomaly, default=DEFAULT_ANOMALY_1S)
    parser.add_argument("--anomaly-2s", type=_optional_anomaly, default=None)
    args = parser.parse_args()
    configured = partial(solve, alpha=args.alpha, radius=args.radius,
                         dstar_d1_l=args.dstar_d1_l, quadrature_nodes=args.quadrature_nodes,
                         anomaly_1s=args.anomaly_1s, anomaly_2s=args.anomaly_2s)
    run_and_print(parser, configured, args, prefix="eta_c_pp_bw")


if __name__ == "__main__":
    main()
