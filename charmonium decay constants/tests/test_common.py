"""Regression tests for the shared inversion kernels and channel APIs."""
from __future__ import annotations

import inspect
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from eta_c_pseudoscalar import lapacke_info, solve as solve_pp
from inverse_common.numerics import borel_grid, sgrid_hats, trapw
from inverse_common.reconstruction import reconstruct_alpha_family
from inverse_common.validation import validate_window
from Jpsi_vector import solve as solve_jpsi
from eta_c_axial_vector import solve as solve_aa


class GeometryTests(unittest.TestCase):
    def test_hat_grid_retains_exact_lambda_and_short_cell(self):
        grid = sgrid_hats(20.037, step=0.1, s_h=14.99)
        self.assertEqual(grid[0], 14.99)
        self.assertEqual(grid[-1], 20.037)
        self.assertEqual(len(grid), 52)
        self.assertLess(grid[-1] - grid[-2], 0.1)
        self.assertTrue(np.all(np.diff(grid) > 0))

    def test_borel_grid_retains_requested_endpoints(self):
        grid = borel_grid(5.003, 20.007, step=0.01)
        self.assertEqual(grid[0], 5.003)
        self.assertEqual(grid[-1], 20.007)
        self.assertTrue(np.all(np.diff(grid) > 0))

    def test_validation_rejects_nonpositive_ds(self):
        with self.assertRaisesRegex(ValueError, "ds"):
            validate_window(5, 20, 20, masses={"m": (3.0,)}, channel="m",
                            s_h=14.99, ds=0)


class ChannelApiTests(unittest.TestCase):
    def test_all_solve_apis_accept_ds(self):
        for solve in (solve_jpsi, solve_aa, solve_pp):
            self.assertIn("ds", inspect.signature(solve).parameters)

    def test_default_and_explicit_two_pole_models(self):
        cases = (
            (solve_jpsi, (3.096900, 3.686097)),
            (solve_aa, (2.9841, 3.51067)),
            (solve_pp, (2.9841, 3.6378)),
        )
        for solve, masses in cases:
            with self.subTest(solver=solve.__module__):
                default = solve(5.0, 20.0, 20.0)
                explicit = solve(5.0, 20.0, 20.0, channel="2pole")
                for result in (default, explicit):
                    self.assertEqual(set(result["channels"]), {"2pole"})
                    model = result["channels"]["2pole"]
                    np.testing.assert_array_equal(model["masses"], masses)
                    self.assertEqual(model["residues"].shape, (2,))
                    self.assertEqual(model["f"].shape, (2,))
                for key in ("residues", "rho", "alpha_star"):
                    np.testing.assert_array_equal(
                        default["channels"]["2pole"][key],
                        explicit["channels"]["2pole"][key],
                    )

    def test_other_counts_and_model_names_are_rejected(self):
        unsupported = [
            f"{count}{name}"
            for count in (1, 2, 3) for name in ("delta", "pole")
            if (count, name) != (2, "pole")
        ]
        for solve in (solve_jpsi, solve_aa, solve_pp):
            for channel in unsupported:
                with self.subTest(solver=solve.__module__, channel=channel):
                    with self.assertRaisesRegex(ValueError, "channel"):
                        solve(5.0, 20.0, 20.0, channel=channel)


@unittest.skipUnless(lapacke_info()["available"], "LP64 LAPACKE_dggsvd3 unavailable")
class PPReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = solve_pp(5.0, 20.0, 20.0, channel="2pole")
        cls.channel = cls.result["channels"]["2pole"]

    def test_blocked_full_alpha_family_is_independent_of_block_size(self):
        factors = dict(self.channel["factors"])
        boundary = self.channel["boundary"]
        fast = reconstruct_alpha_family(
            factors, boundary, n_poles=2, s_grid=self.channel["s"], block_size=128)
        slow = reconstruct_alpha_family(
            factors, boundary, n_poles=2, s_grid=self.channel["s"], block_size=17)
        for a, b in zip(fast, slow):
            np.testing.assert_array_equal(a, b)

    def test_qoc_matches_central_finite_difference(self):
        factors = dict(self.channel["factors"])
        alpha_star = self.channel["alpha_star"]
        h_log10 = 1e-4
        factors["alpha"] = alpha_star * 10.0 ** np.array([-h_log10, 0.0, h_log10])
        _, rho, qoc = reconstruct_alpha_family(
            factors, self.channel["boundary"], n_poles=2,
            s_grid=self.channel["s"], block_size=3)
        derivative = (rho[:, 2] - rho[:, 0]) / (2.0 * h_log10)
        weights = trapw(self.channel["s"])[1:-1]
        finite_difference = np.sum(weights * derivative[1:-1] ** 2)
        self.assertAlmostEqual(qoc[1], finite_difference, delta=2e-5 * finite_difference)


if __name__ == "__main__":
    unittest.main()
