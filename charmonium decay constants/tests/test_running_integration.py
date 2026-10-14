"""Running-integration tests for RG inputs and the corrected P-P GSVD form."""
from __future__ import annotations

import argparse
import inspect
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from eta_c_axial_vector import solve as solve_aa
from eta_c_pseudoscalar import lapacke_info, solve as solve_pp
from eta_c_pseudoscalar import physics_mu2 as pp_physics
from eta_c_pseudoscalar import raus2024_pp
from inverse_common.cli import add_common_arguments
from inverse_common.gsvd import compact_gsvd
from inverse_common.rg import RGInputs, derive
from inverse_common.reconstruction import reconstruct_alpha_family
from inverse_common.numerics import _dot, trapw
from Jpsi_vector import solve as solve_jpsi

WINDOW = (5.0, 20.0, 20.0)
CUSTOM_RG = RGInputs(mbar_c=1.28)


class RGIntegrationTests(unittest.TestCase):
    def test_all_channel_apis_have_running_arguments(self):
        for solve in (solve_jpsi, solve_aa, solve_pp):
            parameters = inspect.signature(solve).parameters
            for name in ("ds", "mu", "rg_loops", "rg_inputs"):
                self.assertIn(name, parameters)
        self.assertEqual(inspect.signature(solve_pp).parameters["rg_loops"].default, 4)

    def test_three_channels_use_the_same_rg_state_for_custom_anchor(self):
        results = [
            solve_jpsi(*WINDOW, channel="2pole", rg_inputs=CUSTOM_RG),
            solve_aa(*WINDOW, channel="2pole", rg_inputs=CUSTOM_RG),
            solve_pp(*WINDOW, channel="2pole", rg_inputs=CUSTOM_RG),
        ]
        for result in results:
            self.assertEqual(set(result["channels"]), {"2pole"})
        states = [result["metadata"]["rg"] for result in results]
        self.assertEqual(states[0], states[1])
        self.assertEqual(states[0], states[2])
        self.assertEqual(states[0]["inputs"]["mbar_c"], 1.28)
        self.assertEqual(states[0]["mass_anchor"]["mu"], 1.28)
        self.assertEqual(states[0]["mass_anchor"]["mc"], 1.28)
        self.assertTrue(states[0]["mc"] != derive()["mc"])

    def test_custom_anchor_does_not_leak_into_later_default_call(self):
        before = derive(mu=2.0, loops=4)
        custom = solve_pp(*WINDOW, channel="2pole", rg_inputs=CUSTOM_RG)
        after = solve_pp(*WINDOW, channel="2pole")
        self.assertEqual(custom["metadata"]["rg"]["inputs"]["mbar_c"], 1.28)
        self.assertEqual(after["metadata"]["rg"], before)

    def test_mu_and_loop_order_change_physical_pp_result(self):
        base = solve_pp(*WINDOW, channel="2pole")
        other_scale = solve_pp(*WINDOW, channel="2pole", mu=2.1)
        one_loop = solve_pp(*WINDOW, channel="2pole", rg_loops=1)
        self.assertEqual(other_scale["metadata"]["mu"], 2.1)
        self.assertNotEqual(base["metadata"]["mc"], other_scale["metadata"]["mc"])
        self.assertNotEqual(base["metadata"]["alpha_s"], other_scale["metadata"]["alpha_s"])
        self.assertFalse(np.array_equal(base["OPE"]["total"], other_scale["OPE"]["total"]))
        self.assertEqual(one_loop["metadata"]["rg"]["loops"], 1)
        self.assertNotEqual(base["metadata"]["mc"], one_loop["metadata"]["mc"])
        self.assertFalse(np.array_equal(base["OPE"]["total"], one_loop["OPE"]["total"]))

    def test_every_pp_physical_component_is_bound_to_running_inputs(self):
        m2 = np.array([5.0, 10.0, 20.0])
        base = pp_physics.pp_ope(20.0, m2, mc=1.0, alpha_s=0.3, components=True)
        alpha_changed = pp_physics.pp_ope(20.0, m2, mc=1.0, alpha_s=0.31, components=True)
        mass_changed = pp_physics.pp_ope(20.0, m2, mc=1.01, alpha_s=0.3, components=True)
        np.testing.assert_array_equal(base[0], alpha_changed[0])
        self.assertFalse(np.array_equal(base[1], alpha_changed[1]))
        np.testing.assert_array_equal(base[2], alpha_changed[2])
        np.testing.assert_array_equal(base[3], alpha_changed[3])
        for index in range(4):
            self.assertFalse(np.array_equal(base[index], mass_changed[index]))
        boundary_base = pp_physics.pp_rho_pert(np.array([20.0]), mc=1.0, alpha_s=0.3)
        boundary_changed = pp_physics.pp_rho_pert(np.array([20.0]), mc=1.01, alpha_s=0.31)
        self.assertNotEqual(float(boundary_base[0]), float(boundary_changed[0]))

    def test_pp_converter_uses_the_same_derived_mass(self):
        result = solve_pp(*WINDOW, channel="2pole", rg_inputs=CUSTOM_RG)
        channel = result["channels"]["2pole"]
        mc = result["metadata"]["mc"]
        positive = channel["residues"] > 0
        expected = 2.0 * mc * np.sqrt(channel["residues"][positive]) / channel["masses"][positive] ** 2
        np.testing.assert_allclose(channel["f"][positive], expected, rtol=1e-14)

    def test_pp_selector_receives_generalized_singular_values(self):
        result = solve_pp(*WINDOW, channel="2pole")
        channel = result["channels"]["2pole"]
        factors = channel["factors"]
        self.assertEqual(
            channel["alpha_diagnostics"]["sigma_min_sq"],
            float(np.min(factors["sigma_tilde"]) ** 2),
        )
        selected = raus2024_pp.select_row(
            factors["sigma_tilde"], factors["beta_c"],
            np.array([0, len(factors["sigma_tilde"])], dtype=np.int64),
            np.array([factors["null"]], float), factors["alpha"],
            channel["qoc"][None, :],
        )
        self.assertEqual(int(selected["index"][0]), channel["selection_index"])

    def test_common_cli_exposes_running_options(self):
        parser = argparse.ArgumentParser()
        add_common_arguments(parser, channels=("2pole",), ds_default=0.1)
        args = parser.parse_args(["5", "20", "20", "--ds", "0.1", "--mu", "2.1", "--rg-loops", "3"])
        self.assertEqual(args.mu, 2.1)
        self.assertEqual(args.rg_loops, 3)


@unittest.skipUnless(lapacke_info()["available"], "LP64 LAPACKE_dggsvd3 unavailable")
class GeneralizedGSVDTests(unittest.TestCase):
    def test_non_near_one_gsvd_stationarity_and_qoc(self):
        rng = np.random.default_rng(20261014)
        m, n, n_poles = 14, 10, 3
        n_internal = n - n_poles
        A = rng.normal(size=(m, n))
        q, _ = np.linalg.qr(rng.normal(size=(n_internal, n_internal)))
        scales = np.diag(np.logspace(-0.7, 0.8, n_internal))
        L = np.column_stack((np.zeros((n_internal, n_poles)), q @ scales))
        b = rng.normal(size=m)
        decomp = compact_gsvd(A, L, method="lapacke", qr_reduce=False)
        self.assertGreater(float(np.max(np.abs(decomp.mu - 1.0))), 0.1)
        self.assertLess(float(np.min(decomp.mu)), 0.2)
        beta = _dot(decomp.U.T, b)
        factors = {
            "alpha": np.array([1e-4, 1e-3, 1e-2, 1e-1]),
            "sigma": decomp.sigma,
            "mu": decomp.mu,
            "sigma_unreg": decomp.sigma_unregularized,
            "beta_c": beta[:n_internal],
            "beta_n": beta[n_internal:],
            "null": 0.0,
            "X": decomp.X,
            "q": np.zeros(n - n_poles),
        }
        residues, rho, qoc = reconstruct_alpha_family(
            factors, np.zeros(2), n_poles=n_poles,
            s_grid=np.linspace(1.0, 2.0, n_internal + 2), block_size=3)
        for index, alpha in enumerate(factors["alpha"]):
            x = np.vstack((
                residues[:, index, np.newaxis],
                rho[1:-1, index, np.newaxis],
            )).reshape(-1)
            stationarity = _dot(A.T, _dot(A, x) - b) + alpha * _dot(L.T, _dot(L, x))
            scale = np.linalg.norm(_dot(A.T, b)) + alpha * np.linalg.norm(_dot(L.T, L)) * np.linalg.norm(x)
            self.assertLess(np.linalg.norm(stationarity) / scale, 2e-12)

        center = 1e-3
        h = 1e-3
        probe_factors = dict(factors)
        probe_factors["alpha"] = center * 10.0 ** np.array([-h, 0.0, h])
        _, probe_rho, probe_qoc = reconstruct_alpha_family(
            probe_factors, np.zeros(2), n_poles=n_poles,
            s_grid=np.linspace(1.0, 2.0, n_internal + 2), block_size=3)
        derivative = (probe_rho[:, 2] - probe_rho[:, 0]) / (2.0 * h)
        finite_difference = np.sum(
            trapw(np.linspace(1.0, 2.0, n_internal + 2))[1:-1]
            * derivative[1:-1] ** 2
        )
        self.assertAlmostEqual(probe_qoc[1], finite_difference,
                               delta=2e-5 * finite_difference)


if __name__ == "__main__":
    unittest.main()
