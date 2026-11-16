import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad
from scipy.linalg import cholesky, lstsq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
S = importlib.import_module("eta_c_pseudoscalar+BW.solver")
P = importlib.import_module("eta_c_pseudoscalar+BW.physics")
F = importlib.import_module("eta_c_pseudoscalar+BW.factors")
from eta_c_pseudoscalar import physics_mu2 as PP
from inverse_common.cli import save_result
from inverse_common.numerics import hats, h_mass_stiffness


class RunningWidthTest(unittest.TestCase):
    def test_state_inputs_and_pole_sums(self):
        for state, mass, width, count in (("3S", 4.0119, .07374374, 1),
                                           ("4S", 4.2693, .03356257, 4),
                                           ("5S", 4.4508, .07952178, 8)):
            with self.subTest(state=state):
                model = P.RunningWidth(state)
                self.assertEqual(model.mass, mass)
                self.assertEqual(len(model.pole_partial_widths), count)
                self.assertAlmostEqual(float(model.gamma(mass*mass)), width, places=14)
                self.assertAlmostEqual(model.pole_partial_widths.sum(), model.inputs["open_width"], places=14)
                np.testing.assert_allclose(model.partial_widths(mass*mass), model.pole_partial_widths, rtol=2e-15)
                self.assertAlmostEqual(float(model.bw(mass*mass)), 1/(np.pi*mass*width), places=13)
                np.testing.assert_array_equal(model.partial_widths([1.0]), np.zeros((1, count)))
                self.assertEqual(float(model.gamma(1.0)), model.inputs["rest_width"])
                self.assertGreater(float(model.bw(1.0)), 0)

    def test_reference_reweighting(self):
        for state in ("4S", "5S"):
            model = P.RunningWidth(state)
            ref = model.phase_space(model.inputs["reference_mass"]**2)
            ratios = model.pole_partial_widths / (np.asarray(model.inputs["partial_widths"]) *
                                                    model.phi_pole / ref)
            np.testing.assert_allclose(ratios, np.repeat(ratios[0], len(ratios)), rtol=2e-15)

    def test_charge_average_and_threshold(self):
        model = P.RunningWidth("3S")
        s = np.array([14.0, (sum((P.D[0], P.DSTAR[0])))**2, 15.2, model.mass**2, 30.0])
        expected = np.zeros_like(s)
        for a, b in zip(P.D, P.DSTAR):
            q = np.sqrt(np.maximum((s-(a+b)**2)*(s-(a-b)**2), 0))/(2*np.sqrt(s))
            expected += .5*np.where(s > (a+b)**2, q**3/np.sqrt(s)/(1+(q*4)**2), 0)
        np.testing.assert_allclose(model.phase_space(s)[:, 0], expected, rtol=2e-15, atol=1e-18)

    def test_broad_convolution_against_adaptive_quad(self):
        model = P.RunningWidth("5S")
        s = 19.80962064
        for a, b, width, daughter, ell in ((P.D[0], P.D0STAR[0], P.D0_WIDTH, P.D[1], 0),
                                           (P.DSTAR[0], P.D1[0], P.D1_WIDTH, P.DSTAR[1], 2)):
            lo = (daughter+P.PIC)**2
            hi = (np.sqrt(s)-a)**2
            def integrand(u):
                q = np.sqrt(max((s-(a+np.sqrt(u))**2)*(s-(a-np.sqrt(u))**2), 0))/(2*np.sqrt(s))
                z2 = (4*q)**2
                barrier = 1 if ell == 0 else 1/(9+3*z2+z2*z2)
                return b*width/(np.pi*((u-b*b)**2+(b*width)**2))*q**(2*ell+1)/np.sqrt(s)*barrier
            expected = quad(integrand, lo, hi, epsabs=1e-14, epsrel=1e-10,
                            points=[b*b] if lo < b*b < hi else None, limit=300)[0]
            actual = model._convolution(np.array([s]), a, b, width, lo, ell)[0]
            self.assertAlmostEqual(actual/expected, 1.0, delta=2e-8)

    def test_radius_zero_and_partial_wave_options(self):
        s = np.array([15.0, 19.80962064, 34.25])
        for ell in (0, 2):
            zero = P.RunningWidth("5S", radius=0, dstar_d1_l=ell)
            tiny = P.RunningWidth("5S", radius=1e-7, dstar_d1_l=ell)
            np.testing.assert_allclose(zero.gamma(s), tiny.gamma(s), rtol=1e-12)
            self.assertAlmostEqual(float(zero.gamma(zero.mass**2)), .07952178, places=14)
        a = P.RunningWidth("5S", radius=0)
        b = P.RunningWidth("5S", radius=4)
        np.testing.assert_array_equal(a.phase_space(s)[:, [2, 3, 7]], b.phase_space(s)[:, [2, 3, 7]])

    def test_shapes_and_quadrature_convergence(self):
        model = P.RunningWidth("5S")
        s = np.array([[16.0, 19.0], [21.0, 30.0]])
        self.assertEqual(model.bw(s).shape, s.shape)
        self.assertEqual(model.partial_widths(s).shape, (2, 2, 8))
        fine = P.RunningWidth("5S", quadrature_nodes=2048)
        np.testing.assert_allclose(model.gamma(s), fine.gamma(s), rtol=2e-8)

    def test_independent_borel_integral_and_truncated_peak(self):
        model = P.RunningWidth("4S")
        t = np.array([5.0, 10.0, 20.0])
        for upper in (17.0, 25.0):
            column, info = model.borel_column(t, 14.99, upper)
            expected = [quad(lambda x: float(model.bw(x))*np.exp(-x/v)/v,
                             14.99, upper, epsabs=1e-13, epsrel=1e-10,
                             points=[x for x in model.thresholds()+[model.mass**2] if 14.99<x<upper],
                             limit=300)[0] for v in t]
            np.testing.assert_allclose(column, expected, rtol=2e-10, atol=1e-13)
            self.assertFalse(info["area_renormalized"])

    def test_invalid_inputs_and_unavailable_6s(self):
        for kwargs in ({"state": "6S"}, {"state": "7S"}, {"radius": -1}, {"radius": np.nan},
                       {"radius": True}, {"radius": np.bool_(True)}, {"radius": np.array(True)},
                       {"dstar_d1_l": 1}, {"dstar_d1_l": True},
                       {"quadrature_nodes": 0}, {"quadrature_nodes": 16.5}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                P.RunningWidth(**kwargs)
        with patch.object(S, "lapacke_info", side_effect=AssertionError("must reject before backend")):
            with self.assertRaisesRegex(ValueError, "6S.*partial widths"):
                S.solve(5, 20, 20, channel="6S", alpha=1e-6)


class AffineBWTest(unittest.TestCase):
    def test_synthetic_augmented_least_squares(self):
        rng = np.random.default_rng(20261116)
        n, m = 7, 25
        poles, G, k, g = rng.normal(size=(m, 2)), rng.normal(size=(m, n)), rng.normal(size=m), rng.normal(size=m)
        weights = rng.uniform(.1, 1, m)
        V = rng.normal(size=(n, n))
        H = V.T@V + np.eye(n)
        endpoints, boundary = np.array([.02, .03]), np.array([0., 1.2])
        system = F.assemble(poles, G, k, g, weights, H, endpoints, boundary)
        self.assertEqual(system["B"].shape, (n-1, n+1))
        np.testing.assert_array_equal(system["B"][:, :2], 0)
        self.assertGreater(np.linalg.norm(system["B"][:, 2]), 0)
        fac = F.factors_one_window(system, g, weights)
        alpha = np.array([1e-5, .2, 10.0])
        R, c = F.reconstruct_alpha_family(fac, system, alpha)
        E = np.zeros((n, n-1)); E[[0, -1], 0] = -endpoints; E[1:-1, 1:] = np.eye(n-2)
        a = np.zeros(n); a[-1] = boundary[-1]
        K = np.column_stack((poles, k+G@E[:, 0], G[:, 1:-1]))
        LH = cholesky(H, lower=False)
        B = np.column_stack((np.zeros((n, 2)), LH@E))
        for j, aa in enumerate(alpha):
            x = lstsq(np.vstack((np.sqrt(weights)[:, None]*K, np.sqrt(aa)*B)),
                      np.r_[np.sqrt(weights)*(g-G@a), -np.sqrt(aa)*(LH@a)])[0]
            np.testing.assert_allclose(R[:, j], x[:3], rtol=2e-11, atol=2e-12)
            np.testing.assert_allclose(c[:, j], a+E@x[2:], rtol=2e-11, atol=2e-12)
            shifted = x[2:]+system["q"]
            self.assertAlmostEqual(float(c[:, j]@H@c[:, j]),
                                   float(np.linalg.norm(system["T"]@shifted)**2+system["kappa"]), places=10)
        d = F.diagnose_family(poles, G, k, g, weights, H, system, R, c, alpha)
        self.assertLess(d["stationarity_relative_all"].max(), 1e-12)
        np.testing.assert_allclose(c[[0, -1]]+endpoints[:, None]*R[2],
                                   np.repeat(boundary[:, None], len(alpha), axis=1), atol=2e-16)


class PhysicalBWTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = {state: S.solve(5, 20, 20, channel=state, ds=.5, alpha=[1e-10, 1e-6, 1e-2])
                       for state in ("3S", "4S", "5S")}

    def test_original_ope_grid_and_signed_reconstruction(self):
        for state, result in self.results.items():
            c = result["channels"][state]
            meta = result["metadata"]
            components = PP.pp_ope(20, result["M2"], mc=meta["mc"], alpha_s=meta["alpha_s"], components=True)
            for key, expected in zip(("LO", "NLO", "D4", "D6"), components):
                np.testing.assert_array_equal(result["OPE"][key], expected)
            np.testing.assert_allclose(c["rho_all"]*c["s"][:, None], c["c_all"], rtol=2e-15)
            np.testing.assert_allclose(c["rho_total_all"]*c["s"][:, None], c["total_hat_all"], rtol=2e-15)
            self.assertLess(np.max(abs(c["endpoint_error_all"])), 1e-14)
            self.assertLess(c["stationarity_relative_all"].max(), 1e-10)
            self.assertTrue(np.all(c["validation_pass_all"]))
            self.assertTrue(np.any(c["residues_all"] < 0))
            self.assertTrue(np.any(c["c_all"] < 0))
            np.testing.assert_array_equal(np.isnan(c["f_raw_all"]), c["residues_all"] <= 0)
            self.assertNotIn("alpha_star", c)
            self.assertEqual(meta["alpha_selection"], "none")
            good = np.isfinite(c["f_all"])
            shifts = np.broadcast_to(c["anomaly_shift_GeV"][:, None], c["f_all"].shape)
            np.testing.assert_allclose((c["f_all"]-c["f_raw_all"])[good], shifts[good], atol=1e-16)

    def test_original_unshifted_objective(self):
        for state, result in self.results.items():
            c = result["channels"][state]
            s, t, w = c["s"], result["M2"], result["weights"]
            G, H = hats(t, s), h_mass_stiffness(s)
            poles = np.exp(-c["masses"][:2]**2/t[:, None])/t[:, None]
            prediction = np.einsum("ij,jk->ik", poles, c["residues_all"][:2])
            prediction += c["k_B"][:, None]*c["residues_all"][2] + np.einsum("ij,jk->ik", G, c["c_all"])
            np.testing.assert_allclose(prediction, c["backsub_all"], rtol=1e-12, atol=1e-12)
            rss = np.sum(w[:, None]*(prediction-result["OPE"]["total"][:, None])**2, axis=0)
            penalty = np.sum(c["c_all"]*(H@c["c_all"]), axis=0)
            np.testing.assert_allclose(c["objective_all"], rss+c["alpha"]*penalty, rtol=1e-9, atol=1e-12)

    def test_anomaly_does_not_change_inversion(self):
        a = self.results["4S"]["channels"]["4S"]
        other = S.solve(5, 20, 20, channel="4S", ds=.5, alpha=a["alpha"], anomaly_1s=None, anomaly_2s=0)
        b = other["channels"]["4S"]
        for key in ("residues_all", "c_all", "f_raw_all", "backsub_all", "objective_all"):
            np.testing.assert_array_equal(a[key], b[key])
        self.assertEqual(other["metadata"]["anomaly"]["applied"], [False, True, False])

    def test_truncated_bw_and_free_dimension(self):
        out = S.solve(5, 20, 18, channel="5S", alpha=1e-6, ds=.5)
        c = out["channels"]["5S"]
        self.assertLess(out["metadata"]["Lambda"], c["masses"][2]**2)
        self.assertGreater(c["k_B"].min(), 0)
        self.assertEqual(c["alpha_diagnostics"]["source"], "user_supplied")
        with self.assertRaisesRegex(ValueError, "Borel points"):
            S.solve(5, 5.5, 20, channel="3S", alpha=1e-6, ds=.1)

    def test_invalid_alpha_and_default_grid(self):
        np.testing.assert_array_equal(S._alpha_grid(None), PP.ALPHAS)
        for value in (True, [], [1e-4, 1e-5], [1e-4, 1e-4], 0, np.nan, np.inf, [1+0j], "1e-4"):
            with self.subTest(alpha=value), self.assertRaises(ValueError):
                S.solve(5, 20, 20, alpha=value)

    def test_full_default_alpha_scan_has_no_selected_point(self):
        out = S.solve(5, 8, 20, channel="3S", ds=.5)
        c = out["channels"]["3S"]
        np.testing.assert_array_equal(c["alpha"], PP.ALPHAS)
        self.assertEqual(c["residues_all"].shape, (3, 2001))
        self.assertNotIn("alpha_star", c)
        self.assertTrue(np.all(c["validation_pass_all"]))

    def test_accepted_running_width_reference_windows(self):
        # References are saved running-width solves; they do not select a platform.
        cases = (
            ("3S", 20.0, 20.0, -8.5, [2.856935614856, -.765277161715, 4.751280141085]),
            ("4S", 20.75, 24.0, -12.4, [2.7401842616951058, .5519016347236629, 4.547692224221237]),
        )
        for state, lo, lam, ell, expected in cases:
            with self.subTest(state=state):
                c = S.solve(lo, 80, lam, channel=state, alpha=10.**ell)["channels"][state]
                np.testing.assert_allclose(c["residues"], expected, rtol=5e-6, atol=1e-8)
                self.assertTrue(np.all(c["validation_pass_all"]))
        c = S.solve(19.375, 80, 34.25, alpha=10.**-13.85)["channels"]["5S"]
        np.testing.assert_allclose(c["f"]*1000, [399.9058667, 152.1886562, 300.6551655],
                                   rtol=5e-6, atol=1e-5)
        self.assertTrue(np.all(c["validation_pass_all"]))

    def test_pickle_free_family_persistence(self):
        with tempfile.TemporaryDirectory() as tmp:
            one = save_result(self.results["5S"], tmp, prefix="bw")
            two = save_result(self.results["5S"], tmp, prefix="bw")
            self.assertNotEqual(one["result_directory"], two["result_directory"])
            self.assertEqual(one["channels"]["5S"]["n_alpha"], 3)
            self.assertIsNone(one["channels"]["5S"]["alpha_selection"])
            path = Path(one["result_directory"])
            with np.load(path/"result.npz", allow_pickle=False) as z:
                for key in z.files:
                    self.assertNotEqual(z[key].dtype.kind, "O")
                np.testing.assert_array_equal(z["channels__5S__residues_all"],
                                              self.results["5S"]["channels"]["5S"]["residues_all"])
            json.loads((path/"summary.json").read_text())

    def test_cli_any_working_directory_and_6s_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            command = [sys.executable, "-B", str(ROOT/"eta_c_pseudoscalar+BW"/"run.py"),
                       "5", "20", "20", "--ds", "0.5", "--alpha", "1e-6", "--output", tmp]
            completed = subprocess.run(command, cwd=tmp, capture_output=True, text=True, check=True)
            summary = json.loads(completed.stdout)
            self.assertEqual(summary["channels"]["5S"]["source"], "user_supplied")
            count = len(list(Path(tmp).iterdir()))
            failed = subprocess.run(command+["--channel", "6S"], cwd=tmp, capture_output=True, text=True)
            self.assertEqual(failed.returncode, 2)
            self.assertIn("partial widths are missing", failed.stderr)
            self.assertEqual(len(list(Path(tmp).iterdir())), count)


if __name__ == "__main__":
    unittest.main()
