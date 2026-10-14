"""Independent invariants of the physical inverse problem and public result format."""
from pathlib import Path
import importlib
import json
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
PACKAGES = ("Jpsi_vector", "eta_c_axial_vector", "eta_c_pseudoscalar")


class InverseProblemContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.outputs = {name: importlib.import_module(name).solve(
            5.003, 12.007, 20.123, ds=0.2) for name in PACKAGES}

    def test_grid_quadrature_and_boundary(self):
        for name, out in self.outputs.items():
            with self.subTest(package=name):
                self.assertEqual(out["M2"][0], 5.003)
                self.assertEqual(out["M2"][-1], 12.007)
                self.assertAlmostEqual(float(out["weights"].sum()), 12.007 - 5.003, places=12)
                self.assertTrue(np.all(out["weights"] > 0))
                for ch in out["channels"].values():
                    self.assertEqual(ch["s"][0], out["metadata"]["S_H"])
                    self.assertEqual(ch["s"][-1], 20.123)
                    self.assertTrue(np.all(np.diff(ch["s"]) > 0))
                    np.testing.assert_array_equal(ch["rho"][[0, -1]], ch["boundary"])
                    self.assertEqual(ch["boundary"][0], 0)

    def test_residue_units_and_pole_identity(self):
        for name, out in self.outputs.items():
            for ch in out["channels"].values():
                positive = ch["residues"] > 0
                factor = (2 * out["metadata"]["mc"] / ch["masses"] ** 2
                          if name == "eta_c_pseudoscalar" else np.ones_like(ch["masses"]))
                np.testing.assert_allclose(ch["f"][positive],
                                           factor[positive] * np.sqrt(ch["residues"][positive]),
                                           rtol=2e-14)
                self.assertTrue(np.isnan(ch["f"][~positive]).all())
        aa = self.outputs["eta_c_axial_vector"]["channels"]["2pole"]
        pp = self.outputs["eta_c_pseudoscalar"]["channels"]["2pole"]
        self.assertEqual(aa["masses"][1], 3.51067)
        self.assertEqual(pp["masses"][1], 3.6378)
        for out in self.outputs.values():
            self.assertEqual(set(out["channels"]), {"2pole"})
            self.assertEqual(len(out["channels"]["2pole"]["masses"]), 2)

    def test_back_substitution_and_first_order_condition(self):
        # Reconstruct the linear system independently from returned physical grids.
        from inverse_common.numerics import hats, h_mass_stiffness
        for name, out in self.outputs.items():
            t, w, g = out["M2"], out["weights"], out["OPE"]["total"]
            for channel, ch in out["channels"].items():
                with self.subTest(package=name, channel=channel):
                    pole_kernel = np.exp(-ch["masses"][None, :] ** 2 / t[:, None]) / t[:, None]
                    hat_kernel = hats(t, ch["s"])
                    predicted = np.dot(pole_kernel, ch["residues"]) + np.dot(hat_kernel, ch["rho"])
                    np.testing.assert_allclose(predicted, ch["backsub"], rtol=2e-13, atol=1e-13)
                    residual = predicted - g
                    np.testing.assert_allclose(residual, ch["residual"], rtol=1e-10, atol=1e-13)
                    self.assertAlmostEqual(float(np.sqrt(np.dot(w, residual ** 2))),
                                           ch["weighted_residual"], places=11)
                    pole_gradient = np.dot(pole_kernel.T, w * residual)
                    gradient = (np.dot(hat_kernel.T, w * residual)
                                + ch["alpha_star"] * np.dot(h_mass_stiffness(ch["s"]), ch["rho"]))
                    scale = max(1.0, float(np.linalg.norm(np.dot(hat_kernel.T, w * g))))
                    self.assertLess(float(np.max(np.abs(pole_gradient))) / scale, 2e-10)
                    self.assertLess(float(np.max(np.abs(gradient[1:-1]))) / scale, 2e-10)

    def test_json_and_npz_output(self):
        from inverse_common.cli import save_result
        for name, out in self.outputs.items():
            with tempfile.TemporaryDirectory() as temporary:
                saved = save_result(out, temporary, prefix=name)
                root = Path(saved["result_directory"])
                raw = (root / "summary.json").read_text()
                self.assertNotIn("NaN", raw)
                self.assertNotIn("Infinity", raw)
                self.assertEqual(json.loads(raw), saved)
                with np.load(root / "result.npz", allow_pickle=False) as arrays:
                    for key in arrays.files:
                        self.assertNotEqual(arrays[key].dtype, object)


class GeneralizedFilterContract(unittest.TestCase):
    def test_nontrivial_gsvd_mu_and_affine_boundary(self):
        """A well-conditioned independent least-squares reference catches missing mu²."""
        from inverse_common.factors import factors_one_window
        from inverse_common.numerics import trapw, h_mass_stiffness
        from inverse_common.reconstruction import reconstruct_alpha_family
        rng = np.random.default_rng(2718)
        s = np.array([15., 15.7, 16.9, 18.])
        t = np.linspace(2., 5., 17)
        K = rng.normal(size=(len(t), 6))
        g = rng.normal(size=len(t))
        H = h_mass_stiffness(s)
        boundary = np.array([0., 0.37])
        factors = factors_one_window(K, g, t, H, 2, boundary, "lapacke")
        self.assertGreater(float(np.max(np.abs(factors["mu"] - 1))), 0.05)
        alpha = np.array([0.01, 0.3, 2.])
        factors["alpha"] = alpha
        residues, rho, _ = reconstruct_alpha_family(
            factors, boundary, n_poles=2, s_grid=s, block_size=2)
        sw = np.sqrt(trapw(t))
        free, fixed = [0, 1, 3, 4], [2, 5]
        Hii, Hif = H[1:-1, 1:-1], H[1:-1][:, [0, -1]]
        q = np.linalg.solve(Hii, np.dot(Hif, boundary))
        R = np.linalg.cholesky(Hii).T
        L = np.column_stack((np.zeros((2, 2)), R))
        A = sw[:, None] * K[:, free]
        b = sw * (g - np.dot(K[:, fixed], boundary)) + np.dot(A[:, 2:], q)
        for j, a in enumerate(alpha):
            reference = np.linalg.lstsq(
                np.vstack((A, np.sqrt(a) * L)), np.r_[b, np.zeros(2)], rcond=None)[0]
            actual = np.r_[residues[:, j], rho[1:-1, j] + q]
            np.testing.assert_allclose(actual, reference, rtol=3e-12, atol=3e-12)


if __name__ == "__main__":
    unittest.main()
