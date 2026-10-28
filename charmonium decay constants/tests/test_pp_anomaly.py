"""Axial Ward identity, unchanged inversion, and persisted PP conventions."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from eta_c_pseudoscalar import solve
from eta_c_pseudoscalar.solver import _pp_f
from inverse_common.cli import save_result


class PPAnomalyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        window = (5.0, 20.0, 20.0)
        cls.default = solve(*window, ds=0.2)
        cls.omitted = solve(*window, ds=0.2, anomaly_1s=None)
        cls.custom = solve(*window, ds=0.2, anomaly_1s=0.0, anomaly_2s=0.02)

    def test_default_shift_and_ward_identity(self):
        out = self.default
        ch = out["channels"]["2pole"]
        shift = np.array([-0.079 / 2.9841**2, 0.0])
        np.testing.assert_allclose(ch["f"] - ch["f_raw"], shift, atol=1e-15)
        np.testing.assert_allclose(ch["f_all"], ch["f_raw_all"] + shift[:, None],
                                   rtol=1e-15, equal_nan=True)
        mc = out["metadata"]["mc"]
        recovered_h = (ch["f"] * ch["masses"]**2 - np.array([-0.079, 0.0])) / (2 * mc)
        np.testing.assert_allclose(recovered_h**2, ch["residues"], rtol=2e-14)
        self.assertEqual(ch["f_eta_c"], ch["f"][0])
        self.assertEqual(out["metadata"]["anomaly"]["applied"], [True, False])
        self.assertIsNone(out["metadata"]["anomaly"]["eta_c_2s_GeV3"])

    def test_anomaly_does_not_change_inversion(self):
        base = self.default
        for other in (self.omitted, self.custom):
            for key in base["OPE"]:
                np.testing.assert_array_equal(base["OPE"][key], other["OPE"][key])
            ch = other["channels"]["2pole"]
            ref = base["channels"]["2pole"]
            for key in ("residues_all", "rho_all", "backsub", "qoc", "f_raw_all"):
                np.testing.assert_array_equal(ch[key], ref[key])
            self.assertEqual(ch["alpha_star"], ref["alpha_star"])
        omitted = self.omitted["channels"]["2pole"]
        np.testing.assert_array_equal(omitted["f"], omitted["f_raw"])

    def test_explicit_zero_differs_from_missing_input(self):
        metadata = self.custom["metadata"]["anomaly"]
        self.assertEqual(metadata["applied"], [True, True])
        self.assertEqual(metadata["eta_c_1s_GeV3"], 0.0)
        ch = self.custom["channels"]["2pole"]
        np.testing.assert_allclose(ch["f"] - ch["f_raw"],
                                   [0, 0.02 / 3.6378**2], atol=1e-15)

    def test_nonpositive_residues_stay_undefined(self):
        r = np.array([[-1.0, 0.0, 4.0], [0.0, -3.0, 1.0]])
        masses = np.array([2.9841, 3.6378])
        raw = _pp_f(r, masses, mc=1.1)
        corrected = raw + np.array([-0.079, 0.02])[:, None] / masses[:, None]**2
        self.assertTrue(np.isnan(corrected[r <= 0]).all())
        self.assertTrue(np.isfinite(corrected[r > 0]).all())

    def test_invalid_anomaly_is_rejected_before_rg_or_inversion(self):
        for value in (True, np.nan, np.inf, -np.inf, [0.0], 1j, "bad"):
            with self.subTest(value=value), patch("eta_c_pseudoscalar.solver.derive") as derive:
                with self.assertRaises(ValueError):
                    solve(5, 20, 20, anomaly_1s=value)
                derive.assert_not_called()

    def test_summary_and_npz_preserve_corrected_and_raw_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            summary = save_result(self.default, temporary, prefix="pp_anomaly")
            ch = summary["channels"]["2pole"]
            self.assertIn("f_raw_GeV", ch)
            self.assertIsNone(summary["metadata"]["anomaly"]["eta_c_2s_GeV3"])
            root = Path(summary["result_directory"])
            self.assertEqual(json.loads((root / "summary.json").read_text()), summary)
            with np.load(root / "result.npz", allow_pickle=False) as arrays:
                np.testing.assert_allclose(arrays["channels__2pole__f"], ch["f_GeV"])
                np.testing.assert_allclose(arrays["channels__2pole__f_raw"], ch["f_raw_GeV"])
                for key in arrays.files:
                    self.assertNotEqual(arrays[key].dtype, object)

    def test_cli_can_omit_ground_anomaly_and_supply_excited_anomaly(self):
        with tempfile.TemporaryDirectory() as temporary:
            proc = subprocess.run([
                sys.executable, "-B", "-m", "eta_c_pseudoscalar", "5", "8", "16",
                "--ds", "0.5", "--anomaly-1s", "none", "--anomaly-2s", "0.02",
                "--output", temporary,
            ], cwd=ROOT, capture_output=True, text=True, check=True)
            summary = json.loads(proc.stdout)
            self.assertEqual(summary["metadata"]["anomaly"]["applied"], [False, True])
            self.assertEqual(summary["metadata"]["anomaly"]["eta_c_2s_GeV3"], 0.02)
            self.assertIn("f_raw_GeV", summary["channels"]["2pole"])


if __name__ == "__main__":
    unittest.main()
