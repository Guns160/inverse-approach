import importlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

try:
    from scipy.io import loadmat
except ImportError:
    loadmat = None


ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ("Jpsi_vector", "eta_c_axial_vector", "eta_c_pseudoscalar")

def matlab_literal(value):
    return "'" + str(value).replace("'", "''") + "'"


class MatlabParityTest(unittest.TestCase):
    def test_current_python_and_matlab_three_and_four_loop(self):
        if os.environ.get("RUN_MATLAB_TESTS") != "1":
            self.skipTest("set RUN_MATLAB_TESTS=1 to launch MATLAB")
        if loadmat is None:
            self.skipTest("scipy.io is unavailable")
        matlab = os.environ.get("MATLAB_EXECUTABLE") or shutil.which("matlab")
        if matlab is None:
            matlab = "/Applications/MATLAB_R2026a.app/bin/matlab"
        matlab = Path(matlab)
        if not matlab.exists():
            self.skipTest("MATLAB executable is unavailable")

        sys.path.insert(0, str(ROOT))
        selected_loops = int(os.environ.get("MATLAB_PARITY_LOOPS", "4"))
        expected = {
            name: importlib.import_module(name).solve(
                5, 20, 20, ds=0.1, rg_loops=selected_loops
            )
            for name in PACKAGES
        }
        paths = [
            ROOT / "matlab",
            ROOT / "Jpsi_vector" / "matlab",
            ROOT / "eta_c_axial_vector" / "matlab",
            ROOT / "eta_c_pseudoscalar" / "matlab",
        ]
        with tempfile.TemporaryDirectory(prefix="charmonium-matlab-") as temporary:
            saved = Path(temporary) / "parity.mat"
            command = "".join(f"addpath({matlab_literal(path)});" for path in paths)
            command += "test_pp_anomaly; "
            command += (
                f"j=solve_jpsi_vector(5,20,20,'DS',0.1,'Mu',2,'RGLoops',{selected_loops},'CheckSolution',true); "
                f"a=solve_eta_c_axial_vector(5,20,20,'DS',0.1,'Mu',2,'RGLoops',{selected_loops},'CheckSolution',true); "
                f"p=solve_eta_c_pseudoscalar(5,20,20,'DS',0.1,'Mu',2,'RGLoops',{selected_loops},'CheckSolution',true); "
                f"save('-v7',{matlab_literal(saved)},'j','a','p');"
            )
            subprocess.run([str(matlab), "-batch", command], check=True)
            actual = loadmat(saved, squeeze_me=True, struct_as_record=False)

        saved_models = {
            "Jpsi_vector": "j",
            "eta_c_axial_vector": "a",
            "eta_c_pseudoscalar": "p",
        }
        for package in PACKAGES:
            want = expected[package]
            got = actual[saved_models[package]]
            with self.subTest(loops=selected_loops, package=package):
                self.assertEqual(float(got.metadata.M2_min), want["metadata"]["M2_min"])
                self.assertEqual(float(got.metadata.rg.loops), selected_loops)
                np.testing.assert_allclose(
                    got.M2, want["M2"], rtol=2e-15, atol=2e-15
                )
                for component, reference_ope in want["OPE"].items():
                    np.testing.assert_allclose(
                        getattr(got.OPE, component),
                        reference_ope,
                        rtol=2e-12,
                        atol=2e-12,
                    )
                np.testing.assert_allclose(
                    [float(got.metadata.mc), float(got.metadata.alpha_s)],
                    [want["metadata"]["mc"], want["metadata"]["alpha_s"]],
                    rtol=2e-11,
                    atol=2e-12,
                )
                self.assertEqual(list(want["channels"]), ["2pole"])
                self.assertEqual(
                    set(got._fieldnames), {"two_pole", "metadata", "M2", "weights", "OPE"}
                )
                self.assertEqual(str(got.metadata.channels), "2pole")
                reference = want["channels"]["2pole"]
                observed = got.two_pole
                self.assertLess(
                    observed.first_order.relative_direct_solution_error_inf, 2e-10
                )
                self.assertLess(
                    observed.first_order.relative_direct_gradient_inf, 2e-10
                )
                alpha_log_error = abs(
                    np.log10(float(observed.alpha_star) / float(reference["alpha_star"]))
                )
                self.assertLessEqual(alpha_log_error, 1e-8)
                np.testing.assert_array_equal(observed.masses, reference["masses"])
                np.testing.assert_allclose(
                    observed.residues, reference["residues"], rtol=3e-5, atol=3e-5
                )
                self.assertAlmostEqual(
                    float(observed.weighted_residual),
                    float(reference["weighted_residual"]),
                    delta=3e-5,
                )
                np.testing.assert_allclose(
                    observed.rho, reference["rho"], rtol=3e-5, atol=3e-5
                )
                np.testing.assert_array_equal(
                    np.isnan(np.asarray(observed.f)), np.isnan(reference["f"])
                )
                np.testing.assert_allclose(
                    observed.f, reference["f"], rtol=3e-5, atol=3e-5, equal_nan=True
                )
                if package == "eta_c_pseudoscalar":
                    np.testing.assert_allclose(observed.f_raw, reference["f_raw"],
                                               rtol=3e-5, atol=3e-5, equal_nan=True)
                    np.testing.assert_allclose(observed.anomaly_shift_GeV,
                                               reference["anomaly_shift_GeV"], rtol=1e-14)
                    np.testing.assert_array_equal(got.metadata.anomaly.applied,
                                                  want["metadata"]["anomaly"]["applied"])
                np.testing.assert_allclose(
                    observed.residual, reference["residual"], rtol=3e-5, atol=3e-5
                )
                residue_error = float(np.max(np.abs(observed.residues - reference["residues"])))
                print(
                    f"loops={selected_loops} {package}/2pole: "
                    f"abs(log10(alpha ratio))={alpha_log_error:.3e}, "
                    f"max_abs_residue_error={residue_error:.3e}"
                )


if __name__ == "__main__":
    unittest.main()
