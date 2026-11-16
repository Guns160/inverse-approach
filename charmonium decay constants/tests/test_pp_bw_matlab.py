"""Opt-in Python/MATLAB parity for the independent two-delta + one-BW channel."""
import importlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]


def matlab_literal(value):
    return "'" + str(value).replace("'", "''") + "'"


class PPBWMatlabParityTest(unittest.TestCase):
    def test_running_widths_and_three_channel_families(self):
        if os.environ.get("RUN_MATLAB_TESTS") != "1":
            self.skipTest("set RUN_MATLAB_TESTS=1 to launch MATLAB")
        matlab = (os.environ.get("MATLAB_EXECUTABLE") or shutil.which("matlab")
                  or "/Applications/MATLAB_R2026a.app/bin/matlab")
        if not Path(matlab).exists():
            self.skipTest("MATLAB executable is unavailable")
        sys.path.insert(0, str(ROOT))
        package = importlib.import_module("eta_c_pseudoscalar+BW")
        physics = importlib.import_module("eta_c_pseudoscalar+BW.physics")
        alpha = np.array([1e-10, 1e-6, 1e-2])
        states = ("3S", "4S", "5S")
        expected = {
            state: package.solve(5, 8, 20, channel=state, ds=1, alpha=alpha)
            for state in states
        }
        truncated = package.solve(5, 8, 16, ds=0.5, alpha=1e-6,
                                  anomaly_1s=None, anomaly_2s=0.02)
        probe = np.array([1.0, 14.99, 15.4, 16.5, 20.0, 25.0, 35.0])
        models = [physics.RunningWidth(state) for state in states]
        models.append(physics.RunningWidth("5S", radius=0, dstar_d1_l=2))
        with tempfile.TemporaryDirectory(prefix="charmonium-pp-bw-matlab-") as temporary:
            saved = Path(temporary) / "bw-parity.mat"
            command = (
                f"addpath({matlab_literal(ROOT / 'matlab')});"
                f"addpath({matlab_literal(ROOT / 'eta_c_pseudoscalar+BW' / 'matlab')});"
                "test_pp_bw;"
                "states={'3S','4S','5S'}; results=cell(1,3);"
                "for j=1:3, results{j}=solve_eta_c_pseudoscalar_bw(5,8,20,"
                "'Channel',states{j},'DS',1,'Alpha',[1e-10,1e-6,1e-2]); end;"
                "truncated=solve_eta_c_pseudoscalar_bw(5,8,16,'DS',0.5,'Alpha',1e-6,"
                "'Anomaly1S',[],'Anomaly2S',0.02);"
                "probe=[1;14.99;15.4;16.5;20;25;35];"
                "width_shapes=cell(1,4); width_gammas=cell(1,4);"
                "width_partials=cell(1,4); width_poles=cell(1,4);"
                "for j=1:4, if j<=3, m=pp_bw_model(states{j});"
                "else, m=pp_bw_model('5S','Radius',0,'DstarD1L',2); end;"
                "[width_shapes{j},width_gammas{j},width_partials{j}]=pp_bw_spectrum(probe,m);"
                "width_poles{j}=m.pole_partial_widths; end;"
                f"save('-v7',{matlab_literal(saved)},'results','truncated',"
                "'width_shapes','width_gammas','width_partials','width_poles');"
            )
            subprocess.run([str(matlab), "-batch", command], check=True, timeout=600)
            actual = loadmat(saved, squeeze_me=True, struct_as_record=False)

        for state, result in zip(states, actual["results"]):
            with self.subTest(channel=state):
                self.compare_result(result, expected[state], state)
        self.compare_result(actual["truncated"], truncated, "5S")
        for i, model in enumerate(models):
            with self.subTest(width_model=i):
                for field, reference in (
                    ("width_shapes", model.bw(probe)),
                    ("width_gammas", model.gamma(probe)),
                    ("width_partials", model.partial_widths(probe)),
                    ("width_poles", model.pole_partial_widths),
                ):
                    observed = np.asarray(actual[field][i]).reshape(np.shape(reference))
                    np.testing.assert_allclose(observed, reference, rtol=3e-11, atol=3e-14)

    def compare_result(self, observed, expected, state):
        self.assertEqual(observed.metadata.schema, "eta_c-pp-bw-matlab-v1")
        self.assertEqual(observed.metadata.channel, state)
        self.assertEqual(observed.metadata.alpha_selection, "none")
        self.assertEqual(observed.metadata.rg.loops, 4)
        channel = getattr(observed, "bw" + state.lower())
        reference = expected["channels"][state]
        self.assertFalse(hasattr(channel, "alpha_star"))
        np.testing.assert_allclose(observed.M2, expected["M2"], rtol=2e-15, atol=2e-15)
        np.testing.assert_allclose(observed.weights, expected["weights"], rtol=2e-15, atol=2e-15)
        for key, value in expected["OPE"].items():
            np.testing.assert_allclose(getattr(observed.OPE, key), value,
                                       rtol=3e-12, atol=3e-12)
        np.testing.assert_allclose(
            [observed.metadata.mc, observed.metadata.alpha_s],
            [expected["metadata"]["mc"], expected["metadata"]["alpha_s"]],
            rtol=3e-11, atol=3e-12)
        for key in ("s", "s_hat", "masses", "alpha", "anomaly_shift_GeV", "boundary_total"):
            got = np.asarray(getattr(channel, key)).reshape(np.shape(reference[key]))
            np.testing.assert_allclose(got, reference[key], rtol=3e-12, atol=3e-12)
        for key in ("bw_shape", "gamma", "k_B"):
            got = np.asarray(getattr(channel, key)).reshape(np.shape(reference[key]))
            np.testing.assert_allclose(got, reference[key], rtol=3e-10, atol=3e-13)
        for key in ("residues_all", "f_all", "f_raw_all", "c_all", "hat_all", "rho_all",
                    "bw_hat_all", "total_hat_all", "rho_total_all", "backsub_all", "residual_all",
                    "weighted_residual_all", "relative_weighted_residual_all",
                    "penalty_all", "objective_all"):
            got = np.asarray(getattr(channel, key)).reshape(np.shape(reference[key]))
            np.testing.assert_allclose(got, reference[key], rtol=3e-8, atol=3e-10, equal_nan=True,
                                       err_msg=f"{state}: {key}")
        np.testing.assert_array_equal(np.atleast_1d(channel.positive_all), reference["positive_all"])
        self.assertLess(np.max(np.abs(channel.endpoint_error_all)), 1e-12)
        self.assertLess(np.max(channel.stationarity_relative_all), 1e-8)
        self.assertEqual(channel.factors.n_unregularized, 2)
        self.assertEqual(np.count_nonzero(np.all(channel.factors.B == 0, axis=0)), 2)
        self.assertEqual(channel.factors.n_regularized, len(reference["s"]) - 1)
        np.testing.assert_allclose(channel.factors.kappa, reference["factors"]["kappa"],
                                   rtol=3e-10, atol=3e-13)
        residue_error = np.max(np.abs(np.asarray(channel.residues_all).reshape(
            reference["residues_all"].shape) - reference["residues_all"]))
        kernel_error = np.max(np.abs(channel.k_B-reference["k_B"]))
        print(f"PP BW {state}, Lambda={expected['metadata']['Lambda']:g}: "
              f"max residue difference={residue_error:.3e}, max k_B difference={kernel_error:.3e}")


if __name__ == "__main__":
    unittest.main()
