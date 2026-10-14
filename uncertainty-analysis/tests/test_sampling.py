from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from uncertainty.sampling import generate_draws, run_sampling, OutputSchemaError


class SamplingTests(unittest.TestCase):
    def setUp(self):
        self.center = {"Q2": 2., "condensate": 1., "fixed": 7.}
        self.variables = {"Q2": {"kind": "normal", "sigma": .2},
                          "condensate": {"kind": "normal", "sigma": .1}}

    def test_single_random_only_changes_one_input(self):
        d = generate_draws(self.center, self.variables, mode="single",
                           parameter="Q2", n=100, seed=9)
        self.assertEqual(len({row["Q2"] for row in d["draws"]}), 100)
        self.assertTrue(all(row["condensate"] == 1 and row["fixed"] == 7 for row in d["draws"]))
        self.assertEqual(d, generate_draws(self.center, self.variables, mode="single",
                                          parameter="Q2", n=100, seed=9))

    def test_joint_covariance_sign_changes_output_variance(self):
        for sign, expected in [(1, .07), (-1, .03)]:
            with self.subTest(sign=sign):
                dependence = {"kind": "normal_covariance", "names": ["Q2", "condensate"],
                              "matrix": [[.04, sign*.01], [sign*.01, .01]]}
                r = run_sampling(lambda p: p["Q2"] + p["condensate"], self.center, self.variables,
                                 mode="joint", n=20000, seed=91, dependence=dependence)
                std = r["summary"]["outputs"]["output"]["std"]
                self.assertAlmostEqual(std**2, expected, delta=.03*expected)

    def test_joint_requires_explicit_dependence(self):
        with self.assertRaisesRegex(ValueError, "explicit dependence"):
            generate_draws(self.center, self.variables, mode="joint", n=3)

    def test_psd_singular_and_small_unit_covariances(self):
        center = {"a": 0., "b": 0.}
        variables = {"a": {"kind": "normal", "sigma": 1e-13},
                     "b": {"kind": "normal", "sigma": 2e-13}}
        good = {"kind": "normal_covariance", "names": ["a", "b"],
                "matrix": [[1e-26, 2e-26], [2e-26, 4e-26]]}
        draws = generate_draws(center, variables, mode="joint", n=20, dependence=good)["draws"]
        np.testing.assert_allclose([row["b"] for row in draws], [2*row["a"] for row in draws], atol=1e-27)
        bad = dict(good, matrix=[[1e-26, 3e-26], [3e-26, 4e-26]])
        with self.assertRaisesRegex(ValueError, "positive semidefinite"):
            generate_draws(center, variables, mode="joint", n=20, dependence=bad)

    def test_bad_covariance_and_diagonal_are_rejected(self):
        for matrix in ([[.04, .01], [.02, .01]], [[.05, 0], [0, .01]], [[.04, .05], [.05, .01]]):
            with self.assertRaises(ValueError):
                generate_draws(self.center, self.variables, mode="joint", n=3, dependence={
                    "kind": "normal_covariance", "names": ["Q2", "condensate"], "matrix": matrix})

    def test_uniform_and_loguniform_with_declared_independence(self):
        center = {"M2": 2., "alpha": 1e-12}
        variables = {"M2": {"kind": "uniform", "low": 1., "high": 3.},
                     "alpha": {"kind": "loguniform", "low": 1e-14, "high": 1e-10}}
        rows = generate_draws(center, variables, mode="joint", n=3000,
                              dependence={"kind": "independent"})["draws"]
        self.assertTrue(all(1 <= row["M2"] <= 3 and 1e-14 <= row["alpha"] <= 1e-10 for row in rows))
        self.assertAlmostEqual(np.mean([np.log10(row["alpha"]) for row in rows]), -12., delta=.08)

    def test_normal_block_requires_explicit_independent_remainder(self):
        variables = dict(self.variables, fixed={"kind": "uniform", "low": 6., "high": 8.})
        dep = {"kind": "normal_covariance", "names": ["Q2", "condensate"],
               "matrix": [[.04, 0], [0, .01]]}
        with self.assertRaisesRegex(ValueError, "independent_remaining"):
            generate_draws(self.center, variables, mode="joint", n=3, dependence=dep)
        generate_draws(self.center, variables, mode="joint", n=3, dependence=dict(dep, independent_remaining=True))

    def test_supplied_pairs_and_single_marginal(self):
        pairs = [{"Q2": float(i), "condensate": -float(i)} for i in range(5)]
        variables = {name: {"kind": "provided"} for name in self.variables}
        joint = run_sampling(lambda p: {"sum": p["Q2"]+p["condensate"]}, self.center, variables,
                             mode="joint", n=5, supplied_draws=pairs)
        self.assertTrue(all(row["values"] == [0.] for row in joint["samples"]))
        self.assertFalse(joint["design"]["seed_used"])
        single = generate_draws(self.center, variables, mode="single", parameter="Q2",
                                n=5, supplied_draws=pairs)["draws"]
        self.assertEqual([row["Q2"] for row in single], list(range(5)))
        self.assertTrue(all(row["condensate"] == 1. for row in single))
        self.assertEqual(joint["design"]["dependence_used"], "provided_paired_rows")

    def test_supplied_draws_reject_misspelled_or_conflicting_dependence(self):
        pairs = [{"Q2": 2., "condensate": 1.}] * 3
        for dependence in ({"kind": "typo"}, {"kind": "independent"}):
            with self.assertRaisesRegex(ValueError, "paired rows"):
                generate_draws(self.center, self.variables, mode="joint", n=3,
                               supplied_draws=pairs, dependence=dependence)
        good = generate_draws(self.center, self.variables, mode="joint", n=3,
                              supplied_draws=pairs, dependence={"kind": "provided_paired_rows"})
        self.assertEqual(good["design"]["sample_input_covariance_names"], list(self.variables))

    def test_three_outputs_and_invalid_pairs_are_not_resampled(self):
        def model(p):
            if p["x"] < 0:
                raise ValueError("outside model domain")
            return {"a": p["x"], "b": 2*p["x"], "c": 3*p["x"]}
        r = run_sampling(model, {"x": 1.}, {"x": {"kind": "provided"}}, mode="joint", n=4,
                         supplied_draws=[{"x": x} for x in [-2, -1, 1, 2]])
        self.assertEqual(r["summary"]["n_valid"], 2)
        self.assertEqual(r["summary"]["n_invalid"], 2)
        self.assertTrue(r["summary"]["conditional_on_valid"])
        self.assertEqual(np.asarray(r["summary"]["covariance"]).shape, (3, 3))
        self.assertEqual(len(r["samples"]), 4)
        self.assertEqual(r["samples"][0]["values"], [None, None, None])

    def test_nonfinite_outputs_invalidate_the_whole_pair(self):
        r = run_sampling(lambda p: [p["x"], np.nan if p["x"] < 0 else 1.],
                         {"x": 0.}, {"x": {"kind": "provided"}}, mode="joint", n=3,
                         supplied_draws=[{"x": x} for x in [-1, 0, 1]])
        self.assertEqual(r["summary"]["n_invalid"], 1)
        self.assertEqual(r["summary"]["outputs"]["output_0"]["mean"], .5)

    def test_all_failed_with_known_output_schema(self):
        def bad(p):
            raise RuntimeError("no solution")
        r = run_sampling(bad, self.center, self.variables, mode="single", parameter="Q2",
                         n=3, output_names=["a", "b", "c"])
        self.assertEqual(r["summary"]["status"], "insufficient_valid_draws")
        self.assertIsNone(r["summary"]["outputs"]["a"]["std"])
        self.assertIsNone(r["summary"]["covariance"])

    def test_output_schema_changes_and_programming_errors_raise(self):
        with self.assertRaises(OutputSchemaError):
            run_sampling(lambda p: [1., 2.], self.center, self.variables,
                         mode="single", parameter="Q2", n=3, output_names=["only_one"])
        def bug(p):
            raise KeyError("coding bug")
        with self.assertRaises(KeyError):
            run_sampling(bug, self.center, self.variables, mode="single", parameter="Q2", n=3)

    def test_skewed_interval_is_not_forced_to_equal_standard_deviation(self):
        r = run_sampling(lambda p: np.exp(p["x"]), {"x": 0.},
                         {"x": {"kind": "normal", "sigma": 1.}},
                         mode="single", parameter="x", n=5000, seed=4)
        s = r["summary"]["outputs"]["output"]
        self.assertGreater(s["plus_from_median"], s["minus_from_median"])
        self.assertGreater(abs(s["std"] - s["plus_from_median"]), .1)

    def test_callback_mutation_does_not_rewrite_input_records(self):
        def mutate(p):
            value = p["fixed"]
            p["fixed"] = 999
            return value
        r = run_sampling(mutate, self.center, self.variables, mode="single",
                         parameter="Q2", n=4)
        self.assertEqual(self.center["fixed"], 7.)
        self.assertTrue(all(row["inputs"]["fixed"] == 7. for row in r["samples"]))

    def test_invalid_configuration_does_not_fall_back(self):
        for n in [1, 0, 3.5, True]:
            with self.assertRaises(ValueError):
                generate_draws(self.center, self.variables, mode="joint", n=n,
                               dependence={"kind": "independent"})
        with self.assertRaises(ValueError):
            generate_draws(self.center, self.variables, mode="joint", n=3,
                           dependence={"kind": "independent", "correlation": .5})
        with self.assertRaises(OutputSchemaError):
            run_sampling(lambda p: {"a": [1.]}, self.center, self.variables,
                         mode="single", parameter="Q2", n=3)


if __name__ == "__main__":
    unittest.main()
