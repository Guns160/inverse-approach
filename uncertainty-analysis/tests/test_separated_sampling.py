from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from uncertainty.sampling import generate_draws, run_sampling, OutputSchemaError


class SeparatedSamplingTests(unittest.TestCase):
    def config(self):
        return {
            "central": {"Q2": 2., "condensate": 1., "m0": .8, "LambdaBar": .55,
                        "Y1": 3., "Y2": 15., "X2": 2., "hx": .01},
            "stat": {
                "variables": {"Q2": {"kind": "normal", "sigma": .2},
                              "condensate": {"kind": "normal", "sigma": .1},
                              "m0": {"kind": "uniform", "low": .6, "high": 1.}},
                "dependence": {"kind": "normal_covariance", "names": ["Q2", "condensate"],
                               "matrix": [[.04, .01], [.01, .01]], "independent_remaining": True},
                "constraints": [{"kind": "positive", "names": ["Q2", "condensate", "m0", "LambdaBar"]},
                                {"kind": "less_than_margin", "left": "LambdaBar", "right": "X2", "margin_name": "hx"}],
            },
            "sys": {
                "variables": {"Y1": {"kind": "uniform", "low": 2., "high": 4.},
                              "Y2": {"kind": "uniform", "low": 10., "high": 15.},
                              "X2": {"kind": "uniform", "low": 2., "high": 3.}},
                "dependence": {"kind": "independent"},
                "constraints": [{"kind": "min_difference", "lower": "Y1", "upper": "Y2", "minimum": 0.},
                                {"kind": "greater_than_margin", "left": "X2", "right": "LambdaBar", "margin_name": "hx"}],
            },
        }

    def statistical(self, c):
        return generate_draws(c["central"], c["stat"]["variables"], ensemble="statistical",
                              n=2000, seed=91, dependence=c["stat"]["dependence"],
                              constraints=c["stat"]["constraints"])

    def systematic(self, c):
        return generate_draws(c["central"], c["sys"]["variables"], ensemble="systematic",
                              n=2000, seed=92, dependence=c["sys"]["dependence"],
                              constraints=c["sys"]["constraints"])

    def test_statistical_fixes_all_systematic_inputs(self):
        rows = self.statistical(self.config())["draws"]
        self.assertEqual(len(rows), 2000)
        for key in ("Y1", "Y2", "X2"):
            self.assertEqual({row[key] for row in rows}, {self.config()["central"][key]})
        self.assertTrue(all(row["LambdaBar"] < row["X2"] - row["hx"] for row in rows))
        self.assertTrue(all(.6 <= row["m0"] <= 1. for row in rows))

    def test_systematic_fixes_all_statistical_inputs(self):
        rows = self.systematic(self.config())["draws"]
        for key in ("Q2", "condensate", "m0"):
            self.assertEqual({row[key] for row in rows}, {self.config()["central"][key]})
        self.assertTrue(all(row["Y2"] > row["Y1"] and row["X2"] > row["LambdaBar"] + row["hx"] for row in rows))

    def test_rejection_is_recorded_and_seed_is_reproducible(self):
        c = self.config()
        c["stat"]["variables"]["Q2"]["sigma"] = 2.0
        c["stat"]["dependence"]["matrix"][0][0] = 4.0
        one = self.statistical(c)
        two = self.statistical(c)
        self.assertEqual(one["draws"], two["draws"])
        self.assertGreaterEqual(one["design"]["n_candidate_inputs_attempted"], 2000)
        self.assertGreater(one["design"]["input_rejection_fraction"], 0)

    def test_model_failures_are_not_input_resampling(self):
        c = self.config()
        supplied = [{"Q2": -1., "condensate": 1., "m0": .8},
                    {"Q2": -2., "condensate": 1., "m0": .8}]
        c["stat"]["variables"] = {name: {"kind": "provided"} for name in supplied[0]}
        c["stat"]["constraints"] = None
        def invalid_for_negative_q2(point):
            if point["Q2"] < 0:
                raise ValueError("negative synthetic input")
            return point["Q2"]

        result = run_sampling(invalid_for_negative_q2, c["central"], c["stat"]["variables"],
                              ensemble="statistical", n=2, dependence={"kind": "provided_paired_rows"},
                              supplied_draws=supplied)
        self.assertEqual(result["summary"]["n_invalid"], 2)
        self.assertEqual(len(result["samples"]), 2)
        self.assertEqual(result["samples"][0]["inputs"]["Y1"], c["central"]["Y1"])

    def test_output_schema_and_programming_errors_are_not_physical_failures(self):
        c = self.config()
        with self.assertRaises(OutputSchemaError):
            run_sampling(lambda p: {"a": [1.]}, c["central"], c["stat"]["variables"],
                         ensemble="statistical", n=2, dependence={"kind": "independent"})
        def bug(p):
            raise KeyError("adapter bug")
        with self.assertRaises(KeyError):
            run_sampling(bug, c["central"], c["stat"]["variables"],
                         ensemble="statistical", n=2, dependence={"kind": "independent"})

    def test_invalid_ensemble_or_hidden_dependence_is_rejected(self):
        c = self.config()
        with self.assertRaises(ValueError):
            generate_draws(c["central"], c["stat"]["variables"], ensemble="total", n=2,
                           dependence={"kind": "independent"})
        with self.assertRaises(ValueError):
            generate_draws(c["central"], c["stat"]["variables"], ensemble="statistical", n=2)


if __name__ == "__main__":
    unittest.main()

