from pathlib import Path
import csv
import json
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from uncertainty.workflow import analyze, load_evaluator, save_analysis


class WorkflowTests(unittest.TestCase):
    def config(self):
        return json.loads((ROOT / "examples/toy_config.json").read_text())

    def test_named_general_outputs_single_and_joint(self):
        config = self.config()
        single = analyze(config, mode="single", parameter="Q2", n=1000, seed=13,
                         base_directory=ROOT / "examples")
        self.assertEqual(len(single["output_names"]), 3)
        for row in single["samples"]:
            self.assertEqual(row["inputs"]["C_example"], 1.0)
            self.assertEqual(row["inputs"]["scale"], 1.0)
        joint = analyze(config, mode="joint", n=10000, seed=13,
                        base_directory=ROOT / "examples")
        a = np.asarray([row["values"] for row in joint["samples"]])
        # Var(x+y)=.04+.01+2*.01=.07; Var(x-y)=.03.
        np.testing.assert_allclose(a.var(axis=0, ddof=1)[:2], [.07, .03], rtol=.06)
        self.assertEqual(joint["summary"]["n_valid"], 10000)
        self.assertEqual(joint["central_result"],
                         {"sum_response": 3.0, "difference_response": 1.0, "scaled_product": 2.0})

    def test_save_every_draw_and_unique_output(self):
        result = analyze(self.config(), mode="single", parameter="Q2", n=8,
                         base_directory=ROOT / "examples")
        with tempfile.TemporaryDirectory() as temporary:
            one = save_analysis(result, temporary)
            two = save_analysis(result, temporary)
            self.assertNotEqual(one["result_directory"], two["result_directory"])
            folder = Path(one["result_directory"])
            with (folder / "samples.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 8)
            self.assertIn("input:Q2", rows[0])
            self.assertIn("output:sum_response", rows[0])
            raw = (folder / "summary.json").read_text()
            self.assertNotIn("NaN", raw)
            self.assertEqual(json.loads(raw), one)

    def test_arbitrary_callable_and_baseline_validation(self):
        config = {"central": {"temperature": 2.0},
                  "variables": {"temperature": {"kind": "uniform", "low": 1, "high": 3}},
                  "dependence": {"kind": "independent"}}
        out = analyze(config, evaluator=lambda p: p["temperature"]**2,
                      mode="joint", n=4)
        self.assertEqual(out["central_result"], 4.0)
        self.assertEqual(len(out["output_names"]), 1)
        with self.assertRaises(ValueError):
            analyze(config, evaluator=lambda p: np.nan, mode="joint", n=4)

    def test_cli_from_other_working_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            command = [sys.executable, "-B", str(ROOT / "scripts/run_uncertainty.py"),
                       "--config", str(ROOT / "examples/toy_config.json"),
                       "--mode", "joint", "--samples", "12", "--seed", "42",
                       "--output", temporary]
            process = subprocess.run(command, cwd=temporary, text=True, capture_output=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            result = json.loads(process.stdout)
            self.assertEqual(result["summary"]["n_valid"], 12)
            self.assertTrue((Path(result["result_directory"]) / "samples.csv").is_file())

    def test_loader_rejects_ambiguous_source(self):
        with self.assertRaises(ValueError):
            load_evaluator({"file": "x.py", "module": "x"})

    def test_all_failed_draws_are_saved_with_known_columns(self):
        config = {"central": {"x": 0.}, "variables": {"x": {"kind": "provided"}}}
        def evaluate(p):
            if p["x"] != 0:
                raise ValueError("outside valid branch")
            return {"a": 1., "b": 2.}
        result = analyze(config, evaluator=evaluate, mode="joint", n=2,
                         supplied_draws=[{"x": 1.}, {"x": 2.}])
        self.assertEqual(result["summary"]["n_valid"], 0)
        self.assertIsNone(result["summary"]["outputs"]["a"]["std"])
        with tempfile.TemporaryDirectory() as temporary:
            saved = save_analysis(result, temporary)
            with (Path(saved["result_directory"]) / "samples.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["output:a"], "")
            self.assertIn("outside valid branch", rows[0]["reason"])


if __name__ == "__main__":
    unittest.main()
