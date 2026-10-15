from pathlib import Path
import csv
import json
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from uncertainty.workflow import analyze, save_analysis


class SeparatedWorkflowTests(unittest.TestCase):
    def config(self):
        return json.loads((ROOT / "examples/toy_config.json").read_text())

    def test_statistical_distribution_and_fixed_systematic_parameters(self):
        result = analyze(self.config(), ensemble="statistical", n=4000,
                         base_directory=ROOT / "examples")
        self.assertEqual(result["design"]["fixed_parameters"],
                         ["scale", "LambdaBar", "Y1", "Y2", "X2", "hx"])
        report = result["uncertainty_reporting"]
        self.assertEqual(report["kind"], "statistical_distribution")
        self.assertIn("std", report["outputs"]["sum_response"])
        self.assertIn("q16", report["outputs"]["sum_response"])
        self.assertIn("no_total_uncertainty_is_computed", result["combination_policy"])
        for row in result["samples"]:
            self.assertEqual(row["inputs"]["Y1"], 3.)

    def test_systematic_envelope_and_fixed_statistical_parameters(self):
        result = analyze(self.config(), ensemble="systematic", n=2000,
                         base_directory=ROOT / "examples")
        report = result["uncertainty_reporting"]
        self.assertEqual(report["kind"], "systematic_envelope")
        self.assertIn("envelope_minus", report["outputs"]["window_width"])
        self.assertNotIn("std", report["outputs"]["window_width"])
        for row in result["samples"]:
            self.assertEqual(row["inputs"]["Q2"], 2.)
            self.assertEqual(row["inputs"]["condensate"], 1.)

    def test_cli_runs_ensembles_separately(self):
        with tempfile.TemporaryDirectory() as temporary:
            outputs = []
            for ensemble, seed in (("statistical", 15001), ("systematic", 15002)):
                command = [sys.executable, "-B", str(ROOT / "scripts/run_uncertainty.py"),
                           "--config", str(ROOT / "examples/toy_config.json"),
                           "--ensemble", ensemble, "--samples", "24", "--seed", str(seed),
                           "--output", temporary]
                process = subprocess.run(command, cwd=temporary, text=True, capture_output=True)
                self.assertEqual(process.returncode, 0, process.stderr)
                summary = json.loads(process.stdout)
                self.assertEqual(summary["ensemble"], ensemble)
                self.assertEqual(summary["summary"]["n_valid"], 24)
                self.assertIn("no_total_uncertainty_is_computed", summary["combination_policy"])
                outputs.append(summary)
            self.assertNotEqual(outputs[0]["result_directory"], outputs[1]["result_directory"])
            with (Path(outputs[0]["result_directory"]) / "samples.csv").open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 24)

    def test_configuration_rejects_old_combined_schema(self):
        config = self.config()
        config["variables"] = config["ensembles"]["statistical"]["variables"]
        config["dependence"] = {"kind": "independent"}
        with self.assertRaisesRegex(ValueError, "configuration"):
            analyze(config, ensemble="statistical", n=3, base_directory=ROOT / "examples")

    def test_all_invalid_results_still_create_known_output_columns(self):
        def evaluate(point):
            if point["Y1"] == 3.0 and point["Y2"] == 15.0 and point["X2"] == 2.0:
                return {"sum_response": 0., "window_width": 12., "scaled_product": 0.}
            raise ValueError("invalid synthetic model domain")

        config = self.config()
        result = analyze(config, ensemble="systematic", n=2, seed=15002,
                         evaluator=evaluate, base_directory=ROOT / "examples")
        self.assertEqual(result["summary"]["n_valid"], 0)
        with tempfile.TemporaryDirectory() as temporary:
            saved = save_analysis(result, temporary)
            with (Path(saved["result_directory"]) / "samples.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["output:window_width"], "")
            self.assertIn("invalid synthetic model domain", rows[0]["reason"])


if __name__ == "__main__":
    unittest.main()
