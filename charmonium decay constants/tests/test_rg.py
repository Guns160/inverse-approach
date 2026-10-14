"""Physical regression tests for RunDec conventions and the charm anchor.

Run with python -m unittest discover -s tests -p test_rg.py -v.
The numeric references use independent 50-digit mpmath quadrature.
"""
from dataclasses import FrozenInstanceError, replace
import math
import json
import os
import subprocess
import tempfile
from pathlib import Path
import sys
import unittest
from scipy.integrate import quad
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inverse_common.rg import (
    RGInputs, beta_coefficients, gamma_coefficients, derive,
    evolve_fixed_nf, decouple_alpha_down,
)

# Independent quadrature with exact rationals and mp.zeta, not DOP853.
# The separately supplied 55-digit physics audit confirms these references.
FOUR_LOOP = {
    "alpha_s": 0.304398404402813591292861525129,
    "mc": 1.093276847184784039788804509884,
    "alpha5_mb": 0.226141845816242853648784456087,
    "alpha4_mb": 0.226374296730276340625973984089,
    "alpha_anchor": 0.391619277160514744523610869383,
}

class TestQCDRG(unittest.TestCase):

    def assertClose(self, actual, expected, *, rel_tol=0.0, abs_tol=0.0):
        """Compare scalars or short sequences with explicit numeric tolerances."""
        if isinstance(expected, (tuple, list)):
            self.assertEqual(len(actual), len(expected))
            for (value, reference) in zip(actual, expected):
                self.assertClose(value, reference, rel_tol=rel_tol, abs_tol=abs_tol)
        else:
            self.assertTrue(math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol), f'{actual!r} != {expected!r} (rel_tol={rel_tol}, abs_tol={abs_tol})')

    def test_three_loop_restores_existing_production_inputs(self):
        r = derive(loops=3)
        self.assertClose(r['alpha_s'], 0.3029450988131294, rel_tol=0, abs_tol=2e-11)
        self.assertClose(r['mc'], 1.0990536546335794, rel_tol=0, abs_tol=2e-11)
        self.assertTrue(r['matching_loops'] == 2)

    def test_four_loop_independent_high_precision_reference(self):
        r = derive()
        self.assertClose(r['alpha_s'], FOUR_LOOP['alpha_s'], rel_tol=0, abs_tol=2e-11)
        self.assertClose(r['mc'], FOUR_LOOP['mc'], rel_tol=0, abs_tol=2e-11)
        self.assertClose(r['threshold']['alpha_s_high'], FOUR_LOOP['alpha5_mb'], rel_tol=0, abs_tol=2e-11)
        self.assertClose(r['threshold']['alpha_s_low'], FOUR_LOOP['alpha4_mb'], rel_tol=0, abs_tol=2e-11)
        self.assertClose(r['mass_anchor']['alpha_s'], FOUR_LOOP['alpha_anchor'], rel_tol=0, abs_tol=2e-11)
        self.assertTrue((r['nf'], r['loops'], r['matching_loops']) == (4, 4, 3))
        self.assertEqual(r["mass_method"], "ode")
        # RunDec Eq. (10) truncates c(a); it is a different higher-order
        # numerical convention, and must not become the ODE regression target.
        truncated_c_mass = 1.092612428731451
        self.assertClose(
            r["mc"] - truncated_c_mass, 0.00066441845333304,
            abs_tol=2e-11,
        )

    def test_coefficients_use_alpha_over_pi_and_log_mu_squared(self):
        for (nf, beta, gamma) in [(4, (2.0833333333333335, 3.2083333333333335, 6.349247685185185, 31.38744695230514), (1.0, 3.6527777777777777, 9.947040261072956, 27.302933634218547)), (5, (1.9166666666666667, 2.4166666666666665, 2.826678240740741, 18.85217315933944), (1.0, 3.513888888888889, 7.419863212143666, 11.03432124206806))]:
            with self.subTest(nf=nf, beta=beta, gamma=gamma):
                self.assertClose(beta_coefficients(nf), beta, rel_tol=2e-14)
                self.assertClose(gamma_coefficients(nf), gamma, rel_tol=2e-14)

    def test_one_loop_against_analytic_solution(self):
        (alpha0, mass0, start, end, nf) = (0.28, 1.2, 2.0, 3.0, 4)
        b0 = (11 - 2 * nf / 3) / 4
        a0 = alpha0 / math.pi
        a1 = a0 / (1 + b0 * a0 * 2 * math.log(end / start))
        expected = (math.pi * a1, mass0 * (a1 / a0) ** (1 / b0))
        self.assertClose(evolve_fixed_nf(alpha0, mass0, start, end, nf, 1), expected, rel_tol=3e-12)

    def test_rg_round_trip_and_semigroup(self):
        for loops in [1, 2, 3, 4]:
            with self.subTest(loops=loops):
                first = evolve_fixed_nf(0.3, 1.1, 2.0, 3.0, 4, loops)
                back = evolve_fixed_nf(*first, 3.0, 2.0, 4, loops)
                self.assertClose(back, (0.3, 1.1), rel_tol=5e-12)
                staged = evolve_fixed_nf(*first, 3.0, 4.18, 4, loops)
                direct = evolve_fixed_nf(0.3, 1.1, 2.0, 4.18, 4, loops)
                self.assertClose(staged, direct, rel_tol=5e-12)

    def test_matching_at_msbar_self_scale_and_flavour_count(self):
        alpha = FOUR_LOOP['alpha5_mb']
        a = alpha / math.pi
        d2 = 11 / 72
        d3 = 0.6334507195485173
        self.assertTrue(decouple_alpha_down(alpha, matching_loops=0) == alpha)
        self.assertTrue(decouple_alpha_down(alpha, matching_loops=1) == alpha)
        self.assertClose(decouple_alpha_down(alpha, matching_loops=2), alpha * (1 + d2 * a * a), rel_tol=0, abs_tol=1e-15)
        self.assertClose(decouple_alpha_down(alpha, matching_loops=3), alpha * (1 + d2 * a * a + d3 * a ** 3), rel_tol=0, abs_tol=1e-15)
        self.assertTrue(decouple_alpha_down(alpha, nf=4) > decouple_alpha_down(alpha, nf=5))
        r = derive(mu=4.18)
        self.assertClose(r['alpha_s'], FOUR_LOOP['alpha4_mb'], rel_tol=0, abs_tol=2e-11)
        self.assertTrue(r['alpha_s'] > r['threshold']['alpha_s_high'])
        self.assertTrue(r['nf'] == 4)

    def test_mass_anchor_moves_with_input_and_matches_independent_quadrature(self):
        for mcbar in [1.245, 1.275, 1.305]:
            with self.subTest(mcbar=mcbar):
                inp = replace(RGInputs(), mbar_c=mcbar)
                anchored = derive(mu=mcbar, inputs=inp)
                self.assertTrue(anchored['mc'] == mcbar)
                self.assertTrue(anchored['nf'] == anchored['mass_anchor']['nf'] == 4)
                r = derive(inputs=inp)
                # Independent coupling-space quadrature: do not call the
                # module coefficient/evolution routines for this reference.
                b = (2.0833333333333335, 3.2083333333333335, 6.349247685185185, 31.38744695230514)
                g = (1.0, 3.6527777777777777, 9.947040261072956, 27.302933634218547)
                B = lambda a: sum((x * a ** (i + 2) for (i, x) in enumerate(b)))
                G = lambda a: sum((x * a ** (i + 1) for (i, x) in enumerate(g)))
                a_mb = FOUR_LOOP['alpha4_mb'] / math.pi
                residual = lambda a: quad(lambda x: 1 / B(x), a_mb, a, epsabs=1e-12, epsrel=1e-12)[0] + 2 * math.log(mcbar / 4.18)
                a_anchor = brentq(residual, a_mb, 0.3, xtol=5e-15)
                a_target = FOUR_LOOP['alpha_s'] / math.pi
                expected_mass = mcbar * math.exp(quad(lambda x: G(x) / B(x), a_anchor, a_target, epsabs=1e-12, epsrel=1e-12)[0])
                self.assertTrue(r['mass_anchor']['mu'] == mcbar)
                self.assertClose(r['mc'], expected_mass, rel_tol=0, abs_tol=2e-11)
                self.assertClose(r['alpha_s'], FOUR_LOOP['alpha_s'], rel_tol=0, abs_tol=2e-11)

    def test_loop_order_changes_numbers_not_just_metadata(self):
        rows = [derive(loops=n) for n in range(1, 5)]
        self.assertTrue(len({r['mc'] for r in rows}) == 4)
        self.assertTrue(len({r['alpha_s'] for r in rows}) == 4)
        self.assertTrue(rows[2]['mc'] - rows[3]['mc'] > 0.005)
        self.assertTrue(rows[3]['alpha_s'] - rows[2]['alpha_s'] > 0.001)
        self.assertTrue([r['matching_loops'] for r in rows] == [0, 1, 2, 3])

    def test_immutable_inputs_and_fresh_cached_metadata(self):
        inp = RGInputs()
        with self.assertRaises(FrozenInstanceError):
            inp.mbar_c = 1.3
        r = derive(inputs=inp)
        r['inputs']['mbar_c'] = 999
        r['mass_anchor']['alpha_s'] = 999
        r['threshold']['alpha_s_low'] = 999
        fresh = derive(inputs=inp)
        self.assertTrue(fresh['inputs']['mbar_c'] == 1.275)
        self.assertTrue(fresh['mass_anchor']['alpha_s'] < 1)
        self.assertTrue(fresh['threshold']['alpha_s_low'] < 1)
        changed = derive(inputs=replace(inp, alpha_mz=0.119))
        self.assertTrue(changed['alpha_s'] > fresh['alpha_s'])
        self.assertTrue(changed['mc'] < fresh['mc'])
        self.assertTrue(changed['input_source'] == 'Caller-supplied RGInputs')

    def test_rejects_invalid_or_unmatched_target_domains(self):
        for mu in [1.274, 4.181, 100.0, 0.0, -1.0, math.nan, math.inf, True, '2']:
            with self.subTest(mu=mu):
                with self.assertRaises(ValueError):
                    derive(mu=mu)

    def test_rejects_invalid_loop_orders(self):
        for loops in [0, 5, 3.0, True, math.nan]:
            with self.subTest(loops=loops):
                with self.assertRaises(ValueError):
                    derive(loops=loops)

    def test_rejects_invalid_boundary_inputs(self):
        for kwargs in [{'mbar_c': 0}, {'mb': 1}, {'mz': 4}, {'alpha_mz': 0}, {'alpha_mz': math.pi}, {'mbar_c': math.nan}, {'mb': True}]:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    RGInputs(**kwargs)

    def test_low_level_validation_and_nonperturbative_failure(self):
        with self.assertRaises(TypeError):
            derive(inputs={'mbar_c': 1.275})
        for kwargs in ({'nf': 0}, {'nf': 4.0}, {'mass': 0}, {'alpha_s': math.pi}):
            args = dict(alpha_s=0.3, mass=1.1, start=2.0, end=3.0, nf=4)
            args.update(kwargs)
            with self.assertRaises(ValueError):
                evolve_fixed_nf(**args)
        with self.assertRaises(ValueError):
            decouple_alpha_down(0.2, matching_loops=4)
        with self.assertRaises(ValueError):
            decouple_alpha_down(0.2, nf=1)
        with self.assertRaisesRegex(RuntimeError, 'perturbative guard'):
            evolve_fixed_nf(0.3, 1.1, 2.0, 0.01, 4, 4)

class TestRGCLI(unittest.TestCase):
    def run_cli(self, *args, direct=False, cwd=None):
        root = Path(__file__).resolve().parents[1]
        entry = [str(root / 'inverse_common' / 'rg.py')] if direct else ['-m', 'inverse_common.rg']
        return subprocess.run([sys.executable, *entry, *args], cwd=cwd or root,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'),
            text=True, capture_output=True, timeout=30)

    def test_module_defaults_emit_json(self):
        completed = self.run_cli()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stderr, '')
        r = json.loads(completed.stdout)
        self.assertEqual((r['mu'], r['loops'], r['matching_loops']), (2., 4, 3))
        self.assertEqual(r['mass_method'], 'ode')
        self.assertAlmostEqual(r['alpha_s'], FOUR_LOOP['alpha_s'], delta=2e-11)
        self.assertAlmostEqual(r['mc'], FOUR_LOOP['mc'], delta=2e-11)

    def test_module_selects_all_loop_orders(self):
        for loops in range(1, 5):
            with self.subTest(loops=loops):
                completed = self.run_cli('--mu', '2', '--loops', str(loops))
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(json.loads(completed.stdout), derive(loops=loops))

    def test_direct_script_custom_anchors_writes_no_files(self):
        with tempfile.TemporaryDirectory() as directory:
            completed = self.run_cli('--mu', '2', '--loops', '4', '--mcbar', '1.305',
                '--alpha-mz', '0.119', '--mb', '4.2', '--mz', '91.2', direct=True, cwd=directory)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(completed.stderr, '')
            expected = derive(inputs=RGInputs(mbar_c=1.305, alpha_mz=.119, mb=4.2, mz=91.2))
            self.assertEqual(json.loads(completed.stdout), expected)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_input_errors_are_nonzero_clear_and_no_stdout(self):
        cases = [
            (('--loops', '5'), 'invalid choice'),
            (('--loops', '3.5'), 'invalid int value'),
            (('--mu', '5'), 'mu must lie'),
            (('--mcbar', 'nan'), 'mbar_c must be'),
            (('--alpha-mz', '-1'), 'alpha_s must be'),
            (('--mb', '1'), 'require mbar_c < mb < mz'),
            (('--mu', 'not-a-number'), 'invalid float value'),
            (('--MassScale', '1.3'), 'unrecognized arguments'),
        ]
        for args, message in cases:
            with self.subTest(args=args):
                completed = self.run_cli(*args)
                self.assertEqual(completed.returncode, 2)
                self.assertEqual(completed.stdout, '')
                self.assertIn(message, completed.stderr)
                self.assertNotIn('Traceback', completed.stderr)

    def test_integration_failure_returns_code_one(self):
        completed = self.run_cli('--alpha-mz', '0.5')
        self.assertEqual(completed.returncode, 1)
        self.assertEqual(completed.stdout, '')
        self.assertIn('perturbative guard', completed.stderr)
        self.assertNotIn('Traceback', completed.stderr)

    def test_help_has_no_calculation_result(self):
        completed = self.run_cli('--help', direct=True)
        self.assertEqual(completed.returncode, 0)
        self.assertIn('--loops', completed.stdout)
        self.assertIn('--mcbar', completed.stdout)
        self.assertEqual(completed.stderr, '')

if __name__ == '__main__':
    unittest.main()
