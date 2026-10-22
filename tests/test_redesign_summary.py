"""Decision, seed, terminal-barrier and donor/count interpretation guards."""
from pathlib import Path
import copy
import json
import sys
import tempfile
import unittest

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import summarize_redesign_experiment as summary


def comparison():
    result = {m: {'gain_in_better_direction': 0., 'finite_cohort_identical': True} for m in summary.METRICS}
    result[summary.PRIMARY[0]]['gain_in_better_direction'] = .01
    return result


class RedesignSummaryTests(unittest.TestCase):
    def test_all_six_protected_and_primary_gain_exceeds_tolerance(self):
        good = comparison()
        self.assertTrue(summary.point_gate(good)['passed'])
        for metric in summary.METRICS:
            bad = copy.deepcopy(good)
            bad[metric]['gain_in_better_direction'] = -2e-8
            self.assertFalse(summary.point_gate(bad)['passed'])
        for gain in (0., 1e-12, 1e-8):
            tiny = comparison(); tiny[summary.PRIMARY[0]]['gain_in_better_direction'] = gain
            self.assertFalse(summary.point_gate(tiny)['passed'])
        boundary = comparison(); boundary[summary.MSE]['gain_in_better_direction'] = -1e-8
        self.assertTrue(summary.point_gate(boundary)['passed'])

    def test_missing_nonfinite_or_changed_cohort_cannot_pass(self):
        missing = comparison(); missing.pop(summary.MSE)
        with self.assertRaises(ValueError): summary.point_gate(missing)
        for metric in summary.METRICS:
            changed = comparison(); changed[metric]['finite_cohort_identical'] = False
            self.assertFalse(summary.point_gate(changed)['passed'])
            invalid = comparison(); invalid[metric]['gain_in_better_direction'] = np.nan
            self.assertFalse(summary.point_gate(invalid)['passed'])

    def test_both_seeds_required_without_winner_seed_selection(self):
        with self.assertRaises(ValueError): summary.recipe_gate({summary.SEEDS[0]: {'passed': True}})
        self.assertFalse(summary.recipe_gate({summary.SEEDS[0]: {'passed': True}, summary.SEEDS[1]: {'passed': False}}))
        self.assertTrue(summary.recipe_gate({seed: {'passed': True} for seed in summary.SEEDS}))

    def test_all_terminal_barrier_rejects_pending_and_wrong_arm(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            names = [summary.arm_name(r, s) for r in summary.RECIPES for s in summary.SEEDS]
            for name in names:
                folder = root/name; folder.mkdir()
                (folder/'terminal.json').write_text(json.dumps({'exit_code': 0, 'arm': name}))
            self.assertEqual(len(summary.require_all_terminals(root)), 8)
            terminal = root/names[-1]/'terminal.json'
            terminal.write_text(json.dumps({'exit_code': None, 'arm': names[-1]}))
            with self.assertRaisesRegex(ValueError, 'successful terminal'): summary.require_all_terminals(root)
            terminal.write_text(json.dumps({'exit_code': 0, 'arm': 'wrong'}))
            with self.assertRaises(ValueError): summary.require_all_terminals(root)

    def test_donor_rows_and_depths_rederived_from_frozen_seed(self):
        controls = sparse.csr_matrix(np.tile(np.array([2, 3, 5], dtype=np.int32), (500, 1)))
        seed = summary.SEEDS[0]; target = 'T'
        donor_seed = summary.atlas.seed_for('donor:H1:'+target, seed)
        rows = np.random.default_rng(donor_seed).choice(500, 400, replace=False)
        selected = controls[rows].toarray()
        audit = {'donor_seed': donor_seed, 'emission_seed': summary.atlas.seed_for('transport:H1:'+target, seed),
                 'donor_rows': rows.tolist(), 'donor_rows_sha256': summary.array_sha256(rows),
                 'donor_count_sha256': summary.array_sha256(selected)}
        np.testing.assert_array_equal(summary.verify_donor_metadata(audit, target, seed, controls, np.full(400, 10)), selected)
        with self.assertRaisesRegex(ValueError, 'library depths'):
            summary.verify_donor_metadata(audit, target, seed, controls, np.full(400, 11))
        changed = copy.deepcopy(audit); changed['donor_rows'][0] = changed['donor_rows'][1]
        with self.assertRaisesRegex(ValueError, 'selection'):
            summary.verify_donor_metadata(changed, target, seed, controls, np.full(400, 10))

    def test_own_exclusion_after_full_axis_normalization_and_whole_axis_report(self):
        p, r, c = np.array([90., 9., 1.]), np.array([80., 15., 5.]), np.array([70., 20., 10.])
        genes = np.array(['T', 'A', 'B'])
        result = summary.mean_errors(p, summary.logbulk(r), summary.logbulk(c), genes, 'T')
        squared = (summary.logbulk(p)-summary.logbulk(r))**2
        self.assertAlmostEqual(result['own_excluded']['uncorrected_mean_error'], squared[1:].mean())
        self.assertAlmostEqual(result['whole_axis']['uncorrected_mean_error'], squared.mean())
        self.assertNotAlmostEqual(result['own_excluded']['uncorrected_mean_error'],
                                 ((summary.logbulk(p[1:])-summary.logbulk(r[1:]))**2).mean())


if __name__ == '__main__':
    unittest.main()
