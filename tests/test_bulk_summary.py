"""Independent guard tests for bulk-ablation endpoint interpretation."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import summarize_bulk_ablation as summary


def bundle(nmae=(.8, .9), numerator=(1., 4.), denominator=(2., 8.)):
    rows = []
    series = {m: np.asarray([.5, .6]) for m in summary.PROTECTED}
    series['de_wilcoxon_lfc_nmae'] = np.asarray(nmae)
    series['expr_mse_unbiased_capped'] = np.asarray(numerator)
    series['expr_distance_unbiased'] = np.asarray(denominator)
    for metric, values in series.items():
        rows.extend({'perturbation': target, 'metric': metric, 'value': value}
                    for target, value in zip(('A', 'B'), values))
    means = {m: float(np.nanmean(v)) for m, v in series.items()}
    raw = {m: means[m] for m in summary.PROTECTED}
    raw[summary.MSE] = float(np.nansum(numerator) / np.nansum(denominator))
    return {'raw': raw, 'tidy': pd.DataFrame(rows),
            'aggregate': pd.DataFrame([means], index=['mean']),
            'aggregation': pd.DataFrame({'n_used': {m: np.isfinite(v).sum() for m, v in series.items()}})}


class BulkSummaryTests(unittest.TestCase):
    def test_own_gene_exclusion_follows_whole_axis_normalization(self):
        prediction = summary.logbulk([90., 9., 1.])
        real = summary.logbulk([80., 15., 5.])
        control = summary.logbulk([70., 20., 10.])
        mask = np.asarray([False, True, True])
        actual = summary.profile_error(prediction, real, control, mask)
        expected = ((np.log1p([4500., 500.]) - np.log1p([7500., 2500.])) ** 2).mean()
        self.assertAlmostEqual(actual['uncorrected_mean_error'], expected)
        wrongly_renormalized = ((summary.logbulk([9., 1.]) - summary.logbulk([15., 5.])) ** 2).mean()
        self.assertNotAlmostEqual(actual['uncorrected_mean_error'], wrongly_renormalized)

    def test_global_mse_is_ratio_of_sums_and_has_no_target_bootstrap(self):
        baseline = bundle()
        candidate = bundle(numerator=(.5, 4.))
        result = summary.compare_metrics(baseline, candidate, ['A', 'B'])
        mse = result[summary.MSE]
        self.assertAlmostEqual(mse['candidate'], .45)
        self.assertNotAlmostEqual(mse['candidate'], np.mean([.5 / 2, 4 / 8]))
        self.assertIsNone(mse['paired_bootstrap'])
        self.assertTrue(mse['finite_cohort_identical'])

    def test_finite_cohort_change_is_not_hidden_by_intersection(self):
        result = summary.compare_metrics(bundle(), bundle(nmae=(.8, np.nan)), ['A', 'B'])
        item = result['de_wilcoxon_lfc_nmae']
        self.assertFalse(item['finite_cohort_identical'])
        self.assertIsNone(item['paired_bootstrap'])
        self.assertEqual(item['candidate_finite_targets'], ['A'])

    def test_incomplete_mse_budget_panel_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'complete frozen panel'):
            summary.compare_metrics(bundle(), bundle(numerator=(1., np.nan)), ['A', 'B'])

    def test_changed_reference_denominator_is_rejected(self):
        with self.assertRaises(AssertionError):
            summary.compare_metrics(bundle(), bundle(denominator=(2.1, 8.)), ['A', 'B'])

    def test_zero_response_alignment_is_not_fabricated(self):
        actual = summary.profile_error(np.zeros(3), np.ones(3), np.zeros(3), np.ones(3, dtype=bool))
        self.assertIsNone(actual['response_residual_cosine'])
        self.assertEqual(actual['response_residual_dot_per_gene'], 0)


if __name__ == '__main__':
    unittest.main()
