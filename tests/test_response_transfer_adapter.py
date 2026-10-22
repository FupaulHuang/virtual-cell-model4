"""Leakage, source-family, numerical-unit and inference-contract tests."""
from pathlib import Path
import copy
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import response_transfer_adapter as adapter


def bank(family, factor=1.0):
    targets = np.asarray(['T1', 'T2', 'T3', 'T4'])
    genes = np.asarray(['T1', 'T2', 'T3', 'T4', 'A', 'B', 'C', 'D'])
    ctrl = np.asarray([.05, .07, .09, .11, .13, .15, .18, .22])
    patterns = np.asarray([[.05, -.04, .03, -.02, .03, -.02, .01, -.04],
                           [-.03, .05, -.02, .03, -.04, .02, -.01, 0],
                           [.02, -.03, .04, -.01, .01, -.02, .02, -.03],
                           [-.02, .01, -.02, .04, -.03, .03, -.02, .01]])
    probability = ctrl[None, :] + factor * patterns
    probability /= probability.sum(axis=1, keepdims=True)
    counts = np.rint(probability * 100000).astype(np.int64)
    means = probability * 1e6
    return adapter.RawBank('synthetic_' + family, targets, genes, np.array([100, 110, 120, 130]),
                           counts, means, np.broadcast_to(ctrl, means.shape),
                           np.broadcast_to(ctrl * 1e6, means.shape), np.ones(len(genes), dtype=bool))


def banks():
    return [bank(f, factor) for f, factor in zip(adapter.FAMILIES, (1., .8, 1.2))]


class ResponseAdapterTests(unittest.TestCase):
    def test_explicit_arithmetic_and_bulk_units_use_matched_controls(self):
        b = bank('K562')
        actual = b.effects('log2fc', 0)[0]
        expected = np.log2(b.mean_cpm + 1) - np.log2(b.control_mean_cpm + 1)
        np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-7)
        actual_bulk = b.effects('bulk_delta', 0)[0]
        p = b.counts / b.counts.sum(axis=1, keepdims=True)
        np.testing.assert_allclose(actual_bulk, np.log1p(50000*p) - np.log1p(50000*b.control_probability),
                                   rtol=1e-6, atol=1e-7)

    def test_full_source_common_excludes_own_coordinate(self):
        b = bank('K562')
        effect, common, _ = b.effects('log2fc', 0)
        for j, gene in enumerate(b.genes):
            valid = b.targets != gene
            self.assertAlmostEqual(common[j], float(effect[valid, j].mean()), places=6)
        pair = b.pair_data(['T1'], ['A'], 'log2fc', 0, 20)
        self.assertEqual(pair['common'][0], common[b._cols['A']])

    def test_h1_unknown_and_duplicate_source_representations_rejected(self):
        for text in ('H1_2025_full', 'H1_K562', 'CD4', 'HCT116_HEK293T'):
            with self.assertRaises(ValueError):
                adapter._family(text)
        with self.assertRaisesRegex(ValueError, 'H1'):
            adapter._family('K562', '/tmp/H1_bank.npz')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            adapter.validate_families([bank('K562'), bank('K562')])

    def test_outer_heldout_changes_cannot_change_any_training_value(self):
        original = banks()
        first = adapter.fold_training_examples(original, 'K562', 'log2fc', max_pairs=24, prior_counts=0)
        poisoned = banks()
        # Change all held-out responses, controls and source-wide shared components.
        poisoned[0].mean_cpm *= 7
        poisoned[0].counts *= 3
        poisoned[0].control_mean_cpm = poisoned[0].control_mean_cpm * 4
        second = adapter.fold_training_examples(poisoned, 'K562', 'log2fc', max_pairs=24, prior_counts=0)
        for a, b in zip(first, second):
            self.assertNotEqual(a['destination'], 'K562')
            self.assertNotIn('K562', a['input_families'])
            for key in ('x', 'confidence_x', 'y', 'weight', 'confidence_y', 'targets', 'genes'):
                np.testing.assert_array_equal(a[key], b[key])
        self.assertEqual(adapter.fit_ridge(first, .03), adapter.fit_ridge(second, .03))

    def test_destination_response_does_not_enter_its_features_or_sampling(self):
        a, b = bank('K562'), bank('HCT116')
        first = adapter.examples(a, [b], 'bulk_delta', max_pairs=24, prior_counts=0)
        a.counts = a.counts[:, ::-1].copy(); a.mean_cpm *= 9; a._cache.clear()
        second = adapter.examples(a, [b], 'bulk_delta', max_pairs=24, prior_counts=0)
        for key in ('targets', 'genes', 'x', 'confidence_x', 'weight'):
            np.testing.assert_array_equal(first[key], second[key])
        self.assertFalse(np.array_equal(first['y'], second['y']))
        with self.assertRaisesRegex(ValueError, 'own input'):
            adapter.examples(a, [a], 'log2fc', max_pairs=24)

    def test_samples_are_response_blind_own_excluded_and_target_balanced(self):
        a, b = bank('K562'), bank('HCT116')
        task = adapter.examples(a, [b], 'log2fc', max_pairs=23, prior_counts=0)
        self.assertTrue(np.all(task['targets'] != task['genes']))
        self.assertEqual(len(set(zip(task['targets'], task['genes']))), len(task['y']))
        for target in a.targets:
            self.assertAlmostEqual(task['weight'][task['targets'] == target].sum(), .25)
        again = adapter.examples(a, [b], 'log2fc', max_pairs=23, prior_counts=0)
        np.testing.assert_array_equal(task['genes'], again['genes'])

    def test_missing_pairs_and_unmeasured_genes_are_not_zero_response_votes(self):
        a, b = bank('K562'), bank('HCT116')
        b.measured[b._cols['A']] = False
        x, _, n, mean = adapter.build_features([a, b], ['T1', 'missing', 'T1'], ['A', 'A', 'absent'],
                                               [100., 100., 100.], 'log2fc', 0)
        np.testing.assert_array_equal(n, [1, 0, 0])
        self.assertAlmostEqual(mean[0], a.effects('log2fc', 0)[0][0, a._cols['A']])
        np.testing.assert_array_equal(x[1:], 0)

    def test_ridge_zero_response_has_no_intercept_background(self):
        task = adapter.examples(bank('K562'), [bank('HCT116')], 'log2fc', max_pairs=24)
        head = adapter.fit_ridge([task], .03)
        np.testing.assert_array_equal(adapter.apply_head(head, np.zeros((3, len(adapter.FEATURE_NAMES)))), 0)
        with self.assertRaises(ValueError):
            adapter.fit_ridge([task], 0)

    def test_uncertainty_mask_distinguishes_observed_zero_from_missing(self):
        a = bank('K562')
        a.mean_variance = np.zeros_like(a.mean_cpm); a.control_mean_variance = np.zeros_like(a.mean_cpm)
        with_var = a.pair_data(['T1'], ['A'], 'log2fc', 0, 20)
        missing = bank('HCT116').pair_data(['T1'], ['A'], 'log2fc', 0, 20)
        self.assertTrue(with_var['se_known'][0]); self.assertFalse(missing['se_known'][0])
        self.assertEqual(with_var['arithmetic_se'][0], missing['arithmetic_se'][0])

    def test_loso_report_and_inference_preserve_axes_and_serialization(self):
        sources = banks()
        model, report = adapter.fit_loso(sources, lambdas=(.03, .3), max_pairs=24, prior_counts=0)
        model = json.loads(json.dumps(model))
        for space in adapter.SPACES:
            for rows in report['loso'][space]['grid'].values():
                self.assertEqual(len(rows), 3)
                for fold in rows:
                    self.assertNotIn(fold['heldout_context'], fold['training_label_contexts'])
                    self.assertNotIn(fold['heldout_context'], fold['training_input_contexts'])
                    self.assertNotIn(fold['heldout_context'], fold['test_input_contexts'])
        targets, genes = np.asarray(['T4', 'T1', 'unknown']), np.asarray(['D', 'A', 'missing'])
        result, confidence = adapter.predict_effects(model, sources, targets, genes, np.array([1e5, 2e5, 0]), chunk_targets=1)
        other, _ = adapter.predict_effects(model, sources, targets[::-1], genes[::-1], np.array([0, 2e5, 1e5]), chunk_targets=3)
        for space in adapter.SPACES:
            self.assertEqual(result[space][0].shape, (3, 3))
            np.testing.assert_allclose(result[space][0], other[space][0][::-1, ::-1], rtol=1e-6, atol=1e-7)
            np.testing.assert_array_equal(result[space][0][-1], 0)
            np.testing.assert_array_equal(result[space][0][:, -1], 0)
            self.assertTrue(((confidence[space] >= 0) & (confidence[space] <= 1)).all())
        changed = copy.deepcopy(model); changed['source_banks'][0]['sha256'] = 'different'
        with self.assertRaisesRegex(ValueError, 'identity'):
            adapter.predict_effects(changed, sources, targets, genes, [1e5, 2e5, 0])

    def test_exact_unit_inverse_and_zero_effect_identity(self):
        control = np.array([0., 1., 9.]); bulk = np.array([.1, .3, .6])
        wanted = np.array([1., 3., 7.])
        bulk_wanted = np.array([.15, .25, .6])
        effects = {'log2fc': (np.log2((wanted+1)/(control+1))[None], np.ones((1,3))),
                   'bulk_delta': ((np.log1p(50000*bulk_wanted)-np.log1p(50000*bulk))[None], np.ones((1,3)))}
        p, b = adapter.effects_to_profiles(effects, control, bulk)
        np.testing.assert_allclose(p[0], wanted/wanted.sum(), rtol=1e-12)
        np.testing.assert_allclose(b[0], bulk_wanted, rtol=1e-12)
        zero = {s: (np.zeros((2,3)), np.ones((2,3))) for s in adapter.SPACES}
        p, b = adapter.effects_to_profiles(zero, control, bulk)
        np.testing.assert_allclose(p, np.broadcast_to(control/control.sum(), (2,3)), atol=1e-15)
        np.testing.assert_allclose(b, np.broadcast_to(bulk, (2,3)), atol=1e-15)


if __name__ == '__main__':
    unittest.main()
