"""Isolation, incumbent equivalence, and actual-emitter audit tests."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import build_bulk_ablation as bulk


def fixtures():
    genes = np.asarray(['A', 'B', 'C', 'D'])
    targets = np.asarray(['A', 'B', 'C', 'low'])
    control = np.asarray([.4, .3, .2, .1], dtype=np.float32)
    probability = np.asarray([[.24,.39,.25,.12], [.44,.18,.23,.15],
                              [.45,.36,.1,.09], [.1,.1,.1,.7]], dtype=np.float32)
    source = bulk.atlas.Source('synthetic', targets, genes, probability, control,
        np.asarray([100, 120, 110, 2]), probability*1e6, control*1e6,
        np.asarray([True, True, True, False]))
    output = np.asarray(['A', 'B', 'missing'])
    mean = np.asarray([.35,.3,.2,.15])
    bctrl = np.asarray([.36,.29,.21,.14])
    # This is already source-centered, as the make_profiles API requires.
    cd4 = (np.asarray([[.05,-.02,.03,.01],[-.03,.04,.02,-.01],[0,0,0,0]], dtype=np.float32),
           np.asarray([[True]*4,[True]*4,[False]*4]))
    return source, output, genes, mean, bctrl, cd4


class BulkAblationTests(unittest.TestCase):
    def test_incumbent_requested_profiles_equal_frozen_builder_sequence(self):
        source, targets, genes, mean, bctrl, cd4 = fixtures()
        effects = bulk.prepare_effects([source], [2.], targets, genes, 'incumbent')
        actual = bulk.make_profiles(effects, mean, bctrl, cd4, .5, 'incumbent')
        expected = []
        # Independent reproduction of the unchanged builder's public calls.
        for space, ctrl, amplitude in [('log2fc', mean, .6), ('bulk_delta', bctrl, .3)]:
            e, coverage = bulk.atlas.fuse_source_centered([source], [2.], targets, genes,
                                                        space=space, common_subtract=1.)
            e, coverage = bulk.atlas.add_cd4_family(e, coverage, *cd4, .5,
                                                   space=space, control_probability=ctrl)
            expected.append(bulk.atlas.desired_mean(ctrl, e, space=space, amplitude=amplitude, clip=3))
        for observed, old in zip(actual, expected):
            np.testing.assert_array_equal(observed, old)

    def test_all_arms_keep_arithmetic_profile_and_cd4_unchanged(self):
        source, targets, genes, mean, bctrl, cd4 = fixtures()
        original_cd4 = tuple(x.copy() for x in cd4)
        profiles = {}
        for arm in bulk.ARMS:
            effects = bulk.prepare_effects([source], [2.], targets, genes, arm)
            profiles[arm] = bulk.make_profiles(effects, mean, bctrl, cd4, .5, arm)
            np.testing.assert_array_equal(profiles[arm][0], profiles['incumbent'][0])
        for observed, old in zip(cd4, original_cd4):
            np.testing.assert_array_equal(observed, old)
        for arm in ('common05', 'bulk015', 'bulk060'):
            self.assertFalse(np.array_equal(profiles[arm][1], profiles['incumbent'][1]))

    def test_centering_population_is_full_source_and_missing_is_masked(self):
        source, targets, genes, mean, bctrl, cd4 = fixtures()
        one = np.asarray(['A'])
        effect, mask = bulk.prepare_effects([source], [1.], one, genes, 'common05')['bulk_delta']
        full, full_mask = bulk.atlas.aligned_effect(source, source.targets, genes,
                                                   space='bulk_delta', common_subtract=.5)
        np.testing.assert_array_equal(effect[0], full[0])
        np.testing.assert_array_equal(mask[0], full_mask[0])
        panel_only, _ = bulk.atlas.aligned_effect(source, one, genes, space='bulk_delta', common_subtract=.5)
        self.assertFalse(np.array_equal(effect, panel_only))
        e, coverage = bulk.prepare_effects([source], [1.], np.asarray(['missing','low']), genes, 'common05')['bulk_delta']
        np.testing.assert_array_equal(e, np.zeros_like(e))
        np.testing.assert_array_equal(coverage, np.zeros_like(coverage))

    def test_common05_changes_only_raw_contribution_before_amplitude(self):
        source, targets, genes, mean, bctrl, cd4 = fixtures()
        fused, raw = {}, {}
        for arm in ('incumbent', 'common05'):
            raw[arm], coverage = bulk.prepare_effects([source], [2.], targets, genes, arm)['bulk_delta']
            fused[arm], den = bulk.atlas.add_cd4_family(raw[arm], coverage, *cd4, .5,
                                                       space='bulk_delta', control_probability=bctrl)
        expected_change = np.divide((raw['common05'] - raw['incumbent'])*coverage, den,
                                    out=np.zeros_like(fused['incumbent']), where=den>0)
        np.testing.assert_allclose(fused['common05']-fused['incumbent'], expected_change,
                                   rtol=1e-6, atol=2e-8)

    def test_amplitude_arms_have_identical_unscaled_effects(self):
        source, targets, genes, *_ = fixtures()
        reference = bulk.prepare_effects([source], [2.], targets, genes, 'incumbent')
        for arm in ('bulk015','bulk060'):
            current = bulk.prepare_effects([source], [2.], targets, genes, arm)
            for space in reference:
                for a, b in zip(current[space], reference[space]):
                    np.testing.assert_array_equal(a, b)

    def test_audited_emission_bitidentical_and_final_diagnostics_independent(self):
        rng = np.random.default_rng(78)
        template = rng.gamma(2., 1., (32, 12))
        template /= template.sum(axis=1, keepdims=True)
        depths = np.linspace(20000, 80000, len(template)).astype(np.int64)
        p = template.mean(axis=0)
        b = (template*depths[:,None]).sum(axis=0)/depths.sum()
        expected = bulk.atlas.dual_moment_counts(template, p, b, depths=depths, seed=20261021)
        actual, audit, stats = bulk.emit_audited(template, depths, p, b, 20261021)
        self.assertEqual(actual.dtype, expected.dtype)
        self.assertEqual(actual.tobytes(), expected.tobytes())
        np.testing.assert_array_equal(actual.sum(axis=1), depths)
        np.testing.assert_array_equal(stats['mean_probability_realized'], (actual/depths[:,None]).mean(axis=0))
        np.testing.assert_array_equal(stats['bulk_probability_realized'], actual.sum(axis=0)/actual.sum())
        self.assertLessEqual(audit['projection_l1'], .03)
        self.assertLessEqual(max(audit['fitting_mean_l1'], audit['fitting_bulk_l1']), .001)
        np.testing.assert_allclose(stats['mean_cpm_rounding_change'],
            1e6*stats['mean_probability_realized']-1e6*stats['mean_probability_fitted'], atol=0, rtol=0)

    def test_equal_depth_projection_is_captured_not_mistaken_for_request(self):
        template = np.tile([.5,.5], (20,1))
        depths = np.full(20, 10000, dtype=np.int64)
        x, audit, stats = bulk.emit_audited(template, depths, np.array([.5,.5]), np.array([.505,.495]), 4)
        np.testing.assert_array_equal(stats['bulk_probability_projected'], [.5,.5])
        self.assertAlmostEqual(audit['projection_l1'], .01)
        self.assertGreater(audit['bulk_probability_l1'], .009)
        self.assertLess(audit['bulk_probability_l1_to_projected'], 1e-12)

    def test_failure_preserves_original_limits_and_restores_profile(self):
        previous = sys.getprofile()
        def sentinel(frame, event, arg):
            pass
        try:
            sys.setprofile(sentinel)
            with self.assertRaisesRegex(ValueError, 'projection too large') as raised:
                bulk.emit_audited(np.tile([.5,.5], (20,1)), np.full(20,10000),
                                  np.array([.5,.5]), np.array([.6,.4]), 4)
            self.assertIs(sys.getprofile(), sentinel)
            self.assertGreater(raised.exception.bulk_diagnostics['projection_error'], .03)
        finally:
            sys.setprofile(previous)

    def test_incumbent_comparison_does_not_hide_realized_cpm_changes(self):
        original = np.asarray([[30,70],[80,120]],dtype=np.int32)
        changed = np.asarray([[31,69],[80,120]],dtype=np.int32)
        depths = original.sum(axis=1)
        stats = {'bulk_probability_realized':changed.sum(axis=0)/changed.sum()}
        audit = bulk.compare_incumbent(changed, original, depths, stats)
        self.assertFalse(audit['counts_bitidentical_to_incumbent'])
        self.assertEqual(audit['genes_exceeding_dispersion_cpm_bound'], 2)
        self.assertFalse(audit['dispersion_cpm_bound_is_enforced'])
        np.testing.assert_allclose(stats['mean_cpm_change_to_incumbent'], [5000,-5000])
        with self.assertRaisesRegex(ValueError, 'bit-identical'):
            bulk.compare_incumbent(changed, original, depths, stats, require_identity=True)

    def test_wrong_rows_or_fractional_counts_rejected(self):
        original = np.array([[30,70],[80,120]])
        for candidate in (original+.1, original+1):
            with self.assertRaises(ValueError):
                bulk.compare_incumbent(candidate, original, original.sum(axis=1), {})


if __name__ == '__main__':
    unittest.main()
