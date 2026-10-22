"""Evidence-voting invariants for the experimental response confidence gate."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from response_confidence import confidence_gate


class ConfidenceGateTests(unittest.TestCase):
    def test_two_source_rule_handles_agreement_conflict_weak_and_zero(self):
        values = np.array([[[.1, .1, .1, .02, -.3, .2]],
                           [[.2, -.1, .01, .02, -.1, .2]]])
        fused = np.array([[.15, .01, .08, .02, -.2, 0.]])
        gate, detail = confidence_gate(values, np.ones_like(values, dtype=bool), fused)
        np.testing.assert_array_equal(gate, [[True, False, False, True, True, False]])
        self.assertEqual(detail['closed_pairs'], 3)
        self.assertEqual(detail['eligible_with_zero_fused_effect'], 1)
        json.dumps(detail, allow_nan=False)

    def test_missing_data_never_votes_even_with_arbitrary_or_nonfinite_values(self):
        values = np.array([[[.1, .1, np.nan]], [[.2, -1e9, np.inf]], [[np.nan, np.inf, -.1]]])
        mask = np.array([[[True, True, False]], [[True, False, False]], [[False, False, False]]])
        gate, detail = confidence_gate(values, mask, np.array([[.15, .1, 0.]]))
        np.testing.assert_array_equal(gate, [[True, True, True]])
        self.assertEqual(detail['eligible_pairs'], 1)
        self.assertEqual(detail['single_source_pairs'], 1)
        self.assertEqual(detail['unmeasured_pairs'], 1)

    def test_weak_measured_sources_remain_in_agreement_denominator(self):
        values = np.array([[[.2]], [[.1]], [[.001]]])
        measured = np.ones_like(values, dtype=bool)
        gate, _ = confidence_gate(values, measured, np.array([[.2]]))
        self.assertFalse(gate.item())  # 2/3 is below .75, not 2/2.
        measured[2] = False
        gate, _ = confidence_gate(values, measured, np.array([[.2]]))
        self.assertTrue(gate.item())

    def test_three_of_four_votes_pass_at_inclusive_threshold(self):
        values = np.array([[[.1]], [[.2]], [[.3]], [[-.1]]])
        measured = np.ones_like(values, dtype=bool)
        gate, _ = confidence_gate(values, measured, np.array([[.2]]), min_agreement=.75)
        self.assertTrue(gate.item())
        gate, _ = confidence_gate(values, measured, np.array([[.2]]), min_agreement=.751)
        self.assertFalse(gate.item())

    def test_supplied_fused_direction_overrides_unweighted_source_majority(self):
        values = np.array([[[.1]], [[.1]], [[.1]], [[-2.]]])
        mask = np.ones_like(values, dtype=bool)
        gate, _ = confidence_gate(values, mask, np.array([[-.5]]))
        self.assertFalse(gate.item())  # Three positive votes cannot support negative fusion.
        gate, _ = confidence_gate(values, mask, np.array([[.1]]))
        self.assertTrue(gate.item())

    def test_source_permutation_invariance_and_no_input_mutation(self):
        rng = np.random.default_rng(5)
        effects = rng.normal(size=(4, 3, 9))
        measured = rng.random(effects.shape) > .2
        fused = rng.normal(size=(3, 9))
        originals = [x.copy() for x in (effects, measured, fused)]
        gate, detail = confidence_gate(effects, measured, fused)
        order = [3, 1, 0, 2]
        other, other_detail = confidence_gate(effects[order], measured[order], fused)
        np.testing.assert_array_equal(gate, other)
        self.assertEqual(detail, other_detail)
        for original, current in zip(originals, (effects, measured, fused)):
            np.testing.assert_array_equal(original, current)

    def test_zero_sources_leave_gate_open(self):
        gate, detail = confidence_gate(np.empty((0, 2, 3)), np.empty((0, 2, 3), dtype=bool), np.zeros((2, 3)))
        self.assertTrue(gate.all())
        self.assertIsNone(detail['retained_fraction_among_eligible'])

    def test_zero_threshold_does_not_turn_zero_effects_into_directional_votes(self):
        effects = np.zeros((2, 1, 1))
        gate, _ = confidence_gate(effects, np.ones_like(effects, dtype=bool), np.ones((1, 1)), min_abs_effect=0)
        self.assertFalse(gate.item())

    def test_at_least_two_votes_are_required_even_when_fraction_threshold_is_zero(self):
        effects = np.array([[[.1]], [[0.]]])
        gate, _ = confidence_gate(effects, np.ones_like(effects, dtype=bool), np.ones((1, 1)), min_agreement=0)
        self.assertFalse(gate.item())

    def test_nonfinite_measured_or_fused_values_reject(self):
        effects = np.array([[[.1]], [[np.nan]]])
        with self.assertRaisesRegex(ValueError, 'finite'):
            confidence_gate(effects, np.ones_like(effects, dtype=bool), np.ones((1, 1)))
        with self.assertRaisesRegex(ValueError, 'finite'):
            confidence_gate(np.ones((2, 1, 1)), np.ones((2, 1, 1), dtype=bool), np.array([[np.inf]]))

    def test_invalid_axes_masks_and_thresholds_reject(self):
        with self.assertRaisesRegex(ValueError, 'shape'):
            confidence_gate(np.ones((2, 1)), np.ones((2, 1), dtype=bool), np.ones(1))
        with self.assertRaisesRegex(ValueError, 'boolean'):
            confidence_gate(np.ones((2, 1, 1)), np.ones((2, 1, 1)), np.ones((1, 1)))
        for threshold in [-.1, np.nan]:
            with self.assertRaisesRegex(ValueError, 'min_abs_effect'):
                confidence_gate(np.ones((2, 1, 1)), np.ones((2, 1, 1), dtype=bool), np.ones((1, 1)), min_abs_effect=threshold)
        for threshold in [-.1, 1.1, np.nan]:
            with self.assertRaisesRegex(ValueError, 'min_agreement'):
                confidence_gate(np.ones((2, 1, 1)), np.ones((2, 1, 1), dtype=bool), np.ones((1, 1)), min_agreement=threshold)


if __name__ == '__main__':
    unittest.main()
