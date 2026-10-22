import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from prepare_emission_null import split_quotas, round_template, matched_marginal_counts


class NullEngineeringTests(unittest.TestCase):
    def test_matching_preserves_joint_counts_when_feasible(self):
        np.testing.assert_array_equal(matched_marginal_counts(
            ['a|g','a|h','b|g','b|h'],[4,1,1,4],[9,9,9,9]),[4,1,1,4])

    def test_capacity_repair_keeps_both_margins(self):
        np.testing.assert_array_equal(matched_marginal_counts(
            ['a|g','a|h','b|g','b|h'],[4,1,1,4],[2,9,9,9]),[2,3,3,2])

    def test_infeasible_matching_is_rejected(self):
        with self.assertRaises(ValueError):
            matched_marginal_counts(['a|g','a|h','b|g','b|h'],[4,1,1,4],[0,0,9,9])

    def test_balanced_strata(self):
        counts = np.random.default_rng(77).multinomial(4000, np.ones(71) / 71)
        q = split_quotas(counts)
        np.testing.assert_array_equal(q.sum(0), counts)
        np.testing.assert_array_equal(q.sum(1), np.full(10, 400))
        self.assertTrue(np.all(q.max(0) - q.min(0) <= 1))
        np.testing.assert_array_equal(q, split_quotas(counts))

    def test_many_rare_strata(self):
        q = split_quotas(np.ones(4000, dtype=int))
        self.assertTrue(np.all(q.sum(0) == 1))
        self.assertTrue(np.all(q.sum(1) == 400))

    def test_integer_template_identity(self):
        x = np.array([[4, 3, 2, 0], [10, 20, 2, 0], [1, 7, 11, 0]], dtype=np.int32)
        depths = x.sum(1)
        y = round_template(x / depths[:, None], depths, 42)
        np.testing.assert_array_equal(x, y)

    def test_rounding_preserves_depths_and_nearest_column_totals(self):
        rng = np.random.default_rng(9)
        p = rng.dirichlet(np.ones(19), size=31)
        p = np.c_[p, np.zeros(31)]
        depths = rng.integers(400, 900, 31)
        x = round_template(p, depths, 123)
        np.testing.assert_array_equal(x.sum(1), depths)
        self.assertTrue(np.all(x >= 0))
        self.assertTrue(np.all(x[:, -1] == 0))
        expected = (p * depths[:, None]).sum(0)
        self.assertTrue(np.all(np.abs(x.sum(0) - expected) < 1 + 1e-10))
        np.testing.assert_array_equal(x, round_template(p, depths, 123))


if __name__ == '__main__':
    unittest.main()
