"""Independent scientific invariants for the vendored AtlasShift count emitter.

Run with ``python -m unittest discover -s tests -p test_count_emission.py -v``.
These checks use known, jointly feasible moments rather than matching internal
steps of the emitter. They do not establish biological transfer performance.
"""

from pathlib import Path
import json
import contextlib
import io
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "external" / "atlasshift"))
from model import dual_moment_counts, desired_mean, repair_bulk_totals
from predict import CountWriter, control_template
sys.path.insert(0, str(ROOT / "scripts"))
from build_prediction import Writer
from validate_prediction_receipt import validate, sha256
from prepare_metric_input import prepare as prepare_metric_input


def moment_errors(counts, mean_probability, bulk_probability):
    depth = counts.sum(axis=1, dtype=np.float64)
    mean = (counts / depth[:, None]).mean(axis=0)
    bulk = counts.sum(axis=0, dtype=np.float64) / depth.sum()
    return float(np.abs(mean - mean_probability).sum()), float(np.abs(bulk - bulk_probability).sum())


class CountEmissionTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(42)
        self.template = rng.gamma(1.4, 1, (160, 32))
        self.template[:, -1] = 0  # no donor ever expressed this output gene
        self.template /= self.template.sum(axis=1, keepdims=True)
        self.depths = np.linspace(20000, 80000, len(self.template)).astype(np.int64)
        # A realizable response with new expression and preserved heterogeneity.
        response = self.template * np.exp(np.linspace(-0.35, 0.35, 32))
        response[:, -1] = 0.02
        response /= response.sum(axis=1, keepdims=True)
        self.mean = response.mean(axis=0)
        self.bulk = (response * self.depths[:, None]).sum(axis=0) / self.depths.sum()

    def emit(self, **overrides):
        args = dict(depths=self.depths, seed=123)
        args.update(overrides)
        return dual_moment_counts(self.template, self.mean, self.bulk, **args)

    def test_distinct_mean_and_bulk_moments_and_exact_cell_depths(self):
        self.assertGreater(np.abs(self.mean - self.bulk).sum(), 0.001)
        counts = self.emit()
        self.assertTrue(np.issubdtype(counts.dtype, np.integer))
        self.assertTrue(np.isfinite(counts).all())
        self.assertGreaterEqual(counts.min(), 0)
        np.testing.assert_array_equal(counts.sum(axis=1), self.depths)
        mean_error, bulk_error = moment_errors(counts, self.mean, self.bulk)
        self.assertLess(mean_error, 0.001)
        self.assertLess(bulk_error, 0.001)

    def test_can_activate_gene_missing_from_all_donors(self):
        counts = self.emit()
        self.assertEqual(self.template[:, -1].sum(), 0)
        self.assertTrue((counts[:, -1] > 0).all())
        self.assertLess(abs((counts / self.depths[:, None]).mean(axis=0)[-1] - self.mean[-1]), 0.0001)

    def test_zero_requested_gene_stays_zero(self):
        p = self.mean.copy()
        b = self.bulk.copy()
        p[3] = b[3] = 0
        p /= p.sum()
        b /= b.sum()
        counts = dual_moment_counts(self.template, p, b, depths=self.depths, seed=3)
        self.assertEqual(counts[:, 3].sum(), 0)
        np.testing.assert_array_equal(counts.sum(axis=1), self.depths)

    def test_equal_depths_reduce_to_one_moment(self):
        depths = np.full(len(self.template), 40000)
        counts = dual_moment_counts(self.template, self.mean, self.mean, depths=depths, seed=4)
        mean_error, bulk_error = moment_errors(counts, self.mean, self.mean)
        self.assertLess(mean_error, 0.00001)
        self.assertLess(bulk_error, 0.00001)
        np.testing.assert_array_equal(counts.sum(axis=1), depths)

    def test_seed_reproducibility(self):
        np.testing.assert_array_equal(self.emit(), self.emit())
        self.assertFalse(np.array_equal(self.emit(seed=123), self.emit(seed=456)))

    def test_impossible_bulk_is_rejected(self):
        bulk = np.zeros_like(self.mean)
        bulk[0] = 1
        with self.assertRaisesRegex(ValueError, "projection|support"):
            dual_moment_counts(self.template, self.mean, bulk, depths=self.depths, seed=3)

    def test_invalid_counts_and_probabilities_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "sum to one"):
            dual_moment_counts(self.template, self.mean * 2, self.bulk, depths=self.depths, seed=3)
        with self.assertRaisesRegex(ValueError, "integers"):
            self.emit(depths=self.depths + 0.5)
        bad = self.template.copy()
        bad[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            dual_moment_counts(bad, self.mean, self.bulk, depths=self.depths, seed=3)

    def test_bulk_repair_preserves_depths_and_nearest_column_totals(self):
        expected = np.array([[1.2, 2.8, 4.0], [3.4, 5.6, 1.0], [2.1, 1.1, 5.8]])
        integer = np.array([[2, 2, 4], [4, 5, 1], [2, 2, 5]], dtype=np.int32)
        repaired = repair_bulk_totals(integer.copy(), expected)
        np.testing.assert_array_equal(repaired.sum(axis=1), integer.sum(axis=1))
        self.assertTrue((np.abs(repaired.sum(axis=0) - expected.sum(axis=0)) < 1).all())

    def test_bulk_delta_inverse_can_activate_but_log2fc_cannot(self):
        control = np.array([0.5, 0.5, 0.0])
        effect = np.array([0.0, 0.0, 1.0])
        additive = desired_mean(control, effect, space="bulk_delta", amplitude=1, clip=3)
        multiplicative = desired_mean(control, effect, space="log2fc", amplitude=1, clip=3)
        self.assertGreater(additive[-1], 0)
        self.assertEqual(multiplicative[-1], 0)
        self.assertAlmostEqual(additive.sum(), 1)

    def test_count_writer_roundtrip_and_labels(self):
        counts = self.emit()[:4]
        genes = np.array([f"g{i}" for i in range(counts.shape[1])])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "prediction.h5ad"
            with CountWriter(path, np.array(["T1", "T2"]), genes, ["A"], 2) as writer:
                writer.append(counts[:2])
                writer.append(counts[2:])
            obj = ad.read_h5ad(path)
            np.testing.assert_array_equal(obj.X.toarray(), counts)
            self.assertEqual(list(obj.obs.target_gene.astype(str)), ["T1", "T1", "T2", "T2"])
            self.assertEqual(set(obj.obs.context.astype(str)), {"A"})
            np.testing.assert_array_equal(obj.var_names, genes)
            self.assertEqual(obj.X.nnz, np.count_nonzero(counts))

    def test_control_template_returns_both_original_control_moments(self):
        rng = np.random.default_rng(7)
        raw = rng.poisson(np.linspace(1, 8, 20), (40, 20))
        raw[:20, :5] *= 4
        genes = np.array([f"g{i}" for i in range(raw.shape[1])])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "control.h5ad"
            ad.AnnData(sparse.csr_matrix(raw), var=pd.DataFrame(index=genes)).write_h5ad(path)
            template, depths, mean, bulk = control_template(path, genes, 10, 4, 3)
        np.testing.assert_allclose(template.sum(axis=1), 1)
        np.testing.assert_allclose(mean, (raw / raw.sum(axis=1, keepdims=True)).mean(axis=0))
        np.testing.assert_allclose(bulk, raw.sum(axis=0) / raw.sum())
        self.assertGreater(np.abs(mean - bulk).sum(), 0.01)
        self.assertTrue((depths > 0).all())


class RootWriterReceiptTests(unittest.TestCase):
    """Exercise the wrapper's actual H5AD format with the complete gene axis."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        # Exercise the complete dimension without requiring a private/local CSV.
        self.genes = np.asarray([f'synthetic_gene_{i:05d}' for i in range(18533)])
        self.assertEqual(len(self.genes), 18533)
        self.targets = ['T1', 'T2']
        gene_path, target_path = self.root / 'genes.csv', self.root / 'targets.csv'
        pd.DataFrame({'gene_name': self.genes}).to_csv(gene_path, index=False)
        pd.DataFrame({'target_gene': self.targets}).to_csv(target_path, index=False)
        self.path = self.root / 'prediction.h5ad'
        contexts, moments = [], []
        writer = Writer(self.path, self.genes)
        try:
            for ci, context in enumerate(['A', 'B']):
                x = np.zeros((2, len(self.genes)), dtype=np.int32)
                x[:, ci] = [100, 200]
                control_path = self.root / f'context_{context}.h5ad'
                control = ad.AnnData(sparse.csr_matrix(x),
                    obs=pd.DataFrame({'context': [context]*2, 'target_gene': ['non-targeting']*2}),
                    var=pd.DataFrame(index=self.genes))
                control.write_h5ad(control_path)
                contexts.append(f'{context}={control_path}')
                for target in self.targets:
                    writer.append(x, context, target)
                    moments.append({'context': context, 'target': target, 'nnz': 2, 'max_library': 200})
            writer.close()
        except BaseException:
            writer.f.close()
            raise
        self.receipt_path = self.path.with_suffix('.receipt.json')
        self.receipt = {'arguments': {'genes': str(gene_path), 'targets': str(target_path),
                        'context': contexts, 'cells': 2, 'target_limit': None},
                        'sources': [{'path': '/public/K562_bank.npz'}],
                        'targets': self.targets, 'output_genes': len(self.genes),
                        'moments': moments, 'cells': 8, 'stored_entries': 8,
                        'bytes': self.path.stat().st_size, 'output_sha256': sha256(self.path)}
        self.write_receipt()

    def tearDown(self):
        self.temp.cleanup()

    def write_receipt(self):
        self.receipt_path.write_text(json.dumps(self.receipt))

    def test_full_gene_axis_roundtrip_and_receipt(self):
        result = validate(self.path, self.receipt_path, excluded_source_tokens=['H1'])
        self.assertEqual(result['shape'], [8, 18533])
        self.assertEqual(result['groups'], 4)
        self.assertTrue(result['hash_verified'])

    def test_rejects_wrong_block_width_before_write(self):
        writer = Writer(self.root / 'bad.h5ad', self.genes)
        try:
            with self.assertRaises(ValueError):
                writer.append(np.ones((2, 3)), 'A', 'T1')
        finally:
            writer.f.close()

    def test_rejects_noninteger_counts(self):
        writer = Writer(self.root / 'bad.h5ad', self.genes)
        try:
            with self.assertRaisesRegex(ValueError, 'raw counts'):
                writer.append(np.full((2, len(self.genes)), .5), 'A', 'T1')
        finally:
            writer.f.close()

    def test_rejects_missing_group_diagnostics(self):
        self.receipt['moments'].pop()
        self.write_receipt()
        with self.assertRaisesRegex(ValueError, 'moment groups'):
            validate(self.path, self.receipt_path)

    def test_rejects_wrong_gene_order(self):
        pd.DataFrame({'gene_name': self.genes[::-1]}).to_csv(self.receipt['arguments']['genes'], index=False)
        with self.assertRaisesRegex(ValueError, 'gene order|gene axis'):
            validate(self.path, self.receipt_path)

    def test_rejects_wrong_per_group_cell_count(self):
        self.receipt['arguments']['cells'] = 3
        self.write_receipt()
        with self.assertRaisesRegex(ValueError, 'shape|cell counts'):
            validate(self.path, self.receipt_path)

    def test_rejects_receipt_hash_mismatch(self):
        self.receipt['output_sha256'] = '0'*64
        self.write_receipt()
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            validate(self.path, self.receipt_path)

    def test_rejects_destination_response_source(self):
        self.receipt['sources'] = [{'path': '/public/H1_2025_full_statistics.npz'}]
        self.write_receipt()
        with self.assertRaisesRegex(ValueError, 'Forbidden evaluation source'):
            validate(self.path, self.receipt_path, excluded_source_tokens=['H1'])

    def test_rejects_swapped_control_context(self):
        self.receipt['arguments']['context'][0] = self.receipt['arguments']['context'][0].replace('A=', 'C=')
        self.write_receipt()
        with self.assertRaisesRegex(ValueError, 'Control context'):
            validate(self.path, self.receipt_path)

    def test_local_scoring_copy_adds_exact_controls_without_changing_prediction(self):
        source_digest = sha256(self.path)
        output = self.root / 'local_scoring.h5ad'
        control_path = self.root / 'context_A.h5ad'
        with contextlib.redirect_stdout(io.StringIO()):
            prepare_metric_input(self.path, control_path, output)
        self.assertEqual(source_digest, sha256(self.path))
        source, local, controls = ad.read_h5ad(self.path), ad.read_h5ad(output), ad.read_h5ad(control_path)
        self.assertEqual(local.shape, (10, 18533))
        np.testing.assert_array_equal(local.X[:8].toarray(), source.X.toarray())
        np.testing.assert_array_equal(local.X[8:].toarray(), controls.X.toarray())
        self.assertEqual(int((local.obs.target_gene.astype(str) == 'non-targeting').sum()), 2)

    def test_local_scoring_copy_refuses_reference_without_controls(self):
        with self.assertRaisesRegex(ValueError, 'Reference has no controls'):
            prepare_metric_input(self.path, self.path, self.root / 'invalid_scoring.h5ad')


if __name__ == "__main__":
    unittest.main()
