"""Count invariants and adversarial feasibility checks for constrained dispersion."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from constrained_dispersion import ConstraintViolation, _integer_box_projection, redistribute_dispersion


def example(seed=2,n=40,g=30,near=False):
    rng=np.random.default_rng(seed)
    p=rng.dirichlet(np.ones(g)*5)
    depth=np.full(n,10000,dtype=np.int64)+(np.arange(n)*2 if near else 0)
    return np.stack([rng.multinomial(int(d),p) for d in depth])


class ConstrainedDispersionTests(unittest.TestCase):
    def check_independent_constraints(self, before, result, rtol=1e-4, atol=1e-8):
        after=result.counts
        self.assertTrue(np.isfinite(after).all())
        self.assertTrue((after>=0).all())
        self.assertTrue(np.equal(after,np.floor(after)).all())
        np.testing.assert_array_equal(after.sum(0,dtype=np.int64),before.sum(0,dtype=np.int64))
        np.testing.assert_array_equal(after.sum(1,dtype=np.int64),before.sum(1,dtype=np.int64))
        depth=before.sum(1,dtype=np.int64)
        a=np.array([sum(float(before[i,j])/float(depth[i]) for i in range(len(before)))*1e6/len(before)
                    for j in range(before.shape[1])])
        b=np.array([sum(float(after[i,j])/float(depth[i]) for i in range(len(after)))*1e6/len(before)
                    for j in range(before.shape[1])])
        tol=np.maximum(atol,rtol*a)
        self.assertTrue((np.abs(b-a)<=tol+1e-9).all())
        np.testing.assert_allclose(result.gene_stats['mean_cpm_change'],b-a,atol=1e-9,rtol=0)
        self.assertTrue(result.diagnostics['constraints_passed'])

    def test_identity_is_bit_exact_and_does_not_mutate(self):
        x=example().astype(np.int32)[:,::2]
        original=x.copy()
        r=redistribute_dispersion(x,1,seed=20261021)
        self.assertEqual(r.counts.dtype,x.dtype)
        self.assertEqual(r.counts.tobytes(),x.tobytes())
        np.testing.assert_array_equal(x,original)
        self.assertEqual(r.diagnostics['status'],'identity')
        self.assertEqual(r.diagnostics['changed_entries'],0)

    def test_equal_depth_reaches_requested_sd_with_exact_cpm(self):
        x=example()
        for factor in [1.5,2.0]:
            r=redistribute_dispersion(x,factor,seed=31)
            self.check_independent_constraints(x,r)
            self.assertTrue(r.diagnostics['target_reached'])
            self.assertAlmostEqual(r.diagnostics['achieved_sd_multiplier'],factor,delta=.01*factor)
            np.testing.assert_array_equal(r.gene_stats['mean_cpm_change'],0)
            self.assertGreater(r.diagnostics['accepted_pairs'],0)

    def test_near_depth_accumulated_budget_multiple_passes(self):
        x=example(7,near=True)
        for seed in [1,7,13]:
            r=redistribute_dispersion(x,2,seed=seed,cpm_rtol=1e-6,max_passes=8)
            self.check_independent_constraints(x,r,rtol=1e-6)
            self.assertGreater(r.diagnostics['accepted_pairs'],0)
            self.assertEqual(r.diagnostics['cpm_constraint_violations'],0)

    def test_no_near_depth_pairs_reports_unattainable(self):
        x=np.array([[50,50],[100,100],[200,100],[300,100]],dtype=np.int32)
        r=redistribute_dispersion(x,2,max_relative_depth_gap=.001)
        np.testing.assert_array_equal(r.counts,x)
        self.assertEqual(r.diagnostics['status'],'unattainable')
        self.assertFalse(r.diagnostics['target_reached'])
        self.check_independent_constraints(x,r)

    def test_zero_and_single_total_genes_do_not_invent_support(self):
        x=example(8,n=20,g=10)
        one=np.zeros((len(x),1),dtype=np.int64);one[0]=1
        x=np.hstack([x,one,np.zeros_like(one)])
        r=redistribute_dispersion(x,2,seed=11)
        self.check_independent_constraints(x,r)
        self.assertEqual(int(r.counts[:,-2].sum()),1)
        self.assertEqual(int(np.count_nonzero(r.counts[:,-2])),1)
        np.testing.assert_array_equal(r.counts[:,-1],0)

    def test_zero_cpm_budget_with_unequal_depths_is_not_relaxed(self):
        x=np.array([[40,10,31,20],[50,10,23,20],[51,15,21,20],[48,15,26,20]])
        r=redistribute_dispersion(x,2,cpm_rtol=0,cpm_atol=0,max_relative_depth_gap=10)
        np.testing.assert_array_equal(r.counts,x)
        self.assertEqual(r.diagnostics['status'],'unattainable')
        self.assertEqual(r.diagnostics['accepted_pairs'],0)
        self.assertEqual(r.diagnostics['max_abs_cpm_change'],0)

    def test_uniform_tied_input_and_one_cell_are_explicit(self):
        for x in [np.ones((16,8),dtype=np.int16)*10,np.array([[0,2,7]],dtype=np.int32)]:
            r=redistribute_dispersion(x,2)
            np.testing.assert_array_equal(r.counts,x)
            self.assertEqual(r.diagnostics['status'],'unattainable')
            self.assertIsNone(r.diagnostics['achieved_sd_multiplier'])
            self.assertFalse(r.diagnostics['target_reached'])
            self.assertTrue(r.diagnostics['constraints_passed'])
            ident=redistribute_dispersion(x,1)
            self.assertEqual(ident.counts.tobytes(),x.tobytes())

    def test_seed_reproducibility(self):
        x=example(8,n=41,near=True)
        a=redistribute_dispersion(x,2,seed=33)
        b=redistribute_dispersion(x,2,seed=33)
        np.testing.assert_array_equal(a.counts,b.counts)

    def test_narrow_input_dtype_cannot_overflow_output(self):
        x=np.array([[250,20],[20,250],[240,30],[30,240]],dtype=np.uint8)
        r=redistribute_dispersion(x,2)
        self.check_independent_constraints(x,r)
        self.assertEqual(r.counts.dtype,np.dtype('int64'))
        self.assertGreater(int(r.counts.max()),255)
        self.assertFalse(r.diagnostics['target_reached'])
        self.assertEqual(r.diagnostics['status'],'partial')
        self.assertTrue(r.diagnostics['duplicate_cells_increased'])

    def test_realized_variances_and_density_are_recomputed(self):
        x=example(16,near=True)
        r=redistribute_dispersion(x,1.5,seed=2)
        cpm=r.counts*1e6/r.counts.sum(1)[:,None]
        actual=cpm.var(0,ddof=1)
        np.testing.assert_allclose(r.gene_stats['variance_cpm_after'],actual,rtol=1e-12,atol=1e-7)
        self.assertEqual(r.diagnostics['nnz_after'],np.count_nonzero(r.counts))
        self.assertEqual(r.diagnostics['duplicate_cells_after']['unique_cells'],len(np.unique(r.counts,axis=0)))

    def test_density_limit_is_a_hard_failure(self):
        with self.assertRaises(ConstraintViolation):redistribute_dispersion(example(),1,max_nnz=0)

    def test_balanced_rounding_adversarial_integer_boundaries(self):
        rng=np.random.default_rng(8)
        for _ in range(50):
            lo=-rng.integers(0,20,31,dtype=np.int64)
            hi=rng.integers(0,20,31,dtype=np.int64)
            desired=rng.integers(-20,20,31)+rng.choice([-1e-12,0,1e-12],31)
            delta=_integer_box_projection(desired,lo,hi,rng)
            self.assertEqual(int(delta.sum()),0)
            self.assertTrue((delta>=lo).all() and (delta<=hi).all())

    def test_invalid_numeric_inputs_are_rejected_before_integer_cast(self):
        for x in [np.array([[1,np.nan]]),np.array([[1,np.inf]]),np.array([[1,-1]]),
                  np.array([[1,.5]]),np.zeros((2,3)),np.array([[2**63]],dtype=np.uint64),
                  np.array([['1','2']]),np.empty((0,3))]:
            with self.assertRaises(ValueError):redistribute_dispersion(x,1.5)
        for factor in [float('nan'),float('inf'),.5]:
            with self.assertRaises(ValueError):redistribute_dispersion(example(),factor)
        with self.assertRaises(ValueError):redistribute_dispersion(example(),1.5,cpm_rtol=-1)


if __name__=='__main__':unittest.main()
