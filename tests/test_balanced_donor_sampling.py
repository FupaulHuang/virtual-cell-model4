"""Martingale enumeration, inclusion frequencies, cardinality and frozen features."""
from pathlib import Path
import json
import sys
import unittest
import numpy as np
from scipy import sparse

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import balanced_donor_sampling as sampling


def enumerate_design(features,cells):
    """Enumerate both randomized branches for a fixed unit order, without RNG."""
    prepared=sampling.prepare_sampler(features,cells);leaves={}
    def recurse(p,mass):
        active=np.flatnonzero((p>0)&(p<1))
        if len(active)>prepared.rank:
            indices=active[:prepared.rank+1]
            a,b,q,_=sampling.flight_branches(p[indices],prepared.constraint[:,indices])
        elif len(active)>1:
            indices=active[:2];a,b,q,_=sampling.pivotal_branches(*p[indices])
        else:
            if len(active):
                remaining=cells-int(np.count_nonzero(p==1))
                if remaining not in (0,1) or abs(p[active[0]]-remaining)>sampling.TERMINAL_SNAP_TOL:
                    raise AssertionError('Enumeration left a genuinely fractional sample')
                p=p.copy();p[active[0]]=remaining
            key=tuple(np.flatnonzero(p==1))
            if len(key)!=cells:raise AssertionError('Enumeration cardinality violation')
            leaves[key]=leaves.get(key,0)+mass;return
        first=p.copy();second=p.copy();first[indices]=a;second[indices]=b
        if q>0:recurse(first,mass*q)
        if q<1:recurse(second,mass*(1-q))
    recurse(np.full(len(features),cells/len(features)),1.)
    return leaves


class BalancedSamplingTests(unittest.TestCase):
    def test_flight_branches_preserve_conditional_marginals_and_balance(self):
        a=np.array([[1.,1.,1.,1.],[-2.,-.5,.3,1.7],[1.,-1.,2.,0.]])
        p=np.array([.1,.7,.3,.8])
        plus,minus,q,info=sampling.flight_branches(p,a)
        np.testing.assert_allclose(q*plus+(1-q)*minus,p,atol=2e-12,rtol=0)
        np.testing.assert_allclose(a@plus,a@p,atol=2e-12,rtol=0)
        np.testing.assert_allclose(a@minus,a@p,atol=2e-12,rtol=0)
        self.assertTrue(np.any((plus==0)|(plus==1)));self.assertTrue(np.any((minus==0)|(minus==1)))
        self.assertLessEqual(info['null_residual'],sampling.NULL_TOL)

    def test_pivotal_branches_near_boundaries_preserve_expected_values(self):
        for a,b in [(0.1,.2),(.3,.7),(.8,.7),(1e-13,.4),(1-1e-13,1-2e-13)]:
            first,second,q,_=sampling.pivotal_branches(a,b)
            np.testing.assert_allclose(q*first+(1-q)*second,[a,b],atol=2e-12,rtol=0)
            self.assertAlmostEqual(first.sum(),a+b,places=11)
            self.assertAlmostEqual(second.sum(),a+b,places=11)
            self.assertTrue(np.all((first>=0)&(first<=1)));self.assertTrue(np.all((second>=0)&(second<=1)))

    def test_exhaustive_tiny_designs_have_exact_first_order_inclusion(self):
        fixtures=[(np.ones((6,1)),2),
                  (np.column_stack([np.ones(6),[-2,-1,0,0,1,2]]),3),
                  (np.column_stack([np.ones(8),np.arange(8),np.sin(np.arange(8))]),3),
                  (np.column_stack([np.ones(5),np.eye(5)]),2)]
        for features,cells in fixtures:
            outcomes=enumerate_design(features,cells)
            self.assertAlmostEqual(sum(outcomes.values()),1.,places=11)
            inclusion=np.zeros(len(features))
            for indices,probability in outcomes.items():inclusion[list(indices)]+=probability
            np.testing.assert_allclose(inclusion,cells/len(features),rtol=0,atol=2e-11)

    def test_seeded_monte_carlo_inclusion_and_fixed_size(self):
        n=10;features=np.column_stack([np.ones(n),np.linspace(-1,1,n),np.cos(np.arange(n))])
        prepared=sampling.prepare_sampler(features,3);frequency=np.zeros(n);draws=6000
        for seed in range(draws):
            indices,audit=prepared.select(seed);frequency[indices]+=1
            self.assertEqual(len(indices),3);self.assertEqual(len(set(indices)),3)
            self.assertLessEqual(audit['accumulated_snap_l1'],sampling.MAX_SNAP_L1)
        np.testing.assert_allclose(frequency/draws,.3,atol=.035,rtol=0)

    def test_population4000_sample400_reusable_sampler(self):
        rng=np.random.default_rng(81);features=np.column_stack([np.ones(4000),rng.normal(size=(4000,5))])
        prepared=sampling.prepare_sampler(features,400);original=features.copy()
        first,audit=prepared.select(20261026)
        np.testing.assert_array_equal(first,prepared.select(20261026)[0])
        for seed in range(4):
            indices,info=prepared.select(seed)
            self.assertEqual(len(indices),400);self.assertEqual(info['inclusion_probability'],.1)
            self.assertLessEqual(info['fractional_units_before_landing'],info['rank'])
            self.assertLessEqual(max(np.abs(info['flight_retained_constraint_residual'])),sampling.MAX_SNAP_L1)
            self.assertLessEqual(info['accumulated_snap_l1'],sampling.MAX_SNAP_L1)
        np.testing.assert_array_equal(features,original)
        self.assertFalse(prepared.features.flags.writeable)

    def test_rank_deficient_columns_and_approximate_landing_honestly_reported(self):
        x=np.arange(7,dtype=float)
        features=np.column_stack([np.ones(7),x,2*x,np.ones(7),np.zeros(7)])
        prepared=sampling.prepare_sampler(features,2)
        self.assertEqual(prepared.rank,2)
        indices,audit=prepared.select(12)
        np.testing.assert_allclose(audit['final_feature_total_residual'],features[indices].sum(0)-2*features.mean(0),atol=1e-12)
        full_rank=sampling.prepare_sampler(np.column_stack([np.ones(5),np.eye(5)]),2)
        _,record=full_rank.select(3)
        self.assertEqual(record['flight_steps'],0)
        self.assertGreater(np.linalg.norm(record['final_feature_total_residual']),0)
        self.assertIn('Approximate',record['auxiliary_balance_claim'])

    def test_constant_control_counts_produce_padded_zero_features(self):
        raw=np.tile([10,20,30],(20,1))
        features,state=sampling.fit_features(raw)
        np.testing.assert_array_equal(features[:,0],1)
        np.testing.assert_array_equal(features[:,1:],0)
        self.assertEqual(state['pca_components_retained'],0)
        self.assertEqual(sampling.prepare_sampler(features,2).rank,1)
        json.dumps(state,allow_nan=False)

    def test_control_features_reconstruct_from_frozen_state_and_raw_counts(self):
        rng=np.random.default_rng(33);raw=rng.poisson(rng.uniform(3,30,(96,80))).astype(np.int32)
        original=raw.copy();features,state=sampling.fit_features(raw,seed=20261026)
        self.assertEqual(features.shape,(96,6))
        self.assertEqual(state['pca_method'],'seeded_numerical_truncated_arpack_svd')
        indices=np.array(state['selected_gene_indices']);depth=raw.sum(1)
        cpm=1e6*raw/depth[:,None]
        z=(cpm[:,indices]-np.array(state['selected_gene_mean_cpm']))/np.array(state['selected_gene_std_cpm'])
        scores=z@np.array(state['loadings'])
        expected=(scores-np.array(state['pc_score_means']))/np.array(state['pc_score_std'])
        np.testing.assert_allclose(features[:,2:],expected,atol=1e-11,rtol=1e-11)
        self.assertLess(max(state['singular_relative_residual_right']),1e-6)
        self.assertLess(max(state['singular_relative_residual_left']),1e-6)
        np.testing.assert_array_equal(raw,original)
        repeated,other=sampling.fit_features(sparse.csr_matrix(raw),seed=20261026)
        self.assertEqual(state['input_counts_sha256'],other['input_counts_sha256'])
        np.testing.assert_allclose(features,repeated,atol=1e-10,rtol=1e-10)
        self.assertEqual(state['feature_matrix_sha256'],sampling.array_sha256(features))
        json.dumps(state,allow_nan=False)

    def test_sampling_returns_unmodified_actual_rows_and_empty_full_cases(self):
        raw=np.arange(60).reshape(10,6)+1;features,_=sampling.fit_features(raw)
        before=raw.copy();rows,audit=sampling.select_balanced(features,3,cells=4)
        self.assertEqual(raw[rows].dtype,raw.dtype);np.testing.assert_array_equal(raw,before)
        self.assertFalse(audit['raw_cells_modified'])
        self.assertEqual(len(sampling.select_balanced(features,3,cells=0)[0]),0)
        np.testing.assert_array_equal(sampling.select_balanced(features,3,cells=10)[0],np.arange(10))

    def test_invalid_counts_features_and_sizes_rejected(self):
        for raw in [np.zeros((2,3)),np.array([[1.,-.1]]),np.array([[1.2,3.]]),np.array([[np.nan,2.]])]:
            with self.assertRaises(ValueError):sampling.fit_features(raw)
        with self.assertRaises(ValueError):sampling.prepare_sampler(np.zeros((10,2)),3)
        with self.assertRaises(ValueError):sampling.prepare_sampler(np.ones((10,2)),11)
        with self.assertRaises(ValueError):sampling.prepare_sampler(np.ones((10,2)),True)


if __name__=='__main__':unittest.main()
