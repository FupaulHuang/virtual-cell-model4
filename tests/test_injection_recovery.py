"""Independent finite-population and count-law checks for oracle references."""
from pathlib import Path
import ast
import itertools
import sys
import tempfile
import unittest
import numpy as np
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import prepare_injection_recovery as truth


class InjectionRecoveryTests(unittest.TestCase):
    def test_object_string_gene_input_has_pickle_free_roundtrip(self):
        m=dict(p=np.array([.4,.6]),b=np.array([.3,.7]),mean_counts=np.array([3.,7.]),mean_depth=10.)
        bank=truth.profile_bank(np.array(['NULL'],dtype=object),[1],['NULL'],[0.],m,[m],
                                [[False,False]],[[True,True]],np.array(['g1','g2'],dtype=object))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'profiles.npz';np.savez_compressed(path,**bank)
            with np.load(path,allow_pickle=False) as saved:
                for name in saved.files:
                    self.assertNotEqual(saved[name].dtype.kind,'O')
                    np.testing.assert_array_equal(saved[name],bank[name])

    def test_fisher_allocation_matches_exhaustive_weighted_subsets(self):
        for tau in (0.,.7,-1.2):
            subsets=list(itertools.combinations(range(7),3))
            ks=np.asarray([sum(i>=4 for i in s) for s in subsets])
            weights=np.exp(tau*ks);weights/=weights.sum()
            k,p=truth.fisher_allocation(4,3,3,tau)
            expected=np.asarray([weights[ks==v].sum() for v in k])
            np.testing.assert_allclose(p,expected,rtol=1e-14,atol=1e-15)
        k,p=truth.fisher_allocation(4,3,3,0.)
        self.assertAlmostEqual(float(k@p),3*3/7,places=13)

    def test_composition_draw_is_unique_reproducible_and_supported(self):
        high=np.array([False]*8+[True]*5)
        for tau in (0.,2.,-2.):
            x=truth.draw_composition(high,7,tau,33)
            np.testing.assert_array_equal(x,truth.draw_composition(high,7,tau,33))
            self.assertEqual(len(set(x)),7)
            self.assertTrue(np.all((x>=0)&(x<len(high))))
        np.testing.assert_array_equal(np.sort(truth.draw_composition(high,13,1.,3)),np.arange(13))

    def test_mixture_moments_match_exhaustive_raw_subset_expectations(self):
        x=np.array([[1,2,0],[3,0,1],[1,1,4],[0,2,7],[6,1,1]],dtype=np.int32)
        high=np.array([False,False,True,True,True]);n=2;tau=.6
        subsets=list(itertools.combinations(range(len(x)),n))
        w=np.exp(tau*np.array([high[list(s)].sum() for s in subsets]));w/=w.sum()
        expected_p=sum(weight*(x[list(s)]/x[list(s)].sum(1)[:,None]).mean(0) for s,weight in zip(subsets,w))
        expected_c=sum(weight*x[list(s)].mean(0) for s,weight in zip(subsets,w))
        m,_=truth.composition_moments(sparse.csr_matrix(x),high,n,tau)
        np.testing.assert_allclose(m['p'],expected_p,rtol=1e-14)
        np.testing.assert_allclose(m['mean_counts'],expected_c,rtol=1e-14)
        np.testing.assert_allclose(m['b'],expected_c/expected_c.sum(),rtol=1e-14)

    def test_null_uses_own_population_and_bulk_is_not_expected_ratio(self):
        x=sparse.csr_matrix([[1,0],[0,9]],dtype=np.int32)
        m,_=truth.composition_moments(x,np.array([False,True]),1,0.)
        np.testing.assert_allclose(m['p'],[.5,.5])
        np.testing.assert_allclose(m['b'],[.1,.9])
        # Mean of the one-cell observed bulk ratios is(.5,.5), not(.1,.9).
        self.assertFalse(np.allclose(m['p'],m['b']))
        bank=truth.profile_bank(['NULL'],[1],['NULL'],[0.],m,[m],[[False,False]],[[True,True]],['g1','g2'])
        np.testing.assert_array_equal(bank['delta_p'],np.zeros((1,2)))
        np.testing.assert_array_equal(bank['delta_b'],np.zeros((1,2)))
        np.testing.assert_array_equal(bank['truth_lfc_cpm_eps1e9'],np.zeros((1,2)))

    def test_joint_depth_state_moments_are_analytical_and_distinct(self):
        mu=np.array([[.8,.2,0],[.1,.6,.3]])
        state=np.array([0,1,1]);depth=np.array([10,20,80])
        m=truth.dm_population_moments(mu,state,depth)
        expected_p=sum(mu[s] for s in state)/3
        expected_c=sum(d*mu[s] for s,d in zip(state,depth))/3
        np.testing.assert_allclose(m['p'],expected_p)
        np.testing.assert_allclose(m['mean_counts'],expected_c)
        np.testing.assert_allclose(m['b'],expected_c/depth.mean())
        self.assertFalse(np.allclose(m['p'],m['b']))

    def test_program_closure_sign_and_conditional_null_support(self):
        mu=np.array([[.2,.3,.1,.15,.25],[.1,.2,.3,.25,.15]])
        before=mu.copy();up=np.array([0,1]);down=np.array([2,3])
        np.testing.assert_array_equal(truth.transfer_program(mu,up,down,0.),mu)
        for strength in (.1,.3,.6):
            q=truth.transfer_program(mu,up,down,strength)
            np.testing.assert_allclose(q.sum(1),1,rtol=0,atol=1e-15)
            self.assertTrue(np.all(q[:,up]>mu[:,up]))
            self.assertTrue(np.all(q[:,down]<mu[:,down]))
            np.testing.assert_allclose(q[:,4],mu[:,4],rtol=0,atol=1e-15)
            self.assertTrue(np.all(q>=0))
        np.testing.assert_array_equal(mu,before)

    def test_same_aggregate_mean_is_not_a_marginal_null_law(self):
        a=np.array([[.8,.2],[.2,.8]]);b=np.array([[.5,.5],[.5,.5]])
        np.testing.assert_array_equal(a.mean(0),b.mean(0))
        self.assertFalse(np.any(np.all(np.abs(a-b)<=truth.GENE_NULL_TOL,axis=0)))

    def test_simulator_reproducibility_integer_depths_and_zero_support(self):
        mu=np.array([[.7,.2,.1,0],[.1,.5,.4,0]])
        state=np.array([0,1,1]);depth=np.array([30,70,140])
        x,atom=truth.simulate_dm(mu,state,depth,60,20000.,123)
        y,atom2=truth.simulate_dm(mu,state,depth,60,20000.,123)
        np.testing.assert_array_equal(x,y);np.testing.assert_array_equal(atom,atom2)
        self.assertEqual(x.dtype,np.int32);self.assertTrue(np.all(x>=0))
        np.testing.assert_array_equal(x.sum(1),depth[atom]);np.testing.assert_array_equal(x[:,3],0)
        z,_=truth.simulate_dm(mu,state,depth,60,20000.,124)
        self.assertFalse(np.array_equal(x,z))

    def test_dm_nonzero_reference_matches_conditional_mean_and_variance(self):
        mu0=np.array([[.7,.2,.1,0],[.1,.5,.4,0]])
        mu=truth.transfer_program(mu0,[0],[1],.6)
        state=np.array([0,1,1]);depth=np.array([30,70,140]);kappa=50.;n=12000
        x,atoms=truth.simulate_dm(mu,state,depth,n,kappa,410)
        L=depth[atoms].astype(float);p=mu[state[atoms]]
        expected=L[:,None]*p
        variance=L[:,None]*p*(1-p)*(L[:,None]+kappa)/(1+kappa)
        residual=x-expected
        bound=6*np.sqrt(variance.sum(0))/n+1e-12
        self.assertTrue(np.all(np.abs(residual.mean(0))<=bound))
        positive=variance.mean(0)>0
        ratio=(residual*residual).mean(0)[positive]/variance.mean(0)[positive]
        self.assertTrue(np.all((ratio>.9)&(ratio<1.1)))
        prob_error=(x/L[:,None]-p).mean(0)
        prob_bound=6*np.sqrt((variance/L[:,None]**2).sum(0))/n+1e-12
        self.assertTrue(np.all(np.abs(prob_error)<=prob_bound))

    def test_program_selection_is_control_only_disjoint_and_repeatable(self):
        genes=np.asarray([f'g{i:03d}' for i in range(300)])
        mean=np.linspace(10,90,300);detection=np.full(300,.25)
        a=truth.select_program(genes,mean,detection,'BROAD')
        b=truth.select_program(genes,mean,detection,'ZERO_RICH',excluded=np.r_[*a])
        self.assertEqual(len(set(np.r_[*a])),64);self.assertEqual(len(set(np.r_[*b])),64)
        self.assertFalse(set(np.r_[*a])&set(np.r_[*b]))
        for x,y in zip(a,truth.select_program(genes,mean,detection,'BROAD')):np.testing.assert_array_equal(x,y)
        with self.assertRaises(ValueError):truth.select_program(genes[:20],mean[:20],detection[:20],'BROAD')

    def test_state_projection_uses_frozen_inference_parameters(self):
        rng=np.random.default_rng(7);x=sparse.csr_matrix(rng.poisson(5,(90,12)))
        parameters,scores=truth.fit_control_states(x)
        np.testing.assert_allclose(truth.project_states(x,parameters),scores,atol=1e-12,rtol=1e-12)
        old=parameters['split_medians'].copy()
        truth.project_states(sparse.csr_matrix(rng.poisson(30,(30,12))),parameters)
        np.testing.assert_array_equal(parameters['split_medians'],old)

    def test_invalid_counts_and_probabilities_fail_loudly(self):
        for x in ([[1,.5]],[[1,-1]],[[1,np.nan]],[[1,2**40]],[[0,0]]):
            with self.assertRaises(ValueError):truth.checked_counts(x)
        with self.assertRaises(ValueError):truth.fisher_allocation(1,1,3,0)
        with self.assertRaises(ValueError):truth.transfer_program([[.5,.5]],[0,0],[1],.3)
        with self.assertRaises(ValueError):truth.simulate_dm([[.5,.5]],[0],[10],2,0,1)

    def test_seed_domains_and_no_candidate_imports(self):
        names=['dm:inference','dm:comparison']+[f'dm:DM_T{t}_{r}:reference' for t in (1,2) for r in truth.REGIMES]
        seeds=[truth.named_seed(seed,name) for seed in truth.SEEDS for name in names]
        self.assertEqual(len(seeds),len(set(seeds)))
        tree=ast.parse(Path(truth.__file__).read_text())
        imported=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):imported.extend(a.name for a in node.names)
            if isinstance(node,ast.ImportFrom):imported.append(node.module or '')
        self.assertFalse(set(imported)&{'control_preserving_emitter','build_prediction','model','predict','stratified_donor_sampling'})


if __name__=='__main__':unittest.main()
