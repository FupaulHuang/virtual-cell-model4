"""Matched representation, metadata freeze, proxy labels and fusion-based gates."""
from pathlib import Path
import copy
import sys
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import response_transfer_adapter as adapter
import benchmark_matched_response as benchmark


def bank(family, factor=1.):
    targets=np.array(['T1','T2','T3']); genes=np.array(['T1','T2','T3','A','B','C'])
    control=np.array([.1,.12,.14,.16,.2,.28])
    changes=np.array([[.04,-.02,.03,-.01,.01,-.05],[-.01,.02,-.03,.02,-.01,.01],[.01,-.02,.01,-.03,.05,-.02]])*factor
    probability=control+changes
    counts=np.rint(probability*100000).astype(np.int64)
    return adapter.RawBank('synthetic_'+family,targets,genes,np.array([100,110,120]),counts,
        probability*1e6,control,control*1e6,np.ones(len(genes),dtype=bool))


def zero_head():
    return {'coefficient':[0.]*len(adapter.FEATURE_NAMES),'scale':[1.]*len(adapter.FEATURE_NAMES),
            'regularization':.03,'confidence':False}


class MatchedResponseTests(unittest.TestCase):
    def test_evaluation_metadata_includes_uncovered_targets_and_ignores_response(self):
        destination=bank('K562')
        first=benchmark.evaluation_pairs(destination,12,5)
        destination.counts *= 3; destination.mean_cpm *= 20
        second=benchmark.evaluation_pairs(destination,12,5)
        for a,b in zip(first,second): np.testing.assert_array_equal(a,b)
        self.assertEqual(set(first[0]),set(destination.targets))
        self.assertTrue(np.all(first[0]!=first[1]))
        for t in destination.targets:
            self.assertAlmostEqual(first[2][first[0]==t].sum(),1/3)

    def test_matching_fusion_is_specific_component_with_shared_fixed_zero(self):
        destination=bank('K562'); sources=[bank('HCT116',.8),bank('HEK293T',1.2)]
        pair=benchmark.evaluation_pairs(destination,15,10)
        task=benchmark.task_from_pairs(destination,sources,'log2fc',pair)
        actual=benchmark.predictions(task,zero_head())
        expected=[]
        for source in sources:
            data=source.pair_data(pair[0],pair[1],'log2fc',100000.,20)
            expected.append(data['effect']-data['common'])
        # Frozen feature pooling sums float32 effects before the float64 common
        # subtraction; independent center-then-average differs by a few ulps.
        np.testing.assert_allclose(actual['matched_fusion'],np.mean(expected,axis=0),rtol=1e-6,atol=2e-8)
        np.testing.assert_array_equal(actual['zero'],0)
        np.testing.assert_array_equal(actual['adapter_v1'],0)

    def test_destination_response_not_in_features_and_self_family_rejected(self):
        destination=bank('K562'); sources=[bank('HCT116')]
        pair=benchmark.evaluation_pairs(destination,15,10)
        before=benchmark.task_from_pairs(destination,sources,'bulk_delta',pair)
        destination.counts=destination.counts[:,::-1].copy();destination._cache.clear()
        after=benchmark.task_from_pairs(destination,sources,'bulk_delta',pair)
        np.testing.assert_array_equal(before['x'],after['x'])
        self.assertFalse(np.array_equal(before['y'],after['y']))
        with self.assertRaises(ValueError): benchmark.task_from_pairs(destination,[destination],'bulk_delta',pair)

    def test_proxy_thresholds_are_inclusive_and_sign_quality_is_explicit(self):
        targets=np.array(['T']*4)
        truth=np.array([.25,-.3,.1,.4]); predicted=np.array([.25,.3,.3,0.])
        _,p=benchmark.target_primitives(targets,truth,predicted,np.ones(4),.25)
        result=benchmark.summarize_primitives(p)
        self.assertAlmostEqual(result['reference_call_fraction'],.75)
        self.assertAlmostEqual(result['predicted_call_fraction'],.75)
        self.assertAlmostEqual(result['signed_precision'],1/3)
        self.assertAlmostEqual(result['signed_recovery'],1/3)
        self.assertAlmostEqual(result['joint_call_sign_accuracy'],.5)
        self.assertEqual(result['stable_false_call_fraction'],1.)
        self.assertAlmostEqual(result['unweighted_mse'],np.mean((truth-predicted)**2))

    def test_zero_calls_have_undefined_precision_not_perfect_signs(self):
        _,p=benchmark.target_primitives(np.array(['T','T']),np.array([.5,0.]),np.zeros(2),np.ones(2),.25)
        result=benchmark.summarize_primitives(p)
        self.assertIsNone(result['signed_precision']);self.assertIsNone(result['joint_call_sign_accuracy'])
        self.assertEqual(result['signed_recovery'],0.)

    def test_target_balancing_is_not_response_class_balancing(self):
        targets=np.array(['A','A','A','B'])
        truth=np.array([0.,0.,0.,2.]);predicted=np.zeros(4)
        _,p=benchmark.target_primitives(targets,truth,predicted,np.ones(4),.25)
        result=benchmark.summarize_primitives(p)
        self.assertEqual(result['unweighted_mse'],2.)
        self.assertEqual(result['stable_mse'],0.)
        self.assertEqual(result['responsive_mse'],4.)

    def test_bootstrap_withholds_undefined_ratio_intervals_without_dropping(self):
        targets=np.array(['A','B']); truth=np.ones(2)
        _,base=benchmark.target_primitives(targets,truth,np.ones(2),np.ones(2),.25)
        _,candidate=benchmark.target_primitives(targets,truth,np.array([1.,0.]),np.ones(2),.25)
        result=benchmark.paired_intervals(base,candidate,draws=100,seed=2)
        self.assertGreater(result['signed_precision']['undefined_draws'],0)
        self.assertIsNone(result['signed_precision']['candidate_minus_fusion_95_ci'])
        self.assertIsNotNone(result['unweighted_mse']['candidate_minus_fusion_95_ci'])
        self.assertEqual(result,benchmark.paired_intervals(base,candidate,draws=100,seed=2))

    def test_advancement_requires_recovery_against_fusion_not_zero(self):
        contexts=[]
        for _ in range(3):
            contexts.append({s:{'matched_fusion':{'signed_precision':.6,'signed_recovery':.4,'unweighted_mse':.2},
                               'adapter_v1':{'signed_precision':.7,'signed_recovery':.1,'unweighted_mse':.1}}
                             for s in adapter.SPACES})
        result=benchmark.advancement(contexts)
        self.assertFalse(result['both_heads_pass']);self.assertTrue(result['under_recovery_established'])
        for c in contexts:
            for s in adapter.SPACES: c[s]['adapter_v1']['signed_recovery']=.5
        self.assertTrue(benchmark.advancement(contexts)['both_heads_pass'])
        contexts[0]['bulk_delta']['adapter_v1']['signed_precision']=.5
        self.assertFalse(benchmark.advancement(contexts)['both_heads_pass'])
        contexts[0]['bulk_delta']['adapter_v1']['signed_precision']=None
        self.assertFalse(benchmark.advancement(contexts)['both_heads_pass'])

    def test_full_axis_inverse_not_sampled_gene_renormalization(self):
        destination=bank('K562'); sources=[bank('HCT116'),bank('HEK293T')]
        heads={s:zero_head() for s in adapter.SPACES}; shapes=[]
        original=adapter.effects_to_profiles
        def observe(effects,mean,bulk,**kwargs):
            shapes.append(mean.shape)
            return original(effects,mean,bulk,**kwargs)
        with patch.object(adapter,'effects_to_profiles',side_effect=observe):
            rows=benchmark.full_axis_diagnostic(destination,sources,heads,['T1'])
        self.assertEqual(shapes,[(1,6)]*3)
        self.assertEqual({r['scored_genes'] for r in rows},{5})
        a=next(r for r in rows if r['method']=='zero');b=next(r for r in rows if r['method']=='adapter_v1')
        self.assertEqual(a['arithmetic_composition_log2_mse'],b['arithmetic_composition_log2_mse'])
        self.assertEqual(a['bulk_log1p_mse'],b['bulk_log1p_mse'])


if __name__=='__main__': unittest.main()
