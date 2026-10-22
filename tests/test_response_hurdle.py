"""Prior weights, monotone recalibration, no-call boundaries and LOSO guards."""
from pathlib import Path
import copy
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import response_transfer_adapter as adapter
import benchmark_matched_response as p1
import benchmark_response_hurdle as p2


def task(context,seed):
    rng=np.random.default_rng(seed);targets=np.repeat([f'T{i}' for i in range(20)],20)
    x=rng.normal(size=(len(targets),18));x[:,0]*=.3;x[:,9]*=.1
    confidence=rng.normal(size=(len(targets),10));confidence[:,0]=1
    responsive=np.abs(x[:,0])>.2
    y=np.where(responsive,np.where(x[:,0]>=0,1.,-1.)*(.3+np.abs(x[:,0])),.04*x[:,1])
    return {'x':x,'confidence_x':confidence,'y':y,'weight':np.full(len(targets),1/len(targets)),
            'targets':targets,'genes':np.tile([f'G{i}' for i in range(20)],20),
            'coverage':np.ones(len(targets)), 'destination':context,'input_families':['other']}


def bank(family,seed):
    rng=np.random.default_rng(seed);targets=np.asarray([f'T{i}' for i in range(20)]);genes=np.asarray([f'G{i}' for i in range(10)])
    c=np.full(10,.1);p=np.exp2(rng.normal(0,.6,(20,10)))*c;p/=p.sum(1,keepdims=True)
    return adapter.RawBank('synthetic_'+family,targets,genes,np.full(20,100),np.rint(p*1e6).astype(np.int64),
                           p*1e6,c,c*1e6,np.ones(10,dtype=bool))


class ResponseHurdleTests(unittest.TestCase):
    def test_partition_is_metadata_only_exact_and_shared_across_contexts(self):
        fit,cal=p2.split_targets([f'T{i}' for i in range(20)])
        self.assertEqual(len(cal),4);self.assertEqual(len(fit),16)
        self.assertFalse(set(fit)&set(cal))
        self.assertEqual((fit,cal),p2.split_targets(list(reversed([f'T{i}' for i in range(20)]))))

    def test_subset_weights_keep_context_and_target_balance(self):
        a,b=task('HCT116',1),task('HEK293T',2)
        # Remove half the rows for one target and repair its original target weight.
        keep=np.ones(len(a['targets']),dtype=bool);keep[:10]=False
        for key in ('x','confidence_x','y','weight','targets','genes','coverage'): a[key]=a[key][keep]
        for target in set(a['targets']): a['weight'][a['targets']==target]=1/(20*np.sum(a['targets']==target))
        combined=p2.combine([a,b],['T0','T1'])
        self.assertAlmostEqual(combined['weight'].sum(),1.)
        np.testing.assert_allclose([c['renormalized_weight_mass'] for c in combined['contexts']],[.5,.5],rtol=0,atol=1e-15)
        for t in ['T0','T1']:self.assertAlmostEqual(combined['weight'][combined['targets']==t].sum(),.5)

    def test_logistic_intercept_obeys_weighted_not_raw_prevalence(self):
        x=np.ones((4,2));x[:,1]=0;y=np.array([0.,0.,0.,1.]);w=np.array([.1,.1,.1,.7])
        model=p2.fit_logistic(x,y,w)
        np.testing.assert_allclose(p2.sigmoid(p2.score(model,x)),.7,atol=1e-7)
        self.assertFalse(model['intercept_penalized'])

    def test_negative_platt_slope_becomes_constant_prior_and_makes_no_calls(self):
        logits=np.linspace(-3,3,100);y=(logits<0).astype(float);w=np.ones(100)/100
        calibrator=p2.fit_platt(logits,y,w)
        self.assertEqual(calibrator['slope'],0.)
        probability=p2.calibrated_probability(calibrator,logits)
        np.testing.assert_array_equal(probability,.5)
        prediction,calls=p2.decide(probability,.5,np.ones(100),np.ones(100),np.ones(100),.25)
        self.assertFalse(calls.any());np.testing.assert_array_equal(prediction,0)

    def test_strict_probability_boundary_sign_floor_cap_and_unsupported_pairs(self):
        prior=.2;eps=p2.POLICY['decision_epsilon']
        probability=np.array([prior,prior+eps,prior+eps+1e-10,.9,.9])
        prediction,calls=p2.decide(probability,prior,np.array([1,1,.5,0,1]),
                                  np.array([0,0,-4,20,20]),np.array([1,1,1,1,0]),.25)
        np.testing.assert_array_equal(calls,[False,False,True,True,False])
        np.testing.assert_allclose(prediction,[0,0,.25,-3,0],rtol=0,atol=1e-15)

    def test_degenerate_classes_fail_without_alternative_model(self):
        with self.assertRaises(p2.DegenerateClasses):p2.fit_logistic(np.ones((5,2)),np.zeros(5),np.ones(5))
        with self.assertRaises(p2.DegenerateClasses):p2.fit_platt(np.arange(5),np.zeros(5),np.ones(5))
        tasks=[task('HCT116',1),task('HEK293T',2)]
        for item in tasks:item['y'][:]=0
        _,cal=p2.split_targets(tasks[0]['targets'])
        with self.assertRaises(p2.DegenerateClasses):p2.fit_hurdle(tasks,cal,'log2fc')

    def test_fit_calibration_targets_disjoint_class_weights_and_prediction_no_truth_access(self):
        tasks=[task('HCT116',4),task('HEK293T',5)];_,cal=p2.split_targets(tasks[0]['targets'])
        model=p2.fit_hurdle(tasks,cal,'log2fc')
        self.assertFalse(set(model['fitting_target_ids'])&set(model['recalibration_target_ids']))
        np.testing.assert_allclose(model['balanced_class_weight_mass'],[.5,.5],atol=1e-12)
        self.assertFalse(model['calibrator']['class_balanced_weights'])
        self.assertGreaterEqual(model['calibrator']['slope'],0)
        evaluation=task('K562',6);first=p2.predict_hurdle(model,evaluation)
        evaluation['y']*=100;second=p2.predict_hurdle(model,evaluation)
        for a,b in zip(first,second):np.testing.assert_array_equal(a,b)

    def test_outer_context_poisoning_cannot_change_fitting_recalibration_or_boundary(self):
        banks=[bank(f,i+10) for i,f in enumerate(adapter.FAMILIES)];heldout=banks[0]
        pairs={};training=[]
        for destination,source in [(banks[1],banks[2]),(banks[2],banks[1])]:
            key=destination.family
            p1.put_pairs(pairs,key,adapter._sample_pairs(destination,[source],200,20,20261021))
            training.append({'destination':destination.family,'inputs':[source.family],'pairs_key':key})
        fold={'heldout':heldout.family,'training_tasks':training}
        original=p2.training_tasks(banks,fold,pairs,'log2fc')
        _,cal=p2.split_targets(set(t for x in original for t in x['targets']))
        model=p2.fit_hurdle(original,cal,'log2fc')
        heldout.counts=heldout.counts[:,::-1].copy();heldout.mean_cpm*=20;heldout.control_mean_cpm=heldout.control_mean_cpm*5
        again=p2.training_tasks(banks,fold,pairs,'log2fc')
        self.assertEqual(model,p2.fit_hurdle(again,cal,'log2fc'))
        leaked=copy.deepcopy(fold);leaked['training_tasks'][0]['inputs']=['K562']
        with self.assertRaisesRegex(ValueError,'held-out'):p2.training_tasks(banks,leaked,pairs,'log2fc')

    def test_reliability_weighted_prior_and_fixed_ten_bins(self):
        result=p2.reliability(np.array([0.,1.]),np.array([0.,1.]),np.array([.9,.1]))
        self.assertEqual(result['brier'],0.);self.assertAlmostEqual(result['weighted_prevalence'],.1)
        self.assertEqual(len(result['bins']),10);self.assertEqual(result['bins'][1]['weight_mass'],0)

    def test_control_expression_strata_preserve_masks_weights_and_are_not_gate(self):
        destination=bank('K562',2)
        destination.control_mean_cpm=destination.control_mean_cpm.copy()
        destination.control_mean_cpm[0,:2]=[5.,6.]
        task_data={'targets':np.array(['T0','T0']),'genes':np.array(['G0','G1']),
                   'y':np.array([.5,-.5]),'weight':np.array([.8,.2]),'coverage':np.ones(2)}
        result=p2.control_expression_strata(destination,task_data,np.array([.5,-.5]),np.zeros(2),.25)
        self.assertEqual(result['control_CPM_le5']['pairs'],1)
        self.assertAlmostEqual(result['control_CPM_le5']['original_weight_mass'],.8)
        self.assertEqual(result['control_CPM_gt5']['methods']['matched_fusion']['signed_recovery'],1.)
        self.assertFalse(result['control_CPM_gt5']['selection_or_advancement_gate'])


if __name__=='__main__':unittest.main()
