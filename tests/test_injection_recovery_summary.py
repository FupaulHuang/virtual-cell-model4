"""Scientific denominator, stage separation and profile-audit tests."""
from pathlib import Path
import sys
import tempfile
import json
import unittest
from unittest.mock import patch
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import summarize_injection_recovery as summary
from prepare_injection_recovery import population_moments,sha


class InjectionSummaryTests(unittest.TestCase):
    def test_failed_de_cannot_be_reclassified_as_incomplete_candidate(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);folder=root/'de'/'composition'/'seed_20261026'/'current';folder.mkdir(parents=True)
            (folder/'status.json').write_text(json.dumps({'status':'failed','family':'composition','seed':20261026,'method':'current'}))
            with patch.object(summary,'TRUTH',root):
                with self.assertRaisesRegex(ValueError,'DE failed'):summary.read_de('composition',20261026,'current',['A'])

    def test_detection_change_and_donor_edits_are_not_counterfactual_truth(self):
        masks={'up':np.array([True,False]),'down':np.array([False,True])}
        d=summary.detection_changes(np.array([.4,.2]),np.array([.2,.3]),np.array([.2,.4]),masks)
        self.assertAlmostEqual(d['up']['reference_detection_change'],.2)
        self.assertAlmostEqual(d['up']['candidate_detection_change'],0.)
        self.assertAlmostEqual(d['down']['candidate_reference_detection_mae'],.1)
        donor=np.array([[0,1],[1,0]]);pred=np.array([[2,0],[1,0]])
        edits=summary.donor_zero_edits(donor,pred,masks)
        self.assertEqual(edits['groups']['up']['zero_to_positive_entries'],1)
        self.assertEqual(edits['groups']['down']['positive_to_zero_entries'],1)
        identity=summary.donor_zero_edits(donor,donor,masks)
        self.assertEqual(identity['groups']['up']['zero_to_positive_entries'],0)

    def test_development_and_confirmation_partitions_are_disjoint_complete(self):
        actual={(s,t):summary.stage_partition(s,t) for s in (20261026,20261027) for t in (1,2)}
        self.assertEqual(actual[(20261026,1)],'development')
        self.assertEqual(set(actual.values()),{'development','heldout_template','heldout_seed','heldout_template_and_seed'})
        with self.assertRaises(ValueError):summary.stage_partition(20261028,1)

    def test_small_effects_are_not_falsely_labeled_null_or_direction_errors(self):
        truth=np.array([0.,.05,.2,-.4,.6]);observed=np.array([2.,-2.,.2,.4,-.6])
        q=np.full(5,.001);gate=np.ones(5,bool);null=np.array([True,False,False,False,False])
        m=summary.side_metrics(truth,observed,q,gate,null)
        self.assertEqual(m['significant_calls'],5)
        self.assertEqual(m['called_below_adjudication_threshold'],2)
        self.assertEqual(m['marginal_null_calls'],1)
        self.assertEqual(m['oracle_adjudicable_gated'],3)
        self.assertEqual(m['called_with_correct_oracle_sign'],1)
        self.assertAlmostEqual(m['signed_recovery'],1/3)
        self.assertAlmostEqual(m['direction_accuracy_among_adjudicable_calls'],1/3)

    def test_strict_q_gate_and_out_of_gate_true_signal_are_explicit(self):
        m=summary.side_metrics([.1,.25,.5],[.1,.25,.5],[.05,.049,.001],[True,True,False],[False]*3)
        self.assertEqual(m['significant_calls'],1)
        self.assertEqual(m['oracle_adjudicable_all_genes'],3)
        self.assertEqual(m['oracle_adjudicable_gated'],2)
        self.assertEqual(m['oracle_adjudicable_outside_gate'],1)
        self.assertEqual([b['oracle_genes'] for b in m['effect_strata']],[0,0,1,1,0])

    def test_no_calls_or_true_effects_keep_undefined_denominators(self):
        m=summary.side_metrics([0,0],[0,0],[1,1],[True,True],[True,True])
        self.assertIsNone(m['signed_recovery'])
        self.assertIsNone(m['direction_accuracy_among_adjudicable_calls'])
        overlap=summary.call_overlap([1,1],[1,1],[True,True])
        self.assertIsNone(overlap['jaccard'])
        self.assertIsNone(overlap['fraction_candidate_calls_in_reference'])

    def test_overlap_uses_observed_reference_without_claiming_oracle_truth(self):
        m=summary.call_overlap([.001,.001,1],[.001,1,.001],[True]*3)
        self.assertEqual(m['intersection'],1);self.assertEqual(m['union'],3)
        self.assertAlmostEqual(m['jaccard'],1/3)
        self.assertAlmostEqual(m['fraction_reference_calls_recovered'],.5)

    def test_profile_errors_are_in_declared_units_and_unavailable_is_not_zero(self):
        p=np.array([.4,.6]);b=np.array([.3,.7])
        m=summary.profile_error(p,b,p,b)
        self.assertEqual(m['arithmetic_cpm_rmse'],0.)
        self.assertEqual(m['bulk_log1p50k_rmse'],0.)
        self.assertIsNone(summary.profile_error([np.nan,np.nan],[np.nan,np.nan],p,b))
        with self.assertRaises(ValueError):summary.profile_error([np.nan,.6],b,p,b)

    def test_profile_audit_checks_actual_integer_counts_and_current_na_stage(self):
        x=np.tile(np.array([[4,3,2,1],[1,4,3,2]],dtype=np.int32),(200,1))
        genes=np.array(['a','b','c','d']);m=population_moments(x)
        bank={'oracle_p':np.array([m['p']]),'oracle_b':np.array([m['b']]),'delta_p':np.zeros((1,4)),'delta_b':np.zeros((1,4))}
        frame=pd.DataFrame({'alias':'TEST','gene':genes})
        for prefix in ('oracle','inference','request','sample_expected','realized','anchored_before_clip','donor_baseline'):
            frame[prefix+'_p']=m['p'];frame[prefix+'_b']=m['b']
        frame['population_fitted_p']=np.nan;frame['population_fitted_b']=np.nan;frame['bulk_projected_b']=m['b']
        obj=ad.AnnData(sparse.csr_matrix(x),obs=pd.DataFrame({'target_gene':['TEST']*400},index=[f'c{i}' for i in range(400)]),var=pd.DataFrame(index=genes))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'profiles.parquet';frame.to_parquet(path,index=False)
            candidate={'gene_profiles_path':str(path),'gene_profiles_sha256':sha(path)}
            result=summary.profile_audit(candidate,'TEST',genes,bank,0,{},obj)
            self.assertIsNone(result['stages']['population_fitted'])
            self.assertEqual(result['stages']['realized']['arithmetic_cpm_rmse'],0.)
            self.assertEqual(result['requested_delta_error_p_l1'],0.)
            obj.X[0,0]+=1
            with self.assertRaises(AssertionError):summary.profile_audit(candidate,'TEST',genes,bank,0,{},obj)


if __name__=='__main__':unittest.main()
