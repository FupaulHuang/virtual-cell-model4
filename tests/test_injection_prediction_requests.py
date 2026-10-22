"""Anchored oracle requests, exact null replay and complete-method provenance."""
from pathlib import Path
import contextlib
import io
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import build_injection_predictions as build


def anchor(inference,oracle,baseline,is_null=False):
    return build.anchor_requests(inference,inference,oracle,oracle,baseline,baseline,
                                 np.asarray(oracle)-baseline,np.asarray(oracle)-baseline,is_null)


class InjectionRequestTests(unittest.TestCase):
    def test_null_uses_own_inference_baseline_without_renormalizing(self):
        own=np.array([.1,.2,.7000000000000001]);truth=np.array([.2,.3,.5])
        table,audit=anchor(own,truth,truth,True)
        np.testing.assert_array_equal(table['request_p'],own)
        np.testing.assert_array_equal(table['request_b'],own)
        self.assertIsNot(table['request_p'],own)
        self.assertTrue(audit['zero_delta']);self.assertTrue(audit['p']['exact_inference_baseline_copy'])
        self.assertFalse(np.array_equal(table['request_p'],truth))

    def test_clipping_and_closure_drift_are_explicit_and_do_not_mutate_inputs(self):
        own=np.array([.1,.9]);baseline=np.array([.4,.6]);oracle=np.array([.1,.9])
        saved=[x.copy() for x in (own,baseline,oracle)]
        table,audit=anchor(own,oracle,baseline)
        np.testing.assert_allclose(table['anchored_before_clip_p'],[-.2,1.2],atol=1e-15)
        np.testing.assert_allclose(table['clipped_negative_p'],[.2,0],atol=1e-15)
        np.testing.assert_array_equal(table['request_p'],[0,1])
        np.testing.assert_allclose(table['request_delta_error_p'],[.2,-.2],atol=1e-15)
        self.assertAlmostEqual(audit['p']['clipped_negative_mass'],.2)
        self.assertAlmostEqual(audit['p']['post_clip_mass'],1.2)
        for old,new in zip(saved,(own,baseline,oracle)):np.testing.assert_array_equal(old,new)

    def test_no_clip_anchors_oracle_delta_on_inference_not_oracle_absolute(self):
        own=np.array([.3,.2,.5]);baseline=np.array([.2,.3,.5]);oracle=np.array([.15,.31,.54])
        table,audit=anchor(own,oracle,baseline)
        np.testing.assert_allclose(table['request_p'],own+oracle-baseline,atol=1e-15)
        np.testing.assert_allclose(table['request_delta_error_p'],0,atol=1e-15)
        self.assertFalse(np.array_equal(table['request_p'],oracle))
        self.assertEqual(audit['p']['clipped_genes'],0)

    def test_nonzero_null_or_tampered_delta_rejected(self):
        own=np.array([.25,.25,.5]);oracle=np.array([.2,.3,.5])
        with self.assertRaisesRegex(ValueError,'NULL'):anchor(own,oracle,own,True)
        with self.assertRaisesRegex(ValueError,'Declared oracle delta'):
            build.anchor_requests(own,own,oracle,oracle,own,own,np.zeros(3),oracle-own,False)
        with self.assertRaises(ValueError):anchor(np.array([-.1,1.1]),np.array([.5,.5]),np.array([.5,.5]))

    def test_seed_domains_pair_samplers_but_separate_aliases_and_templates(self):
        a=build.seed_policy('composition','COMP_T1_NULL',20261026)
        b=build.seed_policy('composition','COMP_T1_STRONG',20261026)
        self.assertEqual(a['template_seed'],b['template_seed'])
        self.assertNotEqual(a['donor_seed'],b['donor_seed'])
        self.assertNotEqual(a['transport_seed'],a['donor_seed'])
        self.assertEqual(a,build.seed_policy('composition','COMP_T1_NULL',20261026))
        self.assertNotEqual(a,build.seed_policy('composition','COMP_T1_NULL',20261027))

    def test_balanced_requires_final_consistent_gate(self):
        self.assertEqual(build.eligible_methods({'combined_gate_passed':False,'injection_balanced_arm_permitted':False}),list(build.BASE_METHODS))
        self.assertIn('preserving_balanced',build.eligible_methods({'combined_gate_passed':True,'injection_balanced_arm_permitted':True}))
        for record in ({},{'combined_gate_passed':True,'injection_balanced_arm_permitted':False},
                       {'combined_gate_passed':'true','injection_balanced_arm_permitted':True}):
            with self.assertRaises(ValueError):build.eligible_methods(record)

    def test_population_null_fit_and_transport_preserve_selected_raw_counts(self):
        rng=np.random.default_rng(4);raw=rng.multinomial(64,np.full(8,1/8),size=32).astype(np.int32)
        p,b=build.population_profiles(raw);before=raw.copy();previous=sys.getprofile()
        fitted,audit,stages=build.fit_population_audited(raw,p,b)
        self.assertTrue(audit['identity']);self.assertIs(sys.getprofile(),previous)
        rows=np.array([0,3,7,11,19,23,28,31]);counts=build.emit_counts(raw[rows],fitted[rows],123)
        np.testing.assert_array_equal(counts,raw[rows]);np.testing.assert_array_equal(raw,before)
        np.testing.assert_array_equal(stages['bulk_projected_b'],b)

    def test_failed_population_observer_restores_previous_profile(self):
        raw=np.tile([2,2,4],(8,1)).astype(np.int32);p,_=build.population_profiles(raw)
        previous=sys.getprofile()
        with self.assertRaises(ValueError):build.fit_population_audited(raw,p,np.array([.6,.2,.2]))
        self.assertIs(sys.getprofile(),previous)

    def test_failure_types_do_not_relabel_input_errors_as_feasibility(self):
        self.assertEqual(build.failure_record(ValueError('Population fit failed: mean=.01, bulk=.02'))['failure_kind'],'FEASIBILITY')
        self.assertEqual(build.failure_record(ValueError('Invalid probability profile'))['failure_kind'],'ENGINEERING_ERROR')

    def test_pickle_free_control_ids_and_partial_methods_are_unscored(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);out=root/'outputs';out.mkdir();family='composition';seed=20261026;context='COMPOSITION'
            genes=np.array(['G1','G2','G3']);aliases=np.array([f'COMP_TEST_{i}' for i in range(8)])
            raw=np.tile([2,2,4],(4000,1)).astype(np.int32)
            obs=pd.DataFrame({'target_gene':'non-targeting','context':context,
                'original_source_cell_id':[f'original_{i}' for i in range(4000)]},index=[f'cell_{i}' for i in range(4000)])
            controls=root/'inference.h5ad';ad.AnnData(sparse.csr_matrix(raw),obs=obs,var=pd.DataFrame(index=genes)).write_h5ad(controls)
            gene_csv=root/'genes.csv';target_csv=root/'targets.csv'
            pd.DataFrame({'gene_name':genes}).to_csv(gene_csv,index=False)
            pd.DataFrame({'target_gene':aliases,'n_cells':400}).to_csv(target_csv,index=False)
            spec={'context':context,'genes_path':str(gene_csv),'targets_path':str(target_csv),'aliases':aliases.tolist(),
                  'seeds':{str(seed):{f'{role}_path':str(controls) for role in ('inference','comparison','reference')}}}
            spec['seeds'][str(seed)].update({f'{role}_sha256':build.sha256(controls) for role in ('inference','comparison','reference')})
            actual,ids,original=build.load_inference(spec,seed,genes)
            self.assertEqual(ids.dtype.kind,'U');self.assertEqual(original.dtype.kind,'U')
            own_p,own_b=build.population_profiles(actual);tables=[];audits=[]
            for alias in aliases:
                table,audit=build.anchor_requests(own_p,own_b,own_p,own_b,own_p,own_b,np.zeros(3),np.zeros(3),True)
                tables.append(table);audits.append({'alias':str(alias),**audit})
            arrays={k:np.stack([t[k] for t in tables]) for k in tables[0]}
            arrays.update(aliases=aliases,genes=genes,is_null=np.ones(8,dtype=bool),inference_cell_ids=ids,original_source_cell_ids=original)
            requests=root/'requests.npz';np.savez_compressed(requests,**arrays)
            with np.load(requests,allow_pickle=False) as d:
                for key in d.files:_=d[key]
            audit_path=root/'request_audit.json';audit_path.write_text(json.dumps(audits))
            s=build.seed_policy(family,str(aliases[0]),seed)
            template,depth,_,_=build.control_template(controls,genes,400,4,s['template_seed'])
            freeze=out/'freeze.json';freeze.write_text('{}')
            truth_path=root/'truth.json';truth_path.write_text('{}')
            eligibility=root/'eligibility.json';eligibility.write_text('{}')
            prepared=out/'prepared.json';prepared.write_text(json.dumps({'freeze_sha256':build.sha256(freeze),'groups':{f'{family}:{seed}':{
                'requests_path':str(requests),'requests_sha256':build.sha256(requests),'request_audit_path':str(audit_path),
                'request_audit_sha256':build.sha256(audit_path),'inference_sha256':build.sha256(controls),
                'raw_counts_array_sha256':build.array_sha256(raw),'template_seed':s['template_seed'],
                'template_sha256':build.array_sha256(template),'template_depths_sha256':build.array_sha256(depth)}}}))
            frozen={'methods':list(build.BASE_METHODS)};truth={'families':{family:spec}}
            original_emit=build.emit_audited;bad_seed=build.seed_policy(family,str(aliases[2]),seed)['transport_seed']
            # Match the frozen emit_audited positional signature without changing kernels.
            def intercept(template,depth,p,b,emission_seed):
                if emission_seed==bad_seed:raise ValueError('Requested bulk projection too large: 0.04')
                return original_emit(template,depth,p,b,emission_seed)
            with patch.multiple(build,ROOT=root,OUT=out,FREEZE=freeze,PREPARED=prepared,TRUTH=truth_path,ELIGIBILITY=eligibility),\
                 patch.object(build,'verify_freeze',return_value=(frozen,truth)),\
                 patch.object(build,'emit_audited',side_effect=intercept),contextlib.redirect_stdout(io.StringIO()):
                build.generate(family,seed)
            current=json.loads((out/family/f'current_s{seed}'/'manifest.json').read_text())
            preserving=json.loads((out/family/f'preserving_uniform_s{seed}'/'manifest.json').read_text())
            self.assertFalse(current['complete']);self.assertFalse(current['DE_permitted']);self.assertIsNone(current['prediction_path'])
            self.assertEqual(current['failures'][0]['alias'],aliases[2]);self.assertEqual(len(current['successful_aliases']),7)
            self.assertTrue(Path(current['partial_prediction_path']).is_file())
            self.assertTrue(preserving['complete']);self.assertTrue(preserving['DE_permitted'])
            pooled=pd.read_csv(current['selected_rows_path'])
            self.assertEqual(len(pooled),8*1600);self.assertEqual(set(pooled.pool_member),{0,1,2,3})
            self.assertEqual(set(pooled.selection_role),{'current_pooled_template'})
            records=json.loads(Path(preserving['alias_audit_path']).read_text())
            self.assertTrue(all(x['zero_delta_raw_identity'] for x in records))
            self.assertEqual(len(records),8)
            current_records=json.loads(Path(current['alias_audit_path']).read_text())
            for a,b in zip(current_records,records):
                self.assertEqual(a['request_p_sha256'],b['request_p_sha256']);self.assertEqual(a['request_b_sha256'],b['request_b_sha256'])


if __name__=='__main__':unittest.main()
