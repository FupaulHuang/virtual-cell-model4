#!/usr/bin/env python3
"""Freeze current and learned response requests before cell-level evaluation."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(key, '1')
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse
from build_prediction import ROOT, atlas, axis_csv, control_template, sha256
from build_bulk_ablation import make_profiles, prepare_effects, verify_inputs
from response_transfer_adapter import load_raw_bank, predict_effects, effects_to_profiles


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--adapter', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists() or args.output.with_suffix('.json').exists():
        raise FileExistsError(args.output)
    config_path = ROOT/'configs/dispersion_h1_400_incumbent.json'
    config = json.loads(config_path.read_text())
    incumbent = ROOT/'outputs/dispersion_h1_400_incumbent/prediction.h5ad'
    receipt = json.loads(incumbent.with_suffix('.receipt.json').read_text())
    verify_inputs(config, incumbent, receipt)
    genes, targets = axis_csv(ROOT/config['genes'], 'gene_name'), axis_csv(ROOT/config['targets'], 'target_gene')
    context_path = ROOT/config['context'][0].split('=',1)[1]
    control = ad.read_h5ad(context_path)
    # Preserve the frozen incumbent's reduction order, not merely mathematically
    # equivalent moments, so current/current can be checked count for count.
    _, _, mean, bulk = control_template(context_path, genes, 400, 4, config['seed'])
    if not np.array_equal(control.var_names, genes): raise ValueError('Control gene order differs')
    sources = [atlas.load_xatlas(ROOT/s['path'], prior_counts=config['prior_counts']) for s in config['source']]
    effects = prepare_effects(sources, [s['weight'] for s in config['source']], targets, genes, 'incumbent')
    cd4 = atlas.aligned_cd4(ROOT/config['cd4'], targets, genes, common_subtract=1., center_scope='source')
    current_p, current_b = make_profiles(effects, mean, bulk, cd4, config['cd4_weight'], 'incumbent')
    arrays = dict(genes=genes, targets=targets, control_mean=mean, control_bulk=bulk,
                  current_mean=current_p, current_bulk=current_b)
    files = [str(config_path.relative_to(ROOT)), str(context_path.relative_to(ROOT)),
             config['genes'], config['targets'], config['cd4'], *[s['path'] for s in config['source']],
             'scripts/prepare_redesign_profiles.py', 'scripts/build_bulk_ablation.py',
             'scripts/build_prediction.py', 'external/atlasshift/model.py', 'external/atlasshift/predict.py']
    info = dict(created_at_utc=datetime.now(timezone.utc).isoformat(), h1_responses_excluded=True,
                targets=targets.tolist(), output_genes=len(genes),
                current_response='Frozen incumbent requested mean/bulk with CD4',
                adapter_response=None, reserved_responses_accessed=False)
    if args.adapter:
        model = json.loads(args.adapter.read_text())
        banks = [load_raw_bank(ROOT/s['path']) for s in config['source']]
        adapted, confidence = predict_effects(model, banks, targets, genes, mean*1e6)
        arrays['adapter_mean'], arrays['adapter_bulk'] = effects_to_profiles(adapted, mean*1e6, bulk)
        for space in adapted:
            arrays['adapter_effect_'+space], arrays['adapter_coverage_'+space] = adapted[space]
            arrays['adapter_sign_confidence_'+space] = confidence[space]
        files.extend([str(args.adapter.resolve().relative_to(ROOT)), 'scripts/response_transfer_adapter.py'])
        info['adapter_response'] = dict(model=str(args.adapter.resolve()), amplitude=1., clip=3.,
            cd4_included=False, exact_pseudocount_inverse=True,
            target_coverage=int((adapted['log2fc'][1] > 0).any(1).sum()),
            gene_target_coverage=float((adapted['log2fc'][1] > 0).mean()),
            interpretation='Response intervention includes learned effects, exact units and CD4 omission')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **arrays)
    info.update(input_hashes={name:sha256(ROOT/name) for name in files}, profiles_sha256=sha256(args.output))
    args.output.with_suffix('.json').write_text(json.dumps(info, indent=2)+'\n')
    print(json.dumps(info), flush=True)


if __name__ == '__main__': main()
