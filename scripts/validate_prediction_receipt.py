#!/usr/bin/env python3
"""Independently validate emitted axes, all counts, group sizes and provenance."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def validate(prediction, receipt_path, *, excluded_source_tokens=(), verify_hash=True):
    prediction, receipt_path = Path(prediction), Path(receipt_path)
    receipt = json.loads(receipt_path.read_text())
    args = receipt['arguments']
    genes = pd.read_csv(args['genes'], dtype=str).gene_name.to_numpy()
    full_targets = pd.read_csv(args['targets'], dtype=str).target_gene.to_numpy()
    targets = full_targets[:args['target_limit']] if args.get('target_limit') else full_targets
    contexts = [spec.split('=', 1)[0] for spec in args['context']]
    for name, labels in [('genes', genes), ('targets', targets), ('contexts', contexts)]:
        if len(labels) != len(set(labels)) or not len(labels):
            raise ValueError(f'Duplicate or empty {name} axis')
    if list(targets) != receipt['targets'] or len(genes) != receipt['output_genes']:
        raise ValueError('Receipt axes disagree with input axes')
    for spec in args['context']:
        context, path = spec.split('=', 1)
        control = ad.read_h5ad(path, backed='r')
        try:
            if set(control.obs.context.astype(str)) != {context}:
                raise ValueError(f'Control context labels do not match {context}')
            if set(control.obs.target_gene.astype(str)) != {'non-targeting'}:
                raise ValueError('Control source contains perturbed cells')
            if not np.array_equal(control.var_names.astype(str), genes):
                raise ValueError('Control gene order differs')
        finally:
            control.file.close()
    for record in receipt['sources']:
        for token in excluded_source_tokens:
            if token.lower() in str(record['path']).lower():
                raise ValueError(f'Forbidden evaluation source {record["path"]}')
    expected_pairs = {(c, t) for c in contexts for t in targets}
    moments = {(m['context'], m['target']): m for m in receipt['moments']}
    if len(moments) != len(receipt['moments']) or set(moments) != expected_pairs:
        raise ValueError('Receipt moment groups are missing, duplicated or unexpected')
    obj = ad.read_h5ad(prediction, backed='r')
    try:
        if not np.array_equal(obj.var_names.astype(str), genes):
            raise ValueError('Prediction gene axis differs from input order')
        if obj.shape != (len(expected_pairs) * args['cells'], len(genes)):
            raise ValueError('Prediction shape disagrees with complete axes')
        row_pairs = list(zip(obj.obs.context.astype(str), obj.obs.target_gene.astype(str)))
        groups = Counter(row_pairs)
        if set(groups) != expected_pairs or any(n != args['cells'] for n in groups.values()):
            raise ValueError('Unexpected context-target groups or per-group cell counts')
        nnz, max_library = 0, 0
        nnz_by_group, max_by_group = Counter(), Counter()
        for start in range(0, obj.n_obs, 1024):
            block = sparse.csr_matrix(obj.X[start:start + 1024])
            if (not np.isfinite(block.data).all() or (block.data <= 0).any()
                    or (block.data != np.floor(block.data)).any()):
                raise ValueError('Stored counts must be finite positive integers without explicit zeros')
            if block.indices.size and (block.indices.min() < 0 or block.indices.max() >= len(genes)):
                raise ValueError('Sparse gene index outside axis')
            depths = np.asarray(block.sum(axis=1, dtype=np.int64)).ravel()
            if (depths <= 0).any() or (depths > 1_000_000).any():
                raise ValueError('Invalid library depth')
            max_library = max(max_library, int(depths.max()))
            nnz += block.nnz
            for offset, (number, depth) in enumerate(zip(np.diff(block.indptr), depths)):
                key = row_pairs[start + offset]
                nnz_by_group[key] += int(number)
                max_by_group[key] = max(max_by_group[key], int(depth))
        if nnz > 4_750_000_000 or obj.n_obs > 400_000:
            raise ValueError('Official size cap exceeded')
        if receipt['cells'] != obj.n_obs or receipt['stored_entries'] != nnz:
            raise ValueError('Receipt counts disagree with actual file')
        for key, audit in moments.items():
            if audit['nnz'] != nnz_by_group[key] or audit['max_library'] != max_by_group[key]:
                raise ValueError('Receipt group diagnostics disagree with actual file')
        shape = list(obj.shape)
    finally:
        obj.file.close()
    if receipt['bytes'] != prediction.stat().st_size:
        raise ValueError('Receipt file size mismatch')
    if verify_hash and receipt['output_sha256'] != sha256(prediction):
        raise ValueError('Receipt file hash mismatch')
    return {'valid': True, 'shape': shape, 'groups': len(expected_pairs),
            'nnz': nnz, 'max_library': max_library, 'hash_verified': verify_hash,
            'excluded_source_tokens': list(excluded_source_tokens)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prediction', type=Path)
    parser.add_argument('--receipt', type=Path)
    parser.add_argument('--exclude-source-token', action='append', default=[])
    args = parser.parse_args()
    result = validate(args.prediction, args.receipt or args.prediction.with_suffix('.receipt.json'),
                      excluded_source_tokens=args.exclude_source_token)
    print(json.dumps(result, indent=2))
