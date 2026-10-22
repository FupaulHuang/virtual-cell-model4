#!/usr/bin/env python3
"""Append measured controls to a local scoring copy for cell-eval2 0.16.0.

Its initial validate_pair requires the control label on both sides, even when
control_source=real subsequently substitutes the real controls during DE.
Never use this local evaluation copy as a challenge submission.
"""
import argparse
import json
from pathlib import Path
import shutil
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from validate_prediction_receipt import sha256


def prepare(prediction, reference, output):
    prediction, reference, output = map(lambda p: Path(p).resolve(), (prediction, reference, output))
    if output.exists():
        raise FileExistsError(output)
    pred, real = ad.read_h5ad(prediction, backed='r'), ad.read_h5ad(reference, backed='r')
    try:
        if not np.array_equal(pred.var_names, real.var_names):
            raise ValueError('Prediction and reference genes differ')
        if 'non-targeting' in set(pred.obs.target_gene.astype(str)):
            raise ValueError('Prediction already contains controls; use it directly for local scoring')
        control_rows = np.flatnonzero(real.obs.target_gene.astype(str).to_numpy() == 'non-targeting')
        if not len(control_rows):
            raise ValueError('Reference has no controls')
        controls = real[control_rows].to_memory()
        raw = sparse.csr_matrix(controls.X)
        if raw.nnz and ((raw.data < 0).any() or (raw.data != np.floor(raw.data)).any()):
            raise ValueError('Reference controls are not raw counts')
        obs = pd.concat([pred.obs.copy(), controls.obs.copy()], axis=0, join='outer')
        obs.index = [f'local_eval_{i}' for i in range(len(obs))]
        # Avoid object columns containing mixed NaNs/strings in HDF5 serialization.
        for column in obs:
            if obs[column].dtype.kind in 'OU' or isinstance(obs[column].dtype, pd.CategoricalDtype):
                obs[column] = pd.Categorical(obs[column].astype('string').fillna(''))
        original_shape = pred.shape
    finally:
        pred.file.close()
        real.file.close()
    output.parent.mkdir(parents=True, exist_ok=True)
    pending = output.with_suffix('.partial.h5ad')
    if pending.exists():
        raise FileExistsError(pending)
    shutil.copyfile(prediction, pending)
    with h5py.File(pending, 'r+') as handle:
        x = handle['X']
        if x.attrs['encoding-type'] != 'csr_matrix':
            raise ValueError('Expected CSR prediction')
        nnz = len(x['data'])
        for key, values in [('data', raw.data), ('indices', raw.indices),
                            ('indptr', raw.indptr[1:].astype(np.int64) + nnz)]:
            original = x[key]
            before = len(original)
            if original.maxshape[0] is not None and original.maxshape[0] < before + len(values):
                # Permit older fixed-size H5ADs without changing their source.
                array = original[:]
                dtype = np.int64 if key == 'indptr' else original.dtype
                del x[key]
                original = x.create_dataset(key, shape=(before + len(values),), dtype=dtype,
                    chunks=True, compression='gzip', compression_opts=1, shuffle=True)
                original[:before] = array
            else:
                original.resize((before + len(values),))
            original[before:] = values
        x.attrs['shape'] = (len(obs), original_shape[1])
        del handle['obs']
        ad.io.write_elem(handle, 'obs', obs)
    pending.rename(output)
    receipt = {'purpose': 'Local evaluation only; contains measured controls and cannot be submitted',
               'prediction': str(prediction), 'prediction_sha256': sha256(prediction),
               'reference': str(reference), 'reference_sha256': sha256(reference),
               'controls_added': len(control_rows), 'shape': [len(obs), original_shape[1]],
               'control_source_for_metrics': 'real'}
    output.with_suffix('.receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prediction', type=Path)
    parser.add_argument('--reference', type=Path, default=Path('evaluations/h1_validation/reference.h5ad'))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    prepare(args.prediction, args.reference, args.output)
