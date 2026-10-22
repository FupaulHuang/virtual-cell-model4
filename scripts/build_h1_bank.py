#!/usr/bin/env python3
"""Build an audited, batch-matched H1 raw-response bank from local 2025 data.

The primary arrays implement AtlasShift's ``load_xatlas`` schema. All 300 public
targets are retained. Repeated control observations in the three public splits
are verified by raw-count fingerprints and included exactly once.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import time

THREADS = int(os.environ.get("H1_BUILDER_THREADS", "2"))
if not 1 <= THREADS <= 8:
    raise ValueError("H1_BUILDER_THREADS must be between 1 and 8")
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = str(THREADS)

import anndata as ad
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
CONTROL = "non-targeting"
SPLITS = (("Training", "train"), ("Validation", "validation"), ("Test", "test"))


def axis_sha(values):
    return hashlib.sha256("\0".join(map(str, values)).encode()).hexdigest()


def file_sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_metadata(input_dir):
    files, target_split, observations, batch_set = [], {}, {}, set()
    genes = None
    for filename_split, public_split in SPLITS:
        path = (input_dir / f"adata_{filename_split}.h5ad").resolve()
        data = ad.read_h5ad(path, backed="r")
        try:
            obs = data.obs
            if not {"target_gene", "batch", "guide_id"}.issubset(obs.columns):
                raise ValueError(f"Missing target/batch/guide metadata: {path}")
            if not obs.index.is_unique or not data.var_names.is_unique:
                raise ValueError(f"Duplicate observation or gene IDs within {path}")
            current_genes = data.var_names.to_numpy(dtype=str)
            if genes is None:
                genes = current_genes
            elif not np.array_equal(genes, current_genes):
                raise ValueError("H1 gene axes differ across splits")
            ids = obs.index.to_numpy(dtype=str)
            labels = obs.target_gene.to_numpy(dtype=str)
            batches = obs.batch.to_numpy(dtype=str)
            guides = obs.guide_id.to_numpy(dtype=str)
            if any(np.any(np.isin(v, ["", "nan", "None"])) for v in (ids, labels, batches)):
                raise ValueError("Empty observation/target/batch label")
            keep = np.ones(len(ids), dtype=bool)
            duplicate_controls = 0
            for i, cell in enumerate(ids):
                metadata = (labels[i], batches[i], guides[i])
                if cell in observations:
                    if observations[cell] != metadata or labels[i] != CONTROL:
                        raise ValueError(f"Conflicting or repeated treated observation: {cell}")
                    keep[i] = False
                    duplicate_controls += 1
                else:
                    observations[cell] = metadata
            for target in sorted(set(labels) - {CONTROL}):
                if target in target_split:
                    raise ValueError(f"Target overlaps public splits: {target}")
                target_split[target] = public_split
            batch_set.update(batches)
            stat = path.stat()
            files.append(dict(
                path=str(path), split=public_split, ids=ids, labels=labels,
                batches=batches, keep=keep,
                provenance=dict(
                    path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                    n_cells=len(ids), n_genes=len(genes),
                    n_controls=int(np.count_nonzero(labels == CONTROL)),
                    duplicate_controls_omitted=duplicate_controls,
                    obs_ids_sha256=axis_sha(ids), gene_axis_sha256=axis_sha(genes),
                ),
            ))
        finally:
            data.file.close()
    if len(target_split) != 300:
        raise ValueError(f"Expected all 300 H1 targets, found {len(target_split)}")
    return files, genes, target_split, np.asarray(sorted(batch_set), dtype=str)


def build(input_dir, output, chunk_size):
    started = time.time()
    files, genes, target_split, batches = audit_metadata(input_dir)
    targets = np.asarray(list(target_split), dtype=str)
    target_lookup = {target: i for i, target in enumerate(targets)}
    batch_lookup = {batch: i for i, batch in enumerate(batches)}
    nt, nb, ng = len(targets), len(batches), len(genes)
    shape = (nt + nb, ng)
    count_sums = np.zeros(shape, dtype=np.int64)
    cpm_sums = np.zeros(shape, dtype=np.float64)
    cpm_square_sums = np.zeros(shape, dtype=np.float64)
    group_n = np.zeros(nt + nb, dtype=np.int64)
    target_batch_n = np.zeros((nt, nb), dtype=np.int64)
    target_batch_library = np.zeros((nt, nb), dtype=np.int64)
    control_fingerprints = {}
    raw_library_total = 0
    unique_library_total = 0
    duplicates_verified = 0
    for source in files:
        path = Path(source["path"])
        data = ad.read_h5ad(path, backed="r")
        labels, ids, keep = source["labels"], source["ids"], source["keep"]
        batch_codes = np.asarray([batch_lookup[b] for b in source["batches"]])
        group_codes = np.asarray([
            nt + batch_codes[i] if target == CONTROL else target_lookup[target]
            for i, target in enumerate(labels)
        ])
        file_library = 0
        file_nnz = 0
        file_unique_library = 0
        try:
            for left in range(0, len(ids), chunk_size):
                right = min(left + chunk_size, len(ids))
                x = sparse.csr_matrix(data.X[left:right])
                if (not np.isfinite(x.data).all() or (x.data < 0).any()
                        or (x.data != np.floor(x.data)).any()):
                    raise ValueError(f"Invalid raw counts: {path}, rows {left}:{right}")
                x = x.astype(np.int64)
                x.sum_duplicates()
                x.eliminate_zeros()
                x.sort_indices()
                depth = np.asarray(x.sum(axis=1)).ravel()
                if (depth <= 0).any():
                    raise ValueError(f"Zero-depth source cell: {path}")
                file_library += int(depth.sum())
                file_nnz += x.nnz
                for local in np.flatnonzero(labels[left:right] == CONTROL):
                    a, b = x.indptr[local:local + 2]
                    digest = hashlib.blake2b(digest_size=32)
                    digest.update(x.indices[a:b].astype("<i8", copy=False).tobytes())
                    digest.update(x.data[a:b].astype("<i8", copy=False).tobytes())
                    fingerprint = digest.digest()
                    cell = ids[left + local]
                    if cell in control_fingerprints:
                        if control_fingerprints[cell] != fingerprint:
                            raise ValueError(f"Repeated control has different raw counts: {cell}")
                        duplicates_verified += 1
                    else:
                        control_fingerprints[cell] = fingerprint
                selected = keep[left:right]
                x = x[selected]
                depth = depth[selected]
                group = group_codes[left:right][selected]
                batch = batch_codes[left:right][selected]
                if not len(group):
                    continue
                file_unique_library += int(depth.sum())
                groups, inverse = np.unique(group, return_inverse=True)
                membership = sparse.csr_matrix(
                    (np.ones(len(group), dtype=np.int64), (inverse, np.arange(len(group)))),
                    shape=(len(groups), len(group)),
                )
                count_sums[groups] += (membership @ x).toarray()
                normalized = x.multiply((1000000.0 / depth)[:, None]).tocsr()
                cpm_sums[groups] += (membership @ normalized).toarray()
                normalized.data **= 2
                cpm_square_sums[groups] += (membership @ normalized).toarray()
                group_n += np.bincount(group, minlength=nt + nb)
                perturbed = group < nt
                np.add.at(target_batch_n, (group[perturbed], batch[perturbed]), 1)
                np.add.at(target_batch_library, (group[perturbed], batch[perturbed]), depth[perturbed])
                if left % (chunk_size * 16) == 0 or right == len(ids):
                    print(f"{path.name}: {right}/{len(ids)} rows audited; "
                          f"{time.time() - started:.1f}s elapsed", flush=True)
        finally:
            data.file.close()
        source["provenance"].update(
            raw_library_total=file_library, unique_library_total=file_unique_library,
            nonzero_counts=file_nnz, all_counts_finite_nonnegative_integer=True,
        )
        raw_library_total += file_library
        unique_library_total += file_unique_library
    if (group_n == 0).any():
        raise ValueError("Missing target or batch control cells")
    np.testing.assert_equal(int(count_sums.sum()), unique_library_total)
    np.testing.assert_equal(target_batch_n.sum(axis=1), group_n[:nt])
    np.testing.assert_equal(target_batch_library.sum(axis=1), count_sums[:nt].sum(axis=1))
    np.testing.assert_allclose(cpm_sums.sum(axis=1) / group_n, 1e6, rtol=1e-12)
    expected_duplicates = sum(s["provenance"]["duplicate_controls_omitted"] for s in files)
    if duplicates_verified != expected_duplicates:
        raise ValueError("Duplicate controls were not all verified")
    means = cpm_sums / group_n[:, None]
    second_moments = cpm_square_sums / group_n[:, None]
    population_variance = np.maximum(second_moments - means**2, 0)
    sample_variance = population_variance * group_n[:, None] / (group_n[:, None] - 1)
    controls = count_sums[nt:]
    control_n = group_n[nt:]
    control_probability = controls / controls.sum(axis=1, keepdims=True)
    cell_weights = target_batch_n / group_n[:nt, None]
    library_weights = target_batch_library / target_batch_library.sum(axis=1, keepdims=True)
    matched_probability = library_weights @ control_probability
    matched_mean = cell_weights @ means[nt:]
    matched_second = cell_weights @ second_moments[nt:]
    matched_mean_variance = (cell_weights**2) @ (sample_variance[nt:] / control_n[:, None])
    global_probability = controls.sum(axis=0) / controls.sum()
    global_mean = cpm_sums[nt:].sum(axis=0) / control_n.sum()
    np.testing.assert_allclose(matched_probability.sum(axis=1), 1, rtol=1e-12)
    np.testing.assert_allclose(matched_mean.sum(axis=1), 1e6, rtol=1e-12)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = dict(
        source=np.asarray("H1_2025_public_batch_matched"), targets=targets, genes=genes,
        n_cells=group_n[:nt], target_count_sums=count_sums[:nt],
        target_mean_cpm=means[:nt].astype(np.float32),
        matched_control_probability=matched_probability.astype(np.float32),
        matched_control_mean_cpm=matched_mean.astype(np.float32),
        global_control_probability=global_probability.astype(np.float32),
        global_control_mean_cpm=global_mean.astype(np.float32),
        measured_genes=np.ones(ng, dtype=bool),
        public_2025_split=np.asarray([target_split[t] for t in targets]),
        target_variance_cpm=sample_variance[:nt].astype(np.float32),
        target_mean_variance_cpm=(sample_variance[:nt] / group_n[:nt, None]).astype(np.float32),
        matched_control_variance_cpm=np.maximum(matched_second - matched_mean**2, 0).astype(np.float32),
        matched_control_mean_variance_cpm=matched_mean_variance.astype(np.float32),
        batches=batches, target_batch_n=target_batch_n,
        target_batch_library=target_batch_library, control_batch_n=control_n,
        control_batch_count_sums=controls,
        control_batch_mean_cpm=means[nt:].astype(np.float32),
        control_batch_variance_cpm=sample_variance[nt:].astype(np.float32),
    )
    for key, value in result.items():
        if value.dtype.kind in "fiu" and (not np.isfinite(value).all() or (value < 0).any()):
            raise ValueError(f"Nonfinite or negative output array: {key}")
    pending = output.with_suffix(output.suffix + ".partial")
    with pending.open("wb") as handle:
        np.savez_compressed(handle, **result)
    pending.replace(output)
    provenance = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.time() - started, source_files=[s["provenance"] for s in files],
        output=str(output.resolve()), output_bytes=output.stat().st_size,
        output_sha256=file_sha(output), script_sha256=file_sha(__file__),
        n_targets=nt, n_genes=ng, n_batches=nb,
        n_treated_cells=int(group_n[:nt].sum()), n_unique_controls=int(control_n.sum()),
        n_duplicate_control_rows_verified_and_omitted=duplicates_verified,
        raw_library_total_including_duplicate_controls=raw_library_total,
        unique_library_total=unique_library_total,
        target_count_total=int(count_sums[:nt].sum()), control_count_total=int(controls.sum()),
        split_target_counts={s: int(sum(v == s for v in target_split.values())) for _, s in SPLITS},
        control_observation_ids_sha256=axis_sha(sorted(control_fingerprints)),
        control_raw_fingerprints_sha256=hashlib.sha256(b"".join(
            control_fingerprints[c] for c in sorted(control_fingerprints)
        )).hexdigest(),
        threads=THREADS, chunk_size=chunk_size,
        definitions=dict(
            count_sums="Exact integer counts from X; each unique cell contributes once.",
            cpm="Per-cell raw counts / sum of all 18080 measured genes * 1e6; arithmetic cell mean.",
            batch_matching="Control probabilities use target library-mass batch weights; control CPM means use target cell-count batch weights.",
            target_variance_cpm="Unbiased between-cell sample variance (ddof=1).",
            matched_control_variance_cpm="Population variance of the batch-mixture control distribution.",
            matched_control_mean_variance_cpm="Variance of weighted batch means: sum(weight^2 * batch sample variance / batch n).",
            provenance="Source file path, size, mtime, axes fingerprints and audited count totals; source files are not separately whole-file hashed.",
        ),
    )
    receipt = output.with_suffix(".provenance.json")
    receipt.write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"output": str(output), "provenance": str(receipt),
                      "targets": nt, "genes": ng, "unique_controls": int(control_n.sum()),
                      "duplicate_controls_verified": duplicates_verified,
                      "elapsed_seconds": provenance["elapsed_seconds"]}, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=ROOT / "data/h1_2025")
    parser.add_argument("--output", type=Path, default=ROOT / "data/H1_2025_full_statistics.npz")
    parser.add_argument("--chunk-size", type=int, default=2048)
    args = parser.parse_args()
    if args.chunk_size <= 0:
        parser.error("--chunk-size must be positive")
    build(args.input_dir, args.output, args.chunk_size)


if __name__ == "__main__":
    main()
