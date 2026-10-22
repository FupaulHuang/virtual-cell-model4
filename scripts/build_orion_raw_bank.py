#!/usr/bin/env python3
"""Stream locally downloaded X-Atlas Parquet batches into matched count/CPM banks.

Files are accepted only at published size and SHA256. The script can run while
curl downloads the files; checkpoints make each validated batch restartable.
Only officially passing dual-guide cells and the requested target union enter.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

for _name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[_name] = "2"

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy import sparse

from build_orion_de_bank import ROOT, axes, sha256


def reduce_file(path, targets, genes, metadata):
    labels = list(targets) + ["Non-Targeting"]
    target_index = pd.Index(labels)
    token_count = int(metadata.gene_token_id.max()) + 1
    mapping = np.full(token_count, -1, np.int32)
    gene_index = pd.Index(genes)
    mapped = gene_index.get_indexer(metadata.gene_name.astype(str))
    # Distinct Ensembl loci can share a symbol; counts are summed onto the
    # challenge's symbol axis by CSR duplicate-coordinate accumulation below.
    good = mapped >= 0
    mapping[metadata.gene_token_id.to_numpy()[good]] = mapped[good]
    sums = np.zeros((len(labels), len(genes)), np.float64)
    cpm_sums = np.zeros_like(sums)
    n = np.zeros(len(labels), np.int64)
    raw_total = np.zeros(len(labels))
    samples = set()
    observed = 0
    columns = ["gene_token_id", "gene_expression", "gene_target", "pass_guide_filter", "total_counts", "sample"]
    for batch in pq.ParquetFile(path).iter_batches(batch_size=2048, columns=columns):
        observed += batch.num_rows
        tab = pa.Table.from_batches([batch])
        mask = pc.and_(pc.equal(tab["pass_guide_filter"], 1), pc.is_in(tab["gene_target"], value_set=pa.array(labels)))
        tab = tab.filter(mask)
        if not len(tab):
            continue
        samples.update(tab["sample"].to_pylist())
        tok = tab["gene_token_id"].combine_chunks()
        val = tab["gene_expression"].combine_chunks()
        offsets = tok.offsets.to_numpy().astype(np.int64)
        # Arrow list slices may retain offsets; normalize before constructing CSR.
        start, end = int(offsets[0]), int(offsets[-1])
        indices = tok.values.slice(start, end-start).to_numpy().astype(np.int32)
        values = val.values.slice(start, end-start).to_numpy().astype(np.float64)
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError("Non-count raw values")
        x = sparse.csr_matrix((values, indices, offsets-start), shape=(len(tab), token_count))
        depths = np.asarray(x.sum(axis=1)).ravel()
        metadata_depths = tab["total_counts"].to_numpy()
        if not np.allclose(depths, metadata_depths, rtol=0, atol=1e-6) or np.any(depths <= 0):
            raise ValueError("Raw count sums disagree with source library-depth metadata")
        gi = mapping[x.indices]
        xx = x.tocoo()
        use = gi >= 0
        aligned = sparse.csr_matrix((xx.data[use], (xx.row[use], gi[use])), shape=(len(tab),len(genes)))
        groups = target_index.get_indexer(tab["gene_target"].to_pylist())
        membership = sparse.csr_matrix((np.ones(len(groups)), (groups, np.arange(len(groups)))), shape=(len(labels),len(groups)))
        sums += (membership @ aligned).toarray()
        cpm_sums += (membership @ aligned.multiply(1e6/depths[:, None])).toarray()
        n += np.bincount(groups, minlength=len(labels))
        raw_total += np.bincount(groups, weights=depths, minlength=len(labels))
    if len(samples) != 1:
        raise ValueError(f"Expected a single source sample per Parquet batch, got {samples}")
    if not n[-1] or not sums[-1].sum():
        raise ValueError("Batch lacks usable matched controls")
    return dict(count_sums=sums, cpm_sums=cpm_sums, n_cells=n, raw_total=raw_total,
                sample=np.asarray(next(iter(samples))), input_cells=np.asarray(observed))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line", choices=["HCT116", "HEK293T"], required=True)
    parser.add_argument("--wait-seconds", type=int, default=7200)
    args = parser.parse_args()
    started = time.monotonic()
    targets, genes, official, h1_targets = axes()
    listing = json.loads((ROOT/"research/orion_raw_listing.json").read_text())
    records = sorted([x for x in listing if Path(x["path"]).name.startswith(args.line+"_Batch")],
                     key=lambda x:int(x["path"].split("_Batch")[1].split(".")[0]))
    metadata = pd.read_parquet(ROOT/"data/raw/orion_gene_metadata.parquet")
    tmpdir = ROOT/"data/processed/orion_batches"/args.line
    tmpdir.mkdir(parents=True, exist_ok=True)
    n_targets, n_genes = len(targets), len(genes)
    sums = np.zeros((n_targets, n_genes))
    cpms = np.zeros_like(sums)
    matched_cp = np.zeros_like(sums)
    matched_cm = np.zeros_like(sums)
    n_cells = np.zeros(n_targets, np.int64)
    ctrl_sums = np.zeros(n_genes)
    ctrl_cpms = np.zeros(n_genes)
    ctrl_n = 0
    manifest = []
    for i, entry in enumerate(records):
        path = ROOT/"data/raw/orion"/Path(entry["path"]).name
        cache = tmpdir/(path.stem+".npz")
        receipt = cache.with_suffix(".json")
        if cache.exists() and receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved["source_sha256"] != entry["lfs"]["oid"] or saved["axis_sha256"] != hashlib.sha256(("\0".join(targets)+"\n"+"\0".join(genes)).encode()).hexdigest():
                raise ValueError("Stale batch cache")
            with np.load(cache, allow_pickle=False) as z:
                result = dict(z)
        else:
            waiting = time.monotonic()
            last_message = 0
            while not path.exists() or path.stat().st_size < entry["size"]:
                if time.monotonic()-waiting > args.wait_seconds:
                    raise TimeoutError(f"Download did not complete: {path}")
                if time.monotonic()-last_message > 60:
                    print(f"{args.line}: waiting for {path.name}", flush=True)
                    last_message = time.monotonic()
                time.sleep(3)
            if path.stat().st_size != entry["size"] or sha256(path) != entry["lfs"]["oid"]:
                raise ValueError(f"Source size/SHA256 mismatch: {path}")
            result = reduce_file(path, targets, genes, metadata)
            np.savez_compressed(cache, **result)
            receipt.write_text(json.dumps(dict(source_sha256=entry["lfs"]["oid"],
                axis_sha256=hashlib.sha256(("\0".join(targets)+"\n"+"\0".join(genes)).encode()).hexdigest(),
                input=entry["path"], bytes=entry["size"]), indent=2)+"\n")
        counts, mean_sum, n = result["count_sums"], result["cpm_sums"], result["n_cells"]
        cp = counts[-1]/counts[-1].sum()
        cm = mean_sum[-1]/n[-1]
        sums += counts[:-1]
        cpms += mean_sum[:-1]
        n_cells += n[:-1]
        matched_cp += counts[:-1].sum(axis=1)[:, None]*cp[None,:]
        matched_cm += n[:-1,None]*cm[None,:]
        ctrl_sums += counts[-1]
        ctrl_cpms += mean_sum[-1]
        ctrl_n += int(n[-1])
        manifest.append(dict(path=entry["path"],sha256=entry["lfs"]["oid"],
            input_cells=int(result["input_cells"]),selected_cells=int(n.sum()),controls=int(n[-1])))
        print(f"{args.line}: {i+1}/{len(records)} {path.name}; {n_cells.sum()} treated, {ctrl_n} controls",flush=True)
    keep = n_cells > 0
    probability = np.divide(matched_cp, sums.sum(axis=1)[:,None], out=np.zeros_like(matched_cp), where=sums.sum(axis=1)[:,None]>0)
    cpm = np.divide(cpms,n_cells[:,None],out=np.zeros_like(cpms),where=n_cells[:,None]>0)
    matched_cpm = np.divide(matched_cm,n_cells[:,None],out=np.zeros_like(matched_cm),where=n_cells[:,None]>0)
    measured = np.isin(np.arange(n_genes),pd.Index(genes).get_indexer(metadata.gene_name.astype(str)))
    destination = ROOT/f"data/processed/{args.line}_full_statistics.npz"
    np.savez_compressed(destination, source=np.asarray(f"Xaira_Orion_{args.line}"), targets=targets[keep],genes=genes,
        target_count_sums=sums[keep].astype(np.float32), target_mean_cpm=cpm[keep].astype(np.float32),n_cells=n_cells[keep],
        matched_control_probability=probability[keep].astype(np.float32),matched_control_mean_cpm=matched_cpm[keep].astype(np.float32),
        global_control_probability=(ctrl_sums/ctrl_sums.sum()).astype(np.float32),
        global_control_mean_cpm=(ctrl_cpms/ctrl_n).astype(np.float32),measured_genes=measured)
    present=set(targets[keep])
    eligible=set(targets[n_cells>=20])
    duplicates = metadata[metadata.gene_name.isin(set(genes)) & metadata.gene_name.duplicated(keep=False)]
    report=dict(source_url="https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion",
        license="CC-BY-NC-SA-4.0",output=str(destination),targets=int(keep.sum()),genes=n_genes,
        official_coverage=len(set(official)&present),h1_coverage=len((set(h1_targets)-{"non-targeting"})&present),
        official_coverage_min20=len(set(official)&eligible),h1_coverage_min20=len((set(h1_targets)-{"non-targeting"})&eligible),
        official_targets_missing=sorted(set(official)-present),
        official_targets_below20=sorted(set(official)-eligible),
        n_cells_total=int(n_cells.sum()),n_controls=ctrl_n,median_cells=float(np.median(n_cells[keep])),
        minimum_cells=int(n_cells[keep].min()),measured_genes=int(measured.sum()),
        source_filter="pass_guide_filter == 1 (two guides targeting same gene); no target-specific efficacy filter",
        normalization="Arithmetic CPM denominator uses full source library; pseudobulk uses raw summed counts on saved gene axis",
        control_matching="One matched control population per source sample/batch; bulk weighted by treated library mass; CPM weighted by treated cell count",
        duplicate_symbol_tokens_summed=duplicates.to_dict(orient="records"),
        center_scope="The 575-target union of official validation and public H1 targets; raw banks are not centered",
        seconds=time.monotonic()-started,files=manifest)
    (ROOT/f"reports/{args.line}_full_statistics.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k!="files"},indent=2),flush=True)


if __name__ == "__main__":
    main()
