#!/usr/bin/env python3
"""Convert Arc's public PIE Orion DE release to an honest masked response bank.

The published target/ref means are arithmetic CP10k. They do not determine the
library-weighted pseudobulk; effects_bulk_delta is explicitly an approximation.
Full-screen common responses are accumulated before selecting the target union.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT


def axes():
    controls = PARENT / "data"
    h1 = PARENT / "data/h1_2025"
    genes = list(dict.fromkeys(pd.read_csv(controls / "gene_names.csv").iloc[:, 0].astype(str).tolist()
                              + pd.read_csv(h1 / "gene_names.csv").iloc[:, 0].astype(str).tolist()))
    official = pd.read_csv(controls / "pert_counts.csv").target_gene.astype(str).tolist()
    h1_targets = []
    for split in ["Training", "Validation", "Test"]:
        h1_targets += pd.read_csv(h1 / f"pert_counts_{split}.csv").target_gene.astype(str).tolist()
    targets = list(dict.fromkeys(official + [t for t in h1_targets if t != "non-targeting"]))
    return np.asarray(targets), np.asarray(genes), official, h1_targets


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--line", required=True, choices=["hct116", "hek293t"])
    args = parser.parse_args()
    started = time.monotonic()
    targets, genes, official, h1_targets = axes()
    target_index, gene_index = pd.Index(targets), pd.Index(genes)
    path = ROOT / f"data/raw/pie_orion_{args.line}.parquet"
    listing = json.loads((ROOT / "research/pie_orion_de_listing.json").read_text())
    expected = next(x for x in listing if x["path"] == f"de/{args.line}.parquet")
    if path.stat().st_size != expected["size"]:
        raise ValueError("Download size mismatch")
    digest = sha256(path)
    if digest != expected["lfs"]["oid"]:
        raise ValueError("Download SHA256 mismatch")
    shape = (len(targets), len(genes))
    target_cpm = np.zeros(shape, np.float32)
    control_cpm = np.zeros(shape, np.float32)
    p_adj = np.ones(shape, np.float32)
    published_lfc = np.zeros(shape, np.float32)
    measured = np.zeros(shape, bool)
    n_cells = np.zeros(len(targets), np.int32)
    common_lfc_sum = np.zeros(len(genes))
    common_bulk_sum = np.zeros(len(genes))
    common_count = np.zeros(len(genes), np.int64)
    full_targets = set()
    file = pq.ParquetFile(path)
    for ri in range(file.num_row_groups):
        frame = file.read_row_group(ri, columns=["target", "feature", "target_mean", "ref_mean", "target_ncells", "p_adj", "log2_fold_change"]).to_pandas()
        full_targets.update(frame.target.unique())
        columns = gene_index.get_indexer(frame.feature)
        use = columns >= 0
        frame = frame.loc[use].reset_index(drop=True)
        columns = columns[use]
        pt = frame.target_mean.to_numpy() * 100.0
        ct = frame.ref_mean.to_numpy() * 100.0
        if not (np.isfinite(pt).all() and np.isfinite(ct).all() and (pt >= 0).all() and (ct >= 0).all()):
            raise ValueError("Invalid mean expression")
        lfc = np.log2((pt + 1.0) / (ct + 1.0))
        bulk = np.log1p(0.05 * pt) - np.log1p(0.05 * ct)
        nonself = frame.target.to_numpy() != frame.feature.to_numpy()
        common_lfc_sum += np.bincount(columns[nonself], weights=lfc[nonself], minlength=len(genes))
        common_bulk_sum += np.bincount(columns[nonself], weights=bulk[nonself], minlength=len(genes))
        common_count += np.bincount(columns[nonself], minlength=len(genes))
        rows = target_index.get_indexer(frame.target)
        selected = rows >= 0
        rr, cc = rows[selected], columns[selected]
        if measured[rr, cc].any():
            raise ValueError("Duplicate target-feature records")
        target_cpm[rr, cc] = pt[selected]
        control_cpm[rr, cc] = ct[selected]
        p_adj[rr, cc] = frame.p_adj.to_numpy()[selected]
        published_lfc[rr, cc] = frame.log2_fold_change.to_numpy()[selected]
        measured[rr, cc] = True
        cells = frame.target_ncells.to_numpy()[selected]
        previous = n_cells[rr]
        if np.any((previous > 0) & (previous != cells)):
            raise ValueError("Inconsistent cell count for target")
        n_cells[rr] = cells
        if (ri + 1) % 100 == 0:
            print(f"{args.line}: {ri+1}/{file.num_row_groups} row groups", flush=True)
    keep = n_cells > 0
    target_cpm, control_cpm = target_cpm[keep], control_cpm[keep]
    effects_lfc = np.log2((target_cpm + 1.0) / (control_cpm + 1.0))
    effects_bulk = np.log1p(0.05 * target_cpm) - np.log1p(0.05 * control_cpm)
    common_lfc = np.divide(common_lfc_sum, common_count, out=np.zeros_like(common_lfc_sum), where=common_count > 0)
    common_bulk = np.divide(common_bulk_sum, common_count, out=np.zeros_like(common_bulk_sum), where=common_count > 0)
    out = ROOT / "data/processed"
    out.mkdir(parents=True, exist_ok=True)
    destination = out / f"PIE_{args.line.upper()}_DE_statistics.npz"
    np.savez_compressed(destination, source=np.asarray(f"Arc_PIE_Orion_{args.line}_strong_KD"),
        targets=targets[keep], genes=genes, target_mean_cpm=target_cpm,
        control_mean_cpm=control_cpm, log2fc=published_lfc[keep],
        measured=measured[keep], n_cells=n_cells[keep], p_adj=p_adj[keep],
        effects_log2fc=effects_lfc, effects_bulk_delta=effects_bulk,
        source_common_log2fc=common_lfc.astype(np.float32),
        source_common_bulk_delta=common_bulk.astype(np.float32),
        source_common_n=common_count, effects_centered=np.asarray(False),
        bulk_is_approximation=np.asarray(True), source_full_target_count=np.asarray(len(full_targets)))
    present = set(targets[keep])
    report = dict(source_url=f"https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion/blob/main/de/{args.line}.parquet",
        release=json.loads((ROOT / "research/pie_orion_release.json").read_text())["sha"],
        input_path=str(path), input_sha256=digest, output=str(destination), full_source_targets=len(full_targets),
        selected_targets=int(keep.sum()), selected_gene_axis=len(genes), official_coverage=len(set(official)&present),
        h1_coverage=len((set(h1_targets)-{"non-targeting"})&present),
        measured_pairs=int(measured.sum()), cells_min=int(n_cells[keep].min()), cells_median=float(np.median(n_cells[keep])),
        effects="Uncentered log2((target_CPM+1)/(reference_CPM+1)); source_common excludes self genes over the full source screen",
        bulk_limit="Arithmetic-CPM approximation log1p(0.05*target_CPM)-log1p(0.05*reference_CPM), not observed library-weighted pseudobulk",
        missing="Per-target-feature measured mask; missing table records are not zero measurements",
        efficacy_filter="Arc PIE: group residual expression <=0.30, cell residual expression <0.50, >=30 cells",
        seconds=time.monotonic()-started)
    (ROOT / f"reports/PIE_{args.line.upper()}_DE_statistics.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
