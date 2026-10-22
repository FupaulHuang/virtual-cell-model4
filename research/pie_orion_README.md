---
pretty_name: PIE X-Atlas/Orion genome-wide Perturb-seq (HCT116, HEK293T)
tags:
  - biology
  - single-cell
  - perturb-seq
  - crispri
  - differential-expression
size_categories:
  - 100M<n<1B
configs:
  - config_name: hct116
    data_files: de/hct116.parquet
  - config_name: hek293t
    data_files: de/hek293t.parquet
---

# PIE: X-Atlas/Orion

Differential-expression (DE) tables and the preprocessed dir for
[PIE](https://github.com/ArcInstitute/pie), built from X-Atlas/Orion (Huang et al., 2025),
genome-wide CRISPRi Perturb-seq screens in HCT116 and HEK293T. This dataset is a training dataset of
the PIE `replogle_xdataset` experiment.

Datasets and knowledge sources: [PIE collection](https://huggingface.co/collections/arcinstitute/pie-6abc354d475081a4a7cb4492)

## Contents

```
de/<cell_line>.parquet   per-cell-line DE, one row per (perturbation, gene)
preprocessed/            PIE preprocessed dir, with the context map contexts.yaml
```

| Cell line | Cellosaurus | Cells     | Control cells | Perturbed genes | DE rows     |
|-----------|-------------|-----------|---------------|-----------------|-------------|
| HCT116    | CVCL_0291   | 1,684,046 | 165,777       | 8,987           | 102,228,672 |
| HEK293T   | CVCL_0063   | 1,182,507 | 218,838       | 4,297           | 53,646,719  |

Controls are the `non-targeting` cells of each cell line.

### `de/`

One parquet per cell line, from the `pie process` `de` stage (gpudge 0.7.0 on raw counts,
grouped by `gene` against `non-targeting`, CP10k normalization, genes kept at ≥ 5 CPM).

| Column             | Description                                        |
|--------------------|----------------------------------------------------|
| `target`           | perturbed gene                                     |
| `feature`          | measured gene                                      |
| `target_mean`      | mean expression in perturbed cells                 |
| `ref_mean`         | mean expression in control cells                   |
| `target_ncells`    | number of perturbed cells                          |
| `ref_ncells`       | number of control cells                            |
| `log2_fold_change` | log2 fold change, perturbed vs. control            |
| `p_value`          | Mann-Whitney U test p-value                        |
| `Ueffect`          | U-statistic effect size                            |
| `p_adj`            | Benjamini-Hochberg adjusted p-value                |

### `preprocessed/`

The PIE preprocessed dir (`pie prep` output, format version 1) that `pie sources`, `pie train`,
`pie eval` and `pie infer` read: N = 13,284 rows, one per (cell line, perturbation) pair, over
G = 22,182 genes and C = 2 cell lines.

| File               | Shape, dtype      | Description                                                    |
|--------------------|-------------------|----------------------------------------------------------------|
| `meta.json`        |                   | gene axis, cell line and perturbation vocabularies, checksums |
| `contexts.yaml`    |                   | context -> Cellosaurus accession map; `pie sources` reads it for `context_text` |
| `fold_changes.npy` | (N, G) float32    | linear fold change, 0 where untested                           |
| `fdr.npy`          | (N, G) float32    | adjusted p-value, 1 where untested                             |
| `tested.npy`       | (N, G) bool       | whether the gene was tested for the pair                       |
| `lfc_true.npy`     | (N, G) float64    | log2 fold change, NaN where untested                           |
| `delta_p.npy`      | (N, G) float32    | pseudobulk log1p expression of perturbed cells minus control mean |
| `ctrl_means.npy`   | (C, G) float32    | mean log1p expression of control cells per cell line           |
| `ctx_ids.npy`      | (N,) int32        | cell line of each row, index into `context_to_id`              |
| `pert_ids.npy`     | (N,) int32        | perturbation of each row, index into `pert_to_id`              |

## Usage

The `replogle_xdataset` experiment config pins this repo by commit, so `pie train` and `pie eval`
download `preprocessed/` automatically. In your own config, use
`hf://datasets/arcinstitute/PIE_x_atlas_orion@<commit>/preprocessed` in `data.preprocessed_dirs`, or
download a copy:

```bash
hf download arcinstitute/PIE_x_atlas_orion --repo-type dataset \
  --local-dir "$PIE_DATA_ROOT/orion"
```
