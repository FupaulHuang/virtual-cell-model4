# Bulk ablation implementation evidence

The isolated builder passed ten tests and four engineering runs on the same two targets, ACVR1B and ATM (400 cells each, 18,080 genes). Every output passed the unchanged strict validator, including every stored count, axes/group sizes, H1 source exclusion, and file SHA-256. No checks were disabled. No full experiment arm or metric job was run.

The new incumbent reproduces the original saved integer counts exactly for both targets. All four arms have identical requested arithmetic-profile, donor-template, row-depth, and centered-CD4 hashes. The variants have different realized arithmetic CPM after fitting and rounding; this is not a claim of exact or 0.01%-bounded CPM preservation.

| Arm | Incumbent counts identical | Maximum projection L1 | Maximum realized gene-CPM drift | Validation |
| --- | --- | ---: | ---: | --- |
| incumbent | True | 0.000098682 | 0.000000 | PASS |
| common05 | False | 0.000097013 | 0.840178 | PASS |
| bulk015 | False | 0.000036559 | 0.849986 | PASS |
| bulk060 | False | 0.000273300 | 0.809282 | PASS |

The original projection L1 limit remains 0.03 and the fitting L1 limit remains 0.001. Requested, projected, pre-round fitted, and final bulk probabilities are saved per gene. The actual frozen emitter supplies the projection diagnostics through a read-only return-frame profiler; its numerical implementation is unchanged, and the previous profiler is restored on success or failure. A synthetic infeasible example retains the original failure threshold.

The ten tests cover incumbent profile/count equivalence, fixed arithmetic/CD4 calculations, full-source centering and missingness, raw-only centering contrast, amplitude-only contrast, actual emitter diagnostics, equal-depth projection, failure/hook restoration, and independent post-round margin/CPM checks. A first attempt with base Python could not import anndata; the successful test command uses the existing cell-eval2 environment.

```text
.venv/bin/python -m unittest discover -s tests -p test_bulk_ablation.py -v
```

Builder SHA-256: `2819d37693e0b0c6a0365644ef441aeb2c855b8a18b65576ed31d01a848d5db6`.

Test SHA-256: `f289c96cd17c823700781695dac552b14705e071152ac2803a96ef91148bb347`.

Each engineering output is under `outputs/bulk_engineering_{incumbent,common05,bulk015,bulk060}/`. Each has `prediction.h5ad`, `prediction.receipt.json`, `profile_audit.json`, and `gene_profile_audit.parquet`. Exact executed generation/validation commands, output/receipt/audit hashes, profile hashes, per-target diagnostics and validation results are in the JSON evidence (local artifact omitted: `bulk_implementation_evidence.json`).

Independent reviewer approval and root freeze are still required before full generation. Full 50-target count identity must be independently checked before reusing incumbent metrics. Engineering success does not establish a score improvement, unchanged DE metrics, or eligibility for promotion.
