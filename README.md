# Virtual Cell Model 4

Source-only release of a cross-context perturbation-response model and its diagnostic experiments for the Virtual Cell Challenge 2026.

The validated model transfers measured responses from public K562, HCT116, HEK293T, H1 and CD4 sources to destination controls. Separate arithmetic-expression and bulk-expression heads feed an integer-count generator. Experimental work investigates control-preserving emission, donor sampling, response transfer and independent injection–recovery.

## Results and current limits

- The validated model's recorded validation score was **0.1470186853**. This is historical performance on the validation panel, not a final-test score.
- A final D/E/F prediction was successfully scored on October 22, 2026. Final scores are hidden. The frozen source banks covered **141/300 final targets**; the remaining 159 used the existing zero-response fallback.
- Balanced sampling reduced H1 bulk-profile sampling error by **4.1%**, below the fixed 5% advancement threshold.
- The experimental preserving emitter reduced zero-effect DE calls, but missed many injected upregulated effects. It was not promoted into the final submission.

Read the [cross-domain study](reports/CROSS_DOMAIN_FINAL_20261021.md), [next research plan](docs/DISTRIBUTION_RECOVERY_NEXT_PLAN_20261021.md), and [sanitized final-submission report](reports/FINAL_SUBMISSION_20261022.md).

## What is included

`scripts/` contains generation, validation, acquisition and experiment code. `tests/` contains small synthetic unit tests. `configs/` contains parameter and path configurations. `docs/`, the top-level plans, and Markdown reports preserve the research rationale and aggregate findings. `external/atlasshift/` contains the MIT-licensed numerical kernels used by the model, with their original license and provenance.

Datasets, source-statistics banks, trained weights, predictions, submission packages, account identifiers, raw API receipts, logs, caches and credentials are **not included**. This is a clean source snapshot with no imported private Git history. Local artifact references in historical reports are explicitly marked when omitted.

## Setup

Use Python 3.12 or newer and an isolated environment:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

The dependency file records the versions used by the project. Full evaluation scripts use `cell-eval2` 0.16.0; dataset acquisition and challenge submission require access to the respective services. The included upstream AtlasShift requirements describe its original environment; the root requirements describe this project's environment.

Run the small synthetic tests without downloading challenge datasets:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m pytest tests -q
```

## Bring your own data

Obtain official controls and target/gene lists through the VCC website or CLI. Obtain public source data under its applicable terms and build the source banks using the acquisition scripts. See [data setup](docs/DATA_SETUP.md) and [GWPS acquisition](docs/GWPS_ACQUISITION.md).

Configurations use project-relative paths. The final example expects controls and metadata under `data/final_20261022/`; source banks are specified explicitly in the configuration. Put public H1 2025 inputs under `data/h1_2025/`. Historical Model 3 comparator artifacts, if used, belong under `data/legacy_model3/`.

After providing the inputs, inspect the command before generating predictions:

```bash
python scripts/run_prediction_config.py configs/def_final_20261022.json --dry-run
python scripts/run_prediction_config.py configs/def_final_20261022.json
python scripts/validate_prediction_receipt.py \
  outputs/m4_multisource_dual_v1_final_20261022/prediction.h5ad
```

This repository does not reconstruct omitted datasets or historical artifacts automatically. Experiment scripts referring to frozen local manifests require those manifests and their matching artifacts. Historical hashes describe the original experiment; privacy and portability edits mean this export is not a byte-identical forensic archive.

## Credentials and publication checks

Authenticate using the official VCC/GitHub tools in your own environment. Do not put tokens in code, configs, command arguments, notebooks or commits. `scripts/vcc_final.sh` uses an installed `vcc` executable, or `VCC_BIN` when explicitly set.

The repository uses an allowlist-based `.gitignore` and a publication scanner. Enable the optional local hook before adding changes:

```bash
git config core.hooksPath .githooks
python scripts/check_public_export.py
```

The hook scans staged contents, including paths, for data artifacts, credentials, private-key blocks, signed URLs and machine-specific home/mount paths. It is a guard, not a substitute for reviewing what you publish. See [publication scope](PUBLICATION.md).

Third-party code retains its original licensing. See [third-party notices](THIRD_PARTY_NOTICES.md). No new license is asserted for the original project code in this snapshot.
