# Data setup

No data is distributed with this source release. Use the official VCC CLI to download the controls appropriate to the phase, preserving the supplied context labels, gene list, target list and manifest. Final controls use D/E/F; validation controls use A/B/C.

The final example configuration expects the unmodified final bundle under `data/final_20261022/`. Validation examples expect controls and their CSVs under `data/`. Public H1 2025 files belong under `data/h1_2025/`.

Source-statistics bank paths and weights are explicit in each `configs/*.json` file. The acquisition and aggregation implementations include `build_h1_bank.py`, `build_orion_raw_bank.py`, `build_orion_de_bank.py`, `gwps_subset_bank.py` and `build_cd4_ranged_bank.py`. Read their arguments and provenance checks before use; some acquisitions require separately downloaded manifests or metadata. Nothing in this repository bypasses source access requirements.

Use `python scripts/run_prediction_config.py CONFIG --dry-run` to review the exact required inputs and generation command. The runner rejects missing banks rather than substituting a different model. Source-target coverage must be checked against the chosen panel: a complete output shape does not mean that every target has a measured source response.

Historical experiment scripts also depend on frozen cohorts, reference DE tables and manifests that are intentionally omitted. Reconstruct and prospectively freeze a new experiment rather than fabricating missing historical receipts. Optional legacy comparator artifacts should be placed under `data/legacy_model3/`.

Generated `data/`, `outputs/`, `evaluations/`, `models/`, `submissions/`, `scratch/` and `logs/` directories are excluded from Git. Never add authentication files or signed download URLs to source control.
