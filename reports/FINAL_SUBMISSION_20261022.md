# Final submission — October 22, 2026

Final entry **`[PRIVATE_IDENTIFIER_REMOVED]`** for **`m4_multisource_dual_v1_final_20261022`** is **published** on partition **test**: upload and scoring completed successfully. The upload size and MD5 match the local package. The server reports no error and has released the in-flight slot. All scores and rank remain hidden, as required in the final round.

This is the unchanged latest validated `m4_multisource_dual_v1` model applied to the newly released final controls D/E/F and the official final perturbation list. The experimental emitter has not been promoted. The frozen source banks cover **141/300 final targets**; the remaining **159 use the existing zero-response fallback**. This limitation is recorded in the submission description and prediction receipt. No final-score improvement is claimed.

| Item | Verified value |
| --- | --- |
| Final panel | `vcc2026-test-1` |
| Contexts | D, E, F, preserved from official control files |
| Prediction | 360,000 cells × 18,533 genes |
| Groups | 900, each with exactly 400 cells |
| Stored entries | 2,341,872,450 |
| Maximum cell library | 53,290 counts |
| Package version | Official released VCC CLI 0.3.0 |
| Package bytes | 3,831,398,400 |
| Package SHA-256 | `f853504d2db858a8e790fdf0611c443fff55aee81357263455d6fe983659d7f6` |
| New entries created | One |
| Terminal status | `published`, partition `test` |
| Scoring anchor | `vcc2026-testD-r1+vcc2026-testE-r1+vcc2026-testF-r1` |
| Daily allowance after completion | 1 of 2 used; no entry in flight |

All source/model hashes matched the validated release. The official final-control archive passed its CRC32C check. The independent validator inspected every count block and all group/axis constraints. Official `vcc prep` then passed with raw counts preserved, no dropped data, and target/count/control/density checks enabled. The package metadata explicitly records D/E/F and the final panel.

The global CLI remains 0.1.0 because its installation is outside writable roots. A checksum-verified official 0.3.0 wheel is installed and used in `.tools/vcc-cli-0.3.0`, with a refreshed VCC skill (local artifact omitted: `../.tools/vcc-skill-0.3.0/vcc/SKILL.md`). Python DNS was unavailable, so the transfer used the working HTTPS transport to the same official production API and its returned storage session. Current account eligibility, phase, deadline, quota, in-flight state, released CLI version and package-context checks were preserved. No endpoint override or quota bypass occurred. The API was called with the official 0.3.0 user agent. Credentials and signed URLs are excluded from these receipts.

For future checks from a terminal with network access, the [workspace CLI wrapper](../scripts/vcc_final.sh) runs the verified release with the existing login:

```bash
bash scripts/vcc_final.sh status [PRIVATE_IDENTIFIER_REMOVED] --json
```

Final scoring is blind: the server can confirm successful scoring while keeping all metric and rank fields null. This is expected under the [official final-phase instructions](https://vcc-cli-wiki.virtualcellchallenge.org/). The validation score is not a prediction of this final panel's score. Only the latest final entry counts for prizes; no second entry is authorized by this execution.

Temporary transport credentials, signed URLs and the local lock were removed after terminal verification. The original credential store and frozen model remain unchanged. Cleanup receipt (local artifact omitted: `final_20261022_transport_cleanup.json`) · Final release record (local artifact omitted: `../models/m4_multisource_dual/final_release_20261022.json`).

Evidence: preflight (local artifact omitted: `final_20261022_preflight.json`), [configuration](../configs/def_final_20261022.json), prediction receipt (local artifact omitted: `../outputs/m4_multisource_dual_v1_final_20261022/prediction.receipt.json`), full validation (local artifact omitted: `../outputs/m4_multisource_dual_v1_final_20261022/validation.json`), official prep (local artifact omitted: `final_20261022_prep.json`), realized count profiles and uncovered targets (local artifact omitted: `final_20261022_realized_profiles.json`), package manifest (local artifact omitted: `../submissions/m4_multisource_dual_v1_final_20261022.manifest.json`), entry (local artifact omitted: `../submissions/m4_multisource_dual_v1_final_20261022_entry.json`), upload verification (local artifact omitted: `../submissions/m4_multisource_dual_v1_final_20261022_upload.json`), launch (local artifact omitted: `../submissions/m4_multisource_dual_v1_final_20261022_launch.json`), latest status (local artifact omitted: `../submissions/m4_multisource_dual_v1_final_20261022_status.json`), execution state (local artifact omitted: `final_20261022_state.json`).
