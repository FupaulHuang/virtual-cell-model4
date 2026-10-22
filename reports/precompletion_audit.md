Pre-completion audit rechecked at `2026-10-20T15:58:28.131869+00:00`. The official receipt now verifies **PUBLISHED**, score **0.1470186853**, rank **565**. The same-panel/anchor improvement over Model3 is **+0.4416689922**. All six official scaled components are present and their mean matches the overall score.

All inspected lightweight model/config/core/artifact/upload integrity checks passed. Full prediction count scans and whole prediction/package hashes were not rerun; their existing exhaustive receipts, current sizes/mtime, and hash chains were checked.

| Requirement | Evidence/result |
| --- | --- |
| Prior-work and current-competition summary | present; current status prose needs refresh; `docs/RESEARCH_SUMMARY.md`, `research/current_methods_20261020.md`, `reports/leaderboard_20261020_late.json`, `reports/leaderboard_late_summary.json` |
| Three explanatory implementation/experiment/submission plans | present; execution-status checklists lag; `PLAN.md`, `EXPERIMENT_PLAN.md`, `SUBMISSION_PLAN.md` |
| Reproducible selected model, configuration, code and source bank identity | verified by lightweight hash/path audit; `models/m4_multisource_dual/model.json`, `models/m4_multisource_dual/environment.json`, `configs/abc_fullsource_candidate.json`, `models/m4_multisource_dual/release.json` |
| Twelve candidates on one fresh five-split CPU reference | verified; `reports/h1_all_candidate_comparison.json`, `models/m4_multisource_dual/selection_crosscheck.json`, `evaluations/h1_validation/real_bundle_cpu/manifest.json` |
| Both reserved confirmations after model freeze | verified; both improve all six raw metrics; `reports/h1_test_confirmation.json`, `reports/h1_test_terminal_statuses.json`, `reports/k562_confirmation.json`, `reports/k562_confirmation_job_state.json` |
| Complete production artifact and strict official prep | verified from prior exhaustive receipts plus current size/hash-chain checks; `outputs/m4_fullsource_candidate/prediction.receipt.json`, `outputs/m4_fullsource_candidate/validation.json`, `outputs/m4_fullsource_candidate/csr_pointer_audit.json`, `reports/m4_fullsource_candidate_prep.json`, `submissions/m4_fullsource_candidate.manifest.json` |
| Authorized submission and verified upload/launch | verified; `submissions/m4_multisource_dual_v1_entry.json`, `submissions/m4_multisource_dual_v1_upload.json`, `submissions/m4_multisource_dual_v1_launch.json`, `reports/submission_transport.json` |
| Persisted published result and official comparison with prior best | verified published receipt and same-scale score comparison; `submissions/m4_multisource_dual_v1_status.json` |

Findings observed while parent updates were in progress:

- **A2 (operational; resolved_by_parent_during_audit): Continuation record mixes current accepted scoring state with obsolete active instructions.** Keep current outcome and remaining actions authoritative; move obsolete acquisition/evaluation/upload statements to history. Remove instructions to submit again, finish the old PUT, or await already completed H1 confirmation. Do not restart jobs from historical handles. Evidence: `reports/continuation_state.json`.
- **A3 (documentation; resolved_by_parent_during_audit): Canonical progress prose and plan gates lag completed work.** Update README and the three plans to actual completed states and terminal outcome. Preserve honest local-versus-official score distinctions; the final outcome report should link both reserved confirmation reports and the published receipt. Evidence: `EXPERIMENT_PLAN.md: Pending complete raw banks`, `SUBMISSION_PLAN.md: - [ ] Complete `reports/MODEL_SELECTION.md``.
- **A4 (documentation; resolved_by_parent_during_audit): Selection record needs a current release/confirmation pointer without rewriting freeze history.** Retain the immutable model.json and historical selection wording. Append current links to release.json, both completed confirmations, and the final official outcome; do not change frozen model/code hashes. Evidence: `reports/MODEL_SELECTION.md`, `models/m4_multisource_dual/model.json`, `models/m4_multisource_dual/release.json`.
- **A5 (minor; open_or_parent_edit_in_progress): Current package manifest retains its prepromotion update timestamp.** Use a current status_updated_utc timestamp, or explicitly distinguish original preparation time from later confirmation/submission/publication events. Evidence: `submissions/m4_fullsource_candidate.manifest.json`.
- **A6 (navigation; resolved_by_parent_during_audit): README omits the actual bounded GWPS/CD4 acquisition path used by the selected model.** Add links to the selected source builders/provenance and to submission_transport.json. The listed full-file aggregation route alone does not reproduce the actual subset/range acquisition path. Evidence: `README.md`, `data/raw/gwps_subset/ACQUISITION.md`, `scripts/gwps_subset_bank.py`, `scripts/build_cd4_ranged_bank.py`.

Freeze snapshots are intentionally historical. Preserve `model.json`, both frozen confirmation snapshots, the prepromotion manifest, and entry-creation receipt. `release.json`, the current package manifest, and the official status receipt are the present release authority.

The filename-only scan covered **117 text files** in logs/reports/models/submissions and found **0 matching filenames** for recognizable PAT/bearer-token or upload-session URL patterns. No matched values were printed. The temporary private runtime directory is absent. Original credential preservation is documented by the transport receipt; this audit did not read that credential store.

Local Markdown links still unresolved at the latest check:

- `README.md` → `docs/NEXT_ITERATION_PLAN.md`
- `PLAN.md` → `docs/NEXT_ITERATION_PLAN.md`
- `EXPERIMENT_PLAN.md` → `docs/NEXT_ITERATION_PLAN.md`
- `SUBMISSION_PLAN.md` → `docs/NEXT_ITERATION_PLAN.md`
- `docs/RESEARCH_SUMMARY.md` → `NEXT_ITERATION_PLAN.md`

No submissions, status-job restarts, large checks, or frozen model/code changes were performed. This audit supplies evidence for the parent completion decision; it does not change the goal status.
