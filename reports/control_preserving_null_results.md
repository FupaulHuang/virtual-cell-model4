# Control-preserving null audit on supplied A/B/C controls

Each context uses disjoint 4,000-cell comparison, measured-null and donor populations, with exact ntc_id proportions. There are ten technical aliases of 400 cells per arm. All source cells are non-targeting controls; no hidden perturbed responses or reserved targets were accessed.

DE uses the pinned official scanpy helper: CPM normalization, arithmetic LFC, reference-control CPM>5 gene gate, tie-corrected Wilcoxon and per-alias BH; calls require q<0.05. Aliases are absent from the gene axis. These are diagnostic call counts, not a competition score.

| Context | Arm | Significant calls median (range) | Median log2 variance ratio | Detection MAE |
| --- | --- | ---: | ---: | ---: |
| A | measured_null | 0 (0–1) | -0.0172 | 0.0163 |
| A | donor_replay | 0 (0–0) | 0.0018 | 0.0160 |
| A | pool4_no_fit | 6004 (5761–6431) | -1.2100 | 0.1035 |
| B | measured_null | 0 (0–0) | -0.2059 | 0.0166 |
| B | donor_replay | 0 (0–1) | -0.2593 | 0.0171 |
| B | pool4_no_fit | 5687 (5294–6126) | -1.4514 | 0.0995 |
| C | measured_null | 0 (0–1) | -0.0307 | 0.0165 |
| C | donor_replay | 0 (0–2) | -0.0232 | 0.0161 |
| C | pool4_no_fit | 5191 (4715–5506) | -1.1777 | 0.0973 |

Variance and detection compare each group with its disjoint comparison controls on the same comparison-mean-CPM>5 gate. The displayed variance statistic averages the per-alias median log2 variance ratio; detailed depth and gene-level records are retained.

New-emitter zero-effect identity: PASS in all three contexts, including every count in all30replay groups. Reusing donor-replay DE is justified by exact count identity.

Pools contain four cells from the same ntc_id and are depth-sorted within guide. This removes guide-mixing as an explanation but differs deliberately from the older H1 global-depth pooling experiment. Shared donors/comparison controls make aliases correlated; they are not independent biological replicates. Passing this null test alone does not establish a better perturbation model.

Full protocol, frozen hashes, group results and identity receipts: JSON report (local artifact omitted: `control_preserving_null_results.json`).
