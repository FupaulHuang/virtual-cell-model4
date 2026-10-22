# Frozen ABC control-only null protocol

**purpose**: Control-only technical null diagnostic. Aliases are not perturbations; never submit.

**contexts**: ["A", "B", "C"]

**seed**: 20261024

**aliases**: ["NULL_TECH_01", "NULL_TECH_02", "NULL_TECH_03", "NULL_TECH_04", "NULL_TECH_05", "NULL_TECH_06", "NULL_TECH_07", "NULL_TECH_08", "NULL_TECH_09", "NULL_TECH_10"]

**cells per alias**: 400

**comparison cells**: 4000

**donor population cells**: 4000

**roles**: Comparison4000, measured-null4000, donor-population4000: mutually disjoint within context.

**ntc stratification**: Fixed seeded 400-cell integer guide quota; identical per-alias quotas and comparison proportions. Each role has10times that quota.

**pooling**: Sample4times each alias guide quota without replacement from donor population; depth-sort WITHIN guide; groups of4. Replay chooses1random raw cell per pool; pool4 averages normalized cells and uses original rounding/repair without moment fitting.

**difference from h1 null**: Pools stay within ntc_id here, unlike prior H1 global depth pooling; this avoids a guide-mixing confound.

**reuse**: Donor cells may recur across technical aliases; aliases and shared-comparison tests are correlated, not independent biological replicates.

**de**: {"backend": "scanpy", "input_type": "counts", "target_sum": 1000000.0, "mean_calc": "arithmetic", "epsilon": 1e-09, "clip_value": null, "filter_gene_min_cpm_cell": 5.0, "fdr_scope": "per_pert", "threads": 1, "groupby": "target_gene", "reference": "non-targeting"}

**calls**: Strict p_adj<0.05; reference-only arithmetic CPM>5 gate and recomputed per-alias BH are provided by the pinned compute_de helper. Aliases do not match genes, so own-target exclusion removes nothing.

**distribution**: Full-row-depth CPM; sample variance ddof1; detection=positive fraction. Gene summaries on comparison-control meanCPM>5, positive-variance subset for log ratios. Depth KS is descriptive.

**identity gate**: fit_population at donor population own arithmetic mean/raw bulk, then emit the exact replay donor rows. Require expected probabilities and every emitted count to match raw donors exactly; no DE rerun if identity is proven.

**identity status at freeze**: Pending implementation/API availability; measured/replay/pool4 diagnostic is independently defined.

**compute bound**: Nine DE jobs, each10aliases versus4000 controls (90 contrasts); at most3 processes concurrently; no scorer or parameter sweep.

**no hidden or reserved response access**: true

Source identities, all selected row identities, code/helper hashes and package versions are recorded in `control_preserving_null_freeze.json` and `evaluations/control_preserving_null_abc/`. No expression outcomes were inspected before the metadata freeze.
