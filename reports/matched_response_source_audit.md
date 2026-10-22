# Matched source-response provenance audit

**No coherent matched DE-label set is available for the three current raw banks.** Use the frozen raw effects for a response-only benchmark, with predeclared magnitude labels explicitly distinguished from DE significance.

| Source | Raw targets with cells / eligible ≥20 | Raw gene axis | Candidate DE targets | Matched raw-population DE labels? |
| --- | ---: | ---: | ---: | --- |
| K562 | 527 / 522 | 8246 | absent | No |
| HCT116 | 574 / 554 | 18536 | 263 | No |
| HEK293T | 574 / 561 | 18536 | 27 | No |

K562 contains genuine selected raw counts and separate arithmetic mean CPM/mean-variance estimates. It samples up to200 cells per target and50 controls per batch. There is no matching p-value table in this bank; means or standard errors must not be relabeled as Wilcoxon/BH confidence.

**HCT116:** 249/263 shared targets have different perturbation cell counts; median DE/raw count ratio 0.8986. On shared, tested, own-excluded pairs, median absolute log2(CPM+1) mean discrepancy is 0.031513 for targets and 0.009207 for matched controls.

**HEK293T:** 20/27 shared targets have different perturbation cell counts; median DE/raw count ratio 0.8910. On shared, tested, own-excluded pairs, median absolute log2(CPM+1) mean discrepancy is 0.017977 for targets and 0.008968 for matched controls.

The raw Orion banks retain guide-passing cells without an efficacy filter and use target-weighted per-batch controls. Published PIE DE is a separately filtered view, with per-cell-line non-targeting controls. The local report attributes strong knockdown filtering; vendored defaults use strict group residual<0.30, cell residual<0.50 and at least30 cells, but exact release cell memberships and override commands are not recorded here. The observed cell-count and matched-control discrepancies independently prevent treating these q-values as labels for our raw-bank responses.

PIE tables contain adjusted p-values and an explicit tested mask; missing records are not measured zeros or negative calls. Their arithmetic CP10k means can be converted to CPM, but the stored `effects_bulk_delta` is an arithmetic approximation and does not establish measured raw-count pseudobulk. Matching family names and target symbols does not resolve these population, control or estimand differences.

For future true DE-confidence labels, build raw effects and DE jointly from identical cells, controls and normalization/test definitions. A batch-weighted control mean may require a batch-aware test, not a global-control q-value. Treat raw/DE views as one correlated source family.

No H1 or reserved response measurements were accessed. This audit does not fit a model or alter frozen files. Full identities, comparison counts, settings and hashes are in JSON (local artifact omitted: `matched_response_source_audit.json`).

Evidence: [raw Orion builder](../scripts/build_orion_raw_bank.py), [PIE conversion](../scripts/build_orion_de_bank.py), [K562 builder](../scripts/gwps_subset_bank.py), [archived publisher README](../research/pie_orion_README.md), PIE Orion configuration (local artifact omitted: `../external/pie/src/pie/configs/process/dataset/orion.yaml`).
