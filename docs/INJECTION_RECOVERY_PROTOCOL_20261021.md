# Independent injection–recovery protocol

This is a synthetic/control-only diagnostic. It does not use actual H1 perturbed responses, the 25-target reserve, or official submission scores. Root and the independent referee approved the design before truth generation. No tested emitter, transport, repair or fitted prediction is used to generate reference counts.

This is specifically a **moment-oracle diagnostic**: candidates receive known arithmetic/bulk response profiles, not the target cell-state weights or complete marginal/joint distributions. Different distributions can share those moments. Poor DE recovery can therefore indicate insufficient two-moment information or an inadequate generator family, rather than a numerical implementation error. The oracle is not a full biological response predictor.

## Frozen scope and stages

There are two families, eight aliases per family and seeds **20261026/20261027**. Each family has two templates and four strengths: NULL, WEAK, MODERATE, STRONG. Each reference file contains 400 cells per alias and 4,000 comparison controls. Each candidate receives 4,000 inference controls. Template 1/seed 20261026 is development; template 2 and seed 20261027 are held out from emitter-development decisions until the prospectively declared confirmation stage. Preparing and auditing these files does not authorize adaptive use of their recovery outcomes.

The same frozen candidates and requests are compared: current emitter and preserving/uniform emitter. A preserving/balanced arm is eligible only if its separate precision and null requirements pass in all four contexts; eligibility is recorded before opening injection outcomes. No arm is added because it looks favorable here. The operational runner uses the pinned official DE helper only, not an official six-metric score or biological target resolver.

## Real-control composition family

Use only the existing 4,000 H1 inference controls, 4,000 old null comparison controls and 4,000 old measured-null controls as the reference pool. Verify their original IDs are mutually disjoint and every source cell is non-targeting. Preserve the complete 18,080-gene axis and source hashes.

On inference controls only, normalize each cell to CPM, select genes with arithmetic mean CPM>5 and positive log1p-CPM variance, standardize log1p CPM, and fit two principal components with fixed numerical settings/seed. Fix each component's sign by its largest absolute loading. Each template splits cells at that component's **inference** median; apply these frozen rules to the reference pool. Reject unsupported bins instead of changing thresholds. PC1/PC2 shifts are deliberate composition positive controls, not broad unseen-response validation and not a general model of cell-intrinsic knockdowns.

For a reference population with N_high and N_low cells, a 400-cell group has high-state count K distributed as

`P(K=k) ∝ choose(N_high,k) × choose(N_low,400−k) × exp(τk)`.

Draw K from this exact finite support, then choose distinct real cells uniformly within each state. At τ=0 this is exactly a simple random subset of 400 unmodified reference cells. The alias table is:

| Template | NULL | WEAK | MODERATE | STRONG |
| --- | ---: | ---: | ---: | ---: |
| PC1 | 0 | +log(1.25) | +log(2) | +log(4) |
| PC2 | 0 | −log(1.25) | −log(2) | −log(4) |

Let w=E[K]/400. Exact arithmetic probability is `w E[X/L|high]+(1−w)E[X/L|low]`. Exact expected raw counts and depth use the same weights; bulk probability is their ratio. It is **not** the expectation of a finite-group pooled ratio. Truth contrasts use the reference pool's own unperturbed population moments; τ=0 has exactly zero contrast. Differences from the separate realized comparison and inference pools are reported as finite-baseline nuisance, never labeled biological effects. Reference-cell reuse across aliases/seeds is reported; within a group all 400 cells are distinct.

## Independent latent Dirichlet–multinomial family

Define four fixed latent states using inference PC1 quartiles. Each of the 4,000 original inference cells supplies a joint `(state, integer library depth)` atom, sampled uniformly; thus depth/state dependence is retained. State mean compositions are arithmetic cell-probability means, shrunk 5% toward the global mean. These parameters use controls only.

For a sampled atom `(s,L)`, independently draw `θ ~ Dirichlet(κ μ_s)` on its positive support and `X ~ Multinomial(L, θ)`, with **κ=20,000** fixed before outcomes. This yields about 3.5 times the multinomial variance at depth 50,000; κ was chosen prospectively for an identifiable baseline diagnostic, not tuned to recovery. It does not assert realistic dispersion for every gene. Generate synthetic inference controls, comparison controls and reference groups from separate named RNG streams. Their same-law NULLs share no generated cells or random draws. Actual H1 inference counts are parameter-estimation data, not the synthetic family's comparison population.

Two programs each have 32 up and 32 down genes selected deterministically from inference-only metadata. BROAD uses mean CPM10–1000; ZERO_RICH uses CPM10–100 and detection≤0.5, excluding BROAD genes. Within each eligible set, pair adjacent mean-expression genes and choose 32 pairs by a frozen hash ordering/sign rule. Fail if support is insufficient. ZERO_RICH probes activation at **individual donor zeros**, not population-zero genes outside the official DE gate.

At strengths `f=0,0.1,0.3,0.6`, transfer `f × min(up-block mass,down-block mass)` from the down block to the up block in each state, proportional to each block's original probabilities. Preserve closure, compute truth from the final normalized conditional probabilities, and retain the directly modified mask separately. Marginal-null classification requires unchanged conditional state probabilities (numerical tolerance1e−14), concentration and the complete state/depth law; equality of an aggregate mean alone is insufficient.

For joint atoms with weights w_i, exact `p=Σw_i μ_state(i)`, expected counts `c=Σw_i L_i μ_state(i)`, expected depth `d=Σw_i L_i`, and `b=c/d`. Inferential/reference sample moments are distinct from these known population quantities.

## Candidate requests and interpretation

Truth artifacts retain absolute oracle p/b and own baseline p0/b0. Root's candidate builder uses the same **anchored probability-delta** requests for every eligible emitter: `normalize(max(empirical_inference_p0 + oracle_p − oracle_p0,0))`, and similarly for bulk. An explicit NULL branch returns the empirical inference baselines unchanged. Record clipping mass, closure changes, requested-versus-oracle deltas, population-fit drift, conditional sample expectation and realized count profiles. This is not exact LFC transfer; clipping can alter weak signs. Inference-baseline estimation error is separated from count realization. The truth generator never uses this request mapping.

Known LFC is computed after closure from arithmetic CPM with epsilon1e−9; bulk effects use log1p(50,000×bulk probability). Sign adjudication requires absolute true LFC≥0.1. Report continuous error and separate bins `(0,0.1)`, `[0.1,0.25)`, `[0.25,0.5)`, `≥0.5`; tiny composition effects are not called stable. There is no biological on-target gene in these aliases.

Use comparison mean CPM>5, unchanged per-alias BH and q<0.05. Report reference call power against the independent oracle, sign recovery, candidate/reference call overlap, true signals outside the realized comparison gate and false calls on genuine NULLs. Additional DM marginal-null false calls use the conditional-law mask above. A population mean effect does not imply deterministic significance in one sample.

Every independent reference NULL group must have **at most five** significant calls. Otherwise mark that family/seed's recovery interpretation inconclusive, retain all outputs, and do not retune. This is pragmatic diagnostic triage, not an FDR guarantee or model-promotion rule. Held-out templates/seeds, two reference families and independent count generation reduce circularity; they do not establish biological transfer or authorize reserve access/submission.
