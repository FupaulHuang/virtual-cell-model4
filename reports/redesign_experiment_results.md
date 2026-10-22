# Response × emitter development experiment

Eight frozen arms; 50 H1 Training targets, 400 real/predicted cells per target and 4,000 controls. H1 responses were excluded from adapter labels and inputs. These are raw metrics on a development cohort.

| Arm | Seed | PDS ↑ | MSE ↓ | NMAE ↓ | Fidelity ↑ | Reach ↑ | Jaccard ↑ | Own-excluded mean error ↓ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| current_current | 20261021 | 0.910204082 | 0.970348453 | 0.962901741 | 0.510809091 | 0.115049906 | 0.092176112 | 0.002221433 |
| current_current | 20261022 | 0.911020408 | 0.970366480 | 0.962250747 | 0.511310124 | 0.112123311 | 0.091994105 | 0.002221344 |
| current_preserving | 20261021 | 0.749795918 | 1.067804776 | 0.967813378 | 0.390880174 | 0.129449879 | 0.045392049 | 0.002580059 |
| current_preserving | 20261022 | 0.803265306 | 1.065414166 | 0.972636760 | 0.383271216 | 0.111333543 | 0.044627611 | 0.002575798 |
| adapter_current | 20261021 | 0.793061224 | 1.006788336 | 0.992287679 | 0.487431835 | 0.080238996 | 0.070286946 | 0.002140571 |
| adapter_current | 20261022 | 0.794285714 | 1.006797133 | 0.992411093 | 0.478557737 | 0.066771050 | 0.067812878 | 0.002140590 |
| adapter_preserving | 20261021 | 0.597142857 | 1.012855144 | 0.991182317 | 0.007700884 | 0.064554271 | 0.000025652 | 0.002450782 |
| adapter_preserving | 20261022 | 0.629387755 | 1.014523076 | 0.996497934 | 0.002535830 | 0.061965588 | 0.000006108 | 0.002447962 |

## current_preserving: FAIL STRICT

- Seed 20261021: FAIL; all-six nonregression=False; fidelity/Jaccard improvement >1e-8=False.
  - PDS: gain -0.160408163, paired 95% [-0.236734694, -0.086530612].
  - Expression MSE: gain -0.097456323; Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap.
  - LFC NMAE: gain -0.004911637, paired 95% [-0.011482597, +0.002095094].
  - Direction fidelity: gain -0.119928916, paired 95% [-0.167186993, -0.071993066].
  - Direction reach: gain +0.014399973, paired 95% [-0.022594219, +0.052120623].
  - Significance Jaccard: gain -0.046784063, paired 95% [-0.069624270, -0.025027798].
  - Own-excluded mean-error gain -0.000358627; paired 95% [-0.00037917652635954883, -0.0003375756171835623].
- Seed 20261022: FAIL; all-six nonregression=False; fidelity/Jaccard improvement >1e-8=False.
  - PDS: gain -0.107755102, paired 95% [-0.177959184, -0.038336735].
  - Expression MSE: gain -0.095047687; Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap.
  - LFC NMAE: gain -0.010386013, paired 95% [-0.018808177, -0.002743579].
  - Direction fidelity: gain -0.128038908, paired 95% [-0.178144263, -0.080811407].
  - Direction reach: gain -0.000789769, paired 95% [-0.040468207, +0.037147503].
  - Significance Jaccard: gain -0.047366494, paired 95% [-0.071438585, -0.024877499].
  - Own-excluded mean-error gain -0.000354454; paired 95% [-0.0003698522498600297, -0.0003393884163221871].

## adapter_current: FAIL STRICT

- Seed 20261021: FAIL; all-six nonregression=False; fidelity/Jaccard improvement >1e-8=False.
  - PDS: gain -0.117142857, paired 95% [-0.164081633, -0.074275510].
  - Expression MSE: gain -0.036439883; Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap.
  - LFC NMAE: gain -0.029385938, paired 95% [-0.043207730, -0.015489493].
  - Direction fidelity: gain -0.023377256, paired 95% [-0.042657786, -0.007062128].
  - Direction reach: gain -0.034810909, paired 95% [-0.107741452, +0.026973118].
  - Significance Jaccard: gain -0.021889166, paired 95% [-0.035596333, -0.009904252].
  - Own-excluded mean-error gain +0.000080862; paired 95% [4.9539260861540154e-05, 0.00010984615791175343].
- Seed 20261022: FAIL; all-six nonregression=False; fidelity/Jaccard improvement >1e-8=False.
  - PDS: gain -0.116734694, paired 95% [-0.163275510, -0.074683673].
  - Expression MSE: gain -0.036430653; Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap.
  - LFC NMAE: gain -0.030160346, paired 95% [-0.044222717, -0.016279564].
  - Direction fidelity: gain -0.032752387, paired 95% [-0.056902700, -0.011299255].
  - Direction reach: gain -0.045352261, paired 95% [-0.114776340, +0.012108026].
  - Significance Jaccard: gain -0.024181227, paired 95% [-0.040026885, -0.010090327].
  - Own-excluded mean-error gain +0.000080754; paired 95% [4.946562899391453e-05, 0.00010982138474196634].

## adapter_preserving: FAIL STRICT

- Seed 20261021: FAIL; all-six nonregression=False; fidelity/Jaccard improvement >1e-8=False.
  - PDS: gain -0.313061224, paired 95% [-0.401224490, -0.231000000].
  - Expression MSE: gain -0.042506691; Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap.
  - LFC NMAE: gain -0.028280576, paired 95% [-0.043142470, -0.011740343].
  - Direction fidelity: gain -0.503108207, paired 95% [-0.516693608, -0.484615810].
  - Direction reach: gain -0.050495634, paired 95% [-0.110832343, -0.002513436].
  - Significance Jaccard: gain -0.092150460, paired 95% [-0.124269680, -0.061455834].
  - Own-excluded mean-error gain -0.000229349; paired 95% [-0.0002605524332528685, -0.00020061703831298014].
- Seed 20261022: FAIL; all-six nonregression=False; fidelity/Jaccard improvement >1e-8=False.
  - PDS: gain -0.281632653, paired 95% [-0.365724490, -0.197510204].
  - Expression MSE: gain -0.044156597; Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap.
  - LFC NMAE: gain -0.034247186, paired 95% [-0.054009076, -0.016553219].
  - Direction fidelity: gain -0.508774294, paired 95% [-0.518023998, -0.498784210].
  - Direction reach: gain -0.050157723, paired 95% [-0.117574765, +0.004695890].
  - Significance Jaccard: gain -0.091987997, paired 95% [-0.123879151, -0.061397388].
  - Own-excluded mean-error gain -0.000226618; paired 95% [-0.0002632273653944978, -0.00019235214830329176].

No recipe passes both seeds: retain incumbent and stop before second-context or reserved-response extraction.

The point gate requires all six raw metrics nonworse within 1e-8, improvement greater than 1e-8 in fidelity or Jaccard, identical finite cohorts, and a pass at both seeds. Intervals are descriptive and cannot turn a failed point gate into a pass. No global normalized-MSE target bootstrap or cross-panel anchors are used.

Adapter response includes learned transfer, exact pseudocount inversion and omission of CD4. The emitter changes donor construction and count transport. Requested profiles are held fixed across emitters; realized arithmetic/bulk drift is measured. The JSON retains donor/axis/hash checks, coverage, mean-error before and after own-target exclusion, correction decomposition, DE calls and per-gene drift artifacts. See receipt metadata notes (local artifact omitted: `redesign_receipt_metadata_notes.json`) for inherited legacy source weights/common subtraction versus actual adapter semantics.

- Previously observed H1 Training development cohort; not untouched destination validation.
- Paired target intervals are descriptive with shared controls/correlated perturbations and no multiplicity correction.
- Two emission seeds do not estimate all generation uncertainty; neither seed is selected after observing outcomes.
- The adapter uses only three non-H1 source contexts and LOSO-selected regularization; its confidence is not DE significance.
- Adapter versus current response also changes CD4 inclusion and exact arithmetic inverse; those are part of the registered response intervention.
- Adapter prediction receipts inherit legacy raw-source weights 2/1/1; actual adapter pooling is equal over valid raw families before learned feature coefficients, as fixed by adapter code/model/profiles.
- Point nonregression is an operational gate, not statistical proof of no harm.
- Global capped MSE has a shared correction budget and receives no naive target bootstrap.
- Native gene panels and this 18,080-gene reference differ from the challenge; no leaderboard gain is inferred.
