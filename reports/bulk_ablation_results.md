# Bulk-response development results

Fixed 50-target H1 Training panel, 400 real/predicted cells per target and 4,000 real controls. Values are raw official metrics. No overall raw score or cross-panel anchor scaling.

| Metric | Incumbent | common05 | bulk015 | bulk060 | Better |
| --- | ---: | ---: | ---: | ---: | --- |
| PDS | 0.910204082 | 0.908979592 | 0.904081633 | 0.908979592 | higher |
| Expression MSE | 0.970348453 | 0.969871583 | 0.993564953 | 1.170655702 | lower |
| LFC NMAE | 0.962901741 | 0.962755307 | 0.963253510 | 0.962400446 | lower |
| Direction fidelity | 0.510809091 | 0.511297521 | 0.510522401 | 0.510937685 | higher |
| Direction reach | 0.115049906 | 0.114257266 | 0.118680434 | 0.109785434 | higher |
| Significance Jaccard | 0.092176112 | 0.092107842 | 0.092129562 | 0.091830550 | higher |

| Profile error (lower better) | Incumbent | common05 | bulk015 | bulk060 |
| --- | ---: | ---: | ---: | ---: |
| requested | 0.002228714 | 0.002227911 | 0.002150873 | 0.002622441 |
| projected | 0.002223869 | 0.002223079 | 0.002149918 | 0.002595662 |
| fitted | 0.002221281 | 0.002220436 | 0.002149540 | 0.002582649 |
| realized | 0.002221433 | 0.002220581 | 0.002149668 | 0.002582763 |

## common05: FAIL STRICT

Realized mean-error gain: +0.000000851; paired 95% interval [+0.000000283, +0.000001439]. All five protected metrics nonworse: False. Identical finite cohorts: True. Eligible for reserved confirmation: False.

| Metric | Gain in better direction | Paired 95% interval |
| --- | ---: | --- |
| PDS | -0.001224490 | [-0.006938776, +0.004897959] |
| Expression MSE | +0.000476870 | Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap. |
| LFC NMAE | +0.000146434 | [-0.000742625, +0.000959763] |
| Direction fidelity | +0.000488430 | [-0.000230311, +0.001236793] |
| Direction reach | -0.000792639 | [-0.009145021, +0.007699603] |
| Significance Jaccard | -0.000068270 | [-0.000337224, +0.000178346] |

Maximum realized CPM drift: 1.01619. Mean uncorrected error change: -8.51319078e-07; capped numerator change: -8.50036279e-07; effective prediction-correction change: -1.28279835e-09. Both mean response and the global correction budget may move.

Mean realized response/residual dot per gene: -0.000112017636 to -0.000111708535. Positive alignment is a local directional diagnostic, not a fitted amplitude or an extra selection gate.

## bulk015: FAIL STRICT

Realized mean-error gain: +0.000071765; paired 95% interval [+0.000050154, +0.000092069]. All five protected metrics nonworse: False. Identical finite cohorts: True. Eligible for reserved confirmation: False.

| Metric | Gain in better direction | Paired 95% interval |
| --- | ---: | --- |
| PDS | -0.006122449 | [-0.011836735, -0.001632653] |
| Expression MSE | -0.023216500 | Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap. |
| LFC NMAE | -0.000351769 | [-0.001231832, +0.000413918] |
| Direction fidelity | -0.000286690 | [-0.001323436, +0.000760692] |
| Direction reach | +0.003630529 | [-0.005673292, +0.013296071] |
| Significance Jaccard | -0.000046550 | [-0.000310281, +0.000196706] |

Maximum realized CPM drift: 0.924975. Mean uncorrected error change: -7.1764665e-05; capped numerator change: +4.13841893e-05; effective prediction-correction change: -0.000113148854. Both mean response and the global correction budget may move.

Mean realized response/residual dot per gene: -0.000112017636 to -1.84137033e-05. Positive alignment is a local directional diagnostic, not a fitted amplitude or an extra selection gate.

## bulk060: FAIL STRICT

Realized mean-error gain: -0.000361331; paired 95% interval [-0.000419770, -0.000308458]. All five protected metrics nonworse: False. Identical finite cohorts: True. Eligible for reserved confirmation: False.

| Metric | Gain in better direction | Paired 95% interval |
| --- | ---: | --- |
| PDS | -0.001224490 | [-0.004897959, +0.002448980] |
| Expression MSE | -0.200307249 | Global ratio with whole-panel correction budget: no naive per-target MSE bootstrap. |
| LFC NMAE | +0.000501295 | [-0.000411872, +0.001541147] |
| Direction fidelity | +0.000128594 | [-0.000926410, +0.001156078] |
| Direction reach | -0.005264471 | [-0.012210574, -0.000610728] |
| Significance Jaccard | -0.000345562 | [-0.000759470, +0.000023944] |

Maximum realized CPM drift: 1.01548. Mean uncorrected error change: +0.000361330583; capped numerator change: +0.000357054389; effective prediction-correction change: +4.27619483e-06. Both mean response and the global correction budget may move.

Mean realized response/residual dot per gene: -0.000112017636 to -0.000517185063. Positive alignment is a local directional diagnostic, not a fitted amplitude or an extra selection gate.

Deterministic reserve selection: none. No candidate qualifies: retain incumbent, leave the reserve unextracted and stop.

The JSON retains actual finite target IDs, source-coverage strata, response/residual alignment, strict DE calls/TP/FP/FN, terminal statuses and independent final-count/provenance checks.

- Previously observed single-context Training development cohort; this is not untouched biological validation.
- Paired target intervals are conditional/descriptive with shared controls and correlated perturbations; no multiplicity adjustment.
- Point nonregression within 1e-8 is not statistical proof of no harm.
- Fixed requested arithmetic profiles do not ensure identical realized means or DE tests.
- Four thousand controls and the native 18,080-gene axis differ from the official challenge; no leaderboard gain is inferred.
- The global MSE correction is not bootstrapped as independent per-target normalized errors; mean residual and correction effects are separated.
- No model is promoted or submitted by this analysis. Closed H1 Test/K562 confirmations remain closed.
