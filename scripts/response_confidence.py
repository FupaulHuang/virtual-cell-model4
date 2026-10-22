"""Experimental source-agreement mask for a measured response ensemble.

This is a modeling hypothesis, not a validated measure of biological certainty.
Pass centered source effects in the same log2-fold-change units and provide the
actual fused effect after the caller's source weighting. Sources are one vote
each; correlated replicates must be aggregated into independent source families
before invoking this function.
"""
from __future__ import annotations

import numpy as np


def confidence_gate(effects, measured, fused_effect, min_abs_effect=0.02, min_agreement=0.75):
    """Return a boolean target-by-gene gate and serializable diagnostics.

    For a pair measured in zero or one source, the gate stays open: agreement
    evidence is absent. Otherwise, at least two measured sources must have an
    absolute effect >= ``min_abs_effect`` and a nonzero sign equal to the actual
    fused effect's sign. Their count divided by *all measured sources*, including
    weak and opposing sources, must be >= ``min_agreement``. A zero fused effect
    has no agreeing sign and closes the gate if at least two sources measured it.

    Missing pairs never contribute a vote or enter the denominator. Their effect
    slots may contain arbitrary values or NaNs. Nonfinite measured effects and
    nonfinite fused effects are rejected, since silently removing them would
    change the evidence denominator. The caller applies this same gate to both
    arithmetic-mean and bulk response heads; this function does not mutate data.
    """
    values = np.asarray(effects)
    mask = np.asarray(measured)
    fused = np.asarray(fused_effect)
    if values.ndim != 3 or any(size == 0 for size in values.shape[1:]):
        raise ValueError('effects must have nonempty target/gene axes and shape [source,target,gene]')
    if mask.shape != values.shape or mask.dtype != np.bool_:
        raise ValueError('measured must be a boolean array matching effects')
    if fused.shape != values.shape[1:]:
        raise ValueError('fused_effect must have shape [target,gene]')
    if not np.isfinite(min_abs_effect) or min_abs_effect < 0:
        raise ValueError('min_abs_effect must be finite and nonnegative')
    if not np.isfinite(min_agreement) or not 0 <= min_agreement <= 1:
        raise ValueError('min_agreement must be finite and between zero and one')
    if not np.isfinite(values[mask]).all() or not np.isfinite(fused).all():
        raise ValueError('Measured and fused effects must be finite')
    # Replace absent slots before sign/abs operations, so NaNs and infinities in
    # genuinely missing source measurements cannot influence voting.
    present = np.where(mask, values, 0.0)
    measured_count = mask.sum(axis=0, dtype=np.int32)
    strong = mask & (np.abs(present) >= min_abs_effect) & (present != 0)
    strong_count = strong.sum(axis=0, dtype=np.int32)
    fused_sign = np.sign(fused)
    agreeing = strong & (np.sign(present) == fused_sign[None, :, :]) & (fused_sign[None, :, :] != 0)
    agreement_count = agreeing.sum(axis=0, dtype=np.int32)
    fraction = np.divide(agreement_count, measured_count,
                         out=np.zeros(fused.shape, dtype=np.float64), where=measured_count > 0)
    eligible = measured_count >= 2
    passed = (agreement_count >= 2) & (fraction >= min_agreement)
    gate = (~eligible) | passed
    diagnostics = {
        'rule': '>=2 nonweak nonzero measured-source votes match supplied fused sign, and matching/measured >= threshold; <=1 measured stays open',
        'min_abs_effect': float(min_abs_effect), 'min_agreement': float(min_agreement),
        'source_count': int(values.shape[0]), 'target_count': int(values.shape[1]),
        'gene_count': int(values.shape[2]), 'total_pairs': int(fused.size),
        'unmeasured_pairs': int((measured_count == 0).sum()),
        'single_source_pairs': int((measured_count == 1).sum()),
        'eligible_pairs': int(eligible.sum()),
        'retained_eligible_pairs': int((eligible & gate).sum()),
        'closed_pairs': int((~gate).sum()),
        'eligible_with_fewer_than_two_nonweak_votes': int((eligible & (strong_count < 2)).sum()),
        'eligible_with_fewer_than_two_agreeing_votes': int((eligible & (agreement_count < 2)).sum()),
        'eligible_below_agreement_fraction': int((eligible & (fraction < min_agreement)).sum()),
        'eligible_with_zero_fused_effect': int((eligible & (fused == 0)).sum()),
        'retained_fraction_among_eligible': float((eligible & gate).sum() / eligible.sum()) if eligible.any() else None,
        'measured_source_count_histogram': {str(int(k)): int(v) for k, v in zip(*np.unique(measured_count, return_counts=True))},
        'eligible_pairs_by_target': eligible.sum(axis=1).astype(int).tolist(),
        'closed_pairs_by_target': (~gate).sum(axis=1).astype(int).tolist(),
        'validation_status': 'Unvalidated until compared using all six held-out-context metrics.',
    }
    return gate, diagnostics
