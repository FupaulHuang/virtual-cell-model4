"""Within-depth integer contrast amplification with explicit moment constraints.

This is a dispersion hypothesis, not a claim of realistic perturbed-cell biology.
The requested multiplier applies to the aggregate CPM residual STANDARD DEVIATION
(sqrt of the summed per-gene sample variances), not to variance. Individual genes
need not achieve that multiplier. Low-count genes and tight CPM constraints can
make a request unattainable; the returned status reports that explicitly.

Each move reallocates counts between two similar-depth cells. Its integer delta
sums to zero, preserving both row depths and every gene's raw count total. The
remaining *accumulated* arithmetic-CPM budget bounds each gene's delta before
rounding. Final diagnostics independently recompute all constraints from counts.
No random expression noise is injected; randomness only selects nearby pairs and
breaks balanced-rounding ties. The identity arm bypasses all changes.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import time

import numpy as np


class ConstraintViolation(ValueError):
    """An emitted result failed a required count/moment invariant."""


@dataclass
class DispersionResult:
    counts: np.ndarray
    diagnostics: dict
    gene_stats: dict[str, np.ndarray]


def _input_counts(counts, max_counts_per_cell):
    a = np.asarray(counts)
    if a.ndim != 2 or min(a.shape) < 1 or a.dtype.kind not in 'iuf':
        raise ValueError('counts must be a nonempty numeric cell-by-gene matrix')
    if (not np.isfinite(a).all() or (a < 0).any()
            or (a > max_counts_per_cell).any() or not np.equal(a, np.floor(a)).all()):
        raise ValueError('counts must be finite nonnegative integers within the cell cap')
    x = a.astype(np.int64, copy=True)
    depths = x.sum(axis=1, dtype=np.int64)
    if (depths <= 0).any() or (depths > max_counts_per_cell).any():
        raise ValueError('each cell needs positive depth within max_counts_per_cell')
    return a, x, depths


def _moments(x, depths):
    cpm = x.astype(np.float64) * (1e6 / depths)[:, None]
    mean = cpm.mean(axis=0, dtype=np.float64)
    variance = cpm.var(axis=0, ddof=1, dtype=np.float64) if len(x) > 1 else np.zeros(x.shape[1])
    return mean, variance


def _duplicates(x):
    seen = {}
    for row in x:
        key = hashlib.blake2b(np.ascontiguousarray(row).view(np.uint8), digest_size=16).digest()
        seen[key] = seen.get(key, 0) + 1
    return {'unique_cells': len(seen), 'duplicate_cells': len(x) - len(seen),
            'largest_identical_cell_group': max(seen.values(), default=0)}


def _integer_box_projection(desired, lower, upper, rng):
    """Nearest box-constrained real vector with sum zero, then balanced rounding.

    Bounds are integers and contain zero. Flooring the real projection and
    rounding up exactly -sum(floor) entries preserves the zero sum and bounds.
    """
    if np.any(lower > 0) or np.any(upper < 0) or np.any(lower > upper):
        raise ConstraintViolation('integer box does not contain the feasible zero move')
    movable = upper > lower
    result = np.zeros(len(desired), dtype=np.int64)
    if not movable.any():
        return result
    v = desired[movable]
    lo = lower[movable]
    hi = upper[movable]
    left = float(np.min(v - hi))
    right = float(np.max(v - lo))
    for _ in range(55):
        midpoint = (left + right) / 2
        projected = np.clip(v - midpoint, lo, hi)
        if projected.sum(dtype=np.float64) > 0:
            left = midpoint
        else:
            right = midpoint
    projected = np.clip(v - (left + right) / 2, lo, hi)
    integer = np.floor(projected).astype(np.int64)
    need = -int(integer.sum(dtype=np.int64))
    eligible = np.flatnonzero(integer < hi)
    if need < 0 or need > len(eligible):
        raise ConstraintViolation('balanced rounding could not preserve the row total')
    if need:
        # Random permutation makes equal fractional ties reproducibly unbiased;
        # the primary ordering remains the actual fractional remainder.
        eligible = rng.permutation(eligible)
        score = projected[eligible] - integer[eligible]
        chosen = eligible[np.argsort(-score, kind='stable')[:need]]
        integer[chosen] += 1
    if integer.sum() != 0 or (integer < lo).any() or (integer > hi).any():
        raise ConstraintViolation('integer projection violated its margins or bounds')
    result[movable] = integer
    return result


def _pairs(order, pass_index, pair_window, rng):
    if pass_index == 0:
        return [(int(order[i]), int(order[i + 1])) for i in range(0, len(order) - 1, 2)]
    # Offset the windows on alternate passes so their boundaries are not fixed.
    shift = (pass_index % pair_window)
    groups = [order[:shift]] if shift else []
    groups += [order[i:i + pair_window] for i in range(shift, len(order), pair_window)]
    result = []
    for group in groups:
        shuffled = rng.permutation(group)
        result.extend((int(shuffled[i]), int(shuffled[i + 1]))
                      for i in range(0, len(shuffled) - 1, 2))
    return result


def redistribute_dispersion(counts, multiplier, *, seed=20261021,
        cpm_rtol=1e-4, cpm_atol=1e-8, max_passes=8, pair_window=4,
        max_relative_depth_gap=0.01, spread_rtol=0.01,
        diagnostic_gene_mask=None, max_counts_per_cell=1_000_000, max_nnz=None):
    """Return ``DispersionResult(counts, diagnostics, gene_stats)``.

    Required per-gene CPM tolerance is max(cpm_atol, cpm_rtol*original_mean).
    Exact row/column sums always hold, including for a partial/unattainable result.
    "Unattainable" means this bounded pair search found no change, not a proof
    that every possible count-matrix construction is infeasible.
    No tolerance is relaxed. ``status`` is identity, target_reached, partial, or
    unattainable; ``constraints_passed`` does not imply ``target_reached``.
    Density and duplicate-cell changes are reported separately as quality flags.
    ``max_nnz``, if supplied, is a hard per-block limit, not the whole-file limit.
    """
    start = time.monotonic()
    if not np.isfinite(multiplier) or multiplier < 1:
        raise ValueError('multiplier must be finite and >=1 (a standard-deviation multiplier)')
    if any(not np.isfinite(v) or v < 0 for v in [cpm_rtol, cpm_atol, max_relative_depth_gap, spread_rtol]):
        raise ValueError('tolerances and depth-gap limit must be finite and nonnegative')
    if not isinstance(max_passes, int) or max_passes < 0 or not isinstance(pair_window, int) or pair_window < 2:
        raise ValueError('max_passes must be >=0 and pair_window must be >=2')
    if not isinstance(max_counts_per_cell, int) or max_counts_per_cell < 1:
        raise ValueError('max_counts_per_cell must be a positive integer')
    if max_nnz is not None and (not isinstance(max_nnz, int) or max_nnz < 0):
        raise ValueError('max_nnz must be a nonnegative integer or None')
    original, x, depths = _input_counts(counts, max_counts_per_cell)
    original_int = x.copy()
    mean0, var0 = _moments(x, depths)
    tolerance = np.maximum(cpm_atol, cpm_rtol * mean0)
    if not np.isfinite(tolerance).all():
        raise ValueError('CPM tolerance overflowed')
    budget = np.maximum(tolerance - (1e-10 + 32*np.finfo(float).eps*np.maximum(mean0, 1)), 0)
    mask = np.ones(x.shape[1], dtype=bool) if diagnostic_gene_mask is None else np.asarray(diagnostic_gene_mask)
    if mask.shape != (x.shape[1],) or mask.dtype != np.bool_ or not mask.any():
        raise ValueError('diagnostic_gene_mask must be a nonempty boolean gene mask')
    base_trace = float(var0[mask].sum())
    current_trace = base_trace
    target_trace = base_trace * float(multiplier)**2
    drift = np.zeros(x.shape[1])
    rng = np.random.default_rng(seed)
    order = np.argsort(depths, kind='stable')
    changes = attempted = skipped_gap = accepted = completed_passes = 0
    pair_gain = np.sqrt(2*float(multiplier)**2 - 1) - 1
    if multiplier > 1 and base_trace > 0 and len(x) > 1:
        for pass_index in range(max_passes):
            for i, j in _pairs(order, pass_index, pair_window, rng):
                if current_trace >= target_trace * (1-spread_rtol)**2:
                    break
                attempted += 1
                if abs(float(depths[i])-float(depths[j])) / min(depths[i], depths[j]) > max_relative_depth_gap:
                    skipped_gap += 1
                    continue
                ai, aj = 1e6/float(depths[i]), 1e6/float(depths[j])
                coefficient = (ai-aj)/len(x)
                lo, hi = -x[i].copy(), x[j].copy()
                if coefficient != 0:
                    a = (-budget-drift)/coefficient
                    b = (budget-drift)/coefficient
                    lo = np.maximum(lo, np.ceil(np.clip(np.minimum(a,b), -max_counts_per_cell, max_counts_per_cell)).astype(np.int64))
                    hi = np.minimum(hi, np.floor(np.clip(np.maximum(a,b), -max_counts_per_cell, max_counts_per_cell)).astype(np.int64))
                    # Current counts remain feasible even if floating accumulation
                    # sits just beyond the deliberately conservative inner budget.
                    lo = np.minimum(lo,0); hi = np.maximum(hi,0)
                residual = x[i] - (float(depths[i])/(float(depths[i])+float(depths[j]))) * (x[i]+x[j])
                gain = pair_gain
                delta = None
                increment = 0.0
                for _ in range(10):
                    candidate = _integer_box_projection(gain*residual, lo, hi, rng)
                    if not candidate.any():
                        break
                    dm = candidate * coefficient
                    dsq = (2*candidate*(x[i]*ai*ai-x[j]*aj*aj)
                           + candidate.astype(float)**2*(ai*ai+aj*aj))
                    dv = (dsq-len(x)*(2*(mean0+drift)*dm+dm*dm))/(len(x)-1)
                    increment = float(dv[mask].sum())
                    if increment <= 0:
                        break
                    if current_trace+increment <= target_trace*(1+spread_rtol)**2:
                        delta = candidate
                        break
                    gain *= 0.5
                if delta is None:
                    continue
                x[i] += delta; x[j] -= delta
                drift += delta*coefficient
                current_trace += increment
                accepted += 1; changes += int(np.abs(delta).sum())
            completed_passes = pass_index+1
            mean_now, var_now = _moments(x,depths)
            drift = mean_now-mean0
            current_trace = float(var_now[mask].sum())
            if (np.abs(drift)>tolerance).any():
                raise ConstraintViolation('recomputed arithmetic CPM exceeded its declared tolerance')
            if current_trace >= target_trace*(1-spread_rtol)**2:
                break
    mean1, var1 = _moments(x,depths)
    delta_mean = mean1-mean0
    row_exact = np.array_equal(x.sum(1,dtype=np.int64),original_int.sum(1,dtype=np.int64))
    column_exact = np.array_equal(x.sum(0,dtype=np.int64),original_int.sum(0,dtype=np.int64))
    cpm_pass = bool((np.abs(delta_mean)<=tolerance).all())
    if not row_exact or not column_exact or not cpm_pass or (x<0).any():
        raise ConstraintViolation('final independent count/moment constraint audit failed')
    nnz0, nnz1 = int(np.count_nonzero(original_int)), int(np.count_nonzero(x))
    if max_nnz is not None and nnz1 > max_nnz:
        raise ConstraintViolation('output exceeds the declared per-block nonzero limit')
    ratio = np.divide(var1,var0,out=np.full_like(var0,np.nan),where=var0>0)
    sd_ratio = np.sqrt(ratio)
    achieved = float(np.sqrt(var1[mask].sum()/base_trace)) if base_trace>0 else None
    reached = multiplier==1 or (achieved is not None and abs(achieved/multiplier-1)<=spread_rtol)
    changed = not np.array_equal(x,original_int)
    status = 'identity' if multiplier==1 else 'target_reached' if reached else 'partial' if changed else 'unattainable'
    dup0, dup1 = _duplicates(original_int), _duplicates(x)
    eligible = mask & (var0>0)
    quantiles = {str(q):float(np.quantile(sd_ratio[eligible],q)) for q in [0,0.1,0.5,0.9,1]} if eligible.any() else {}
    diagnostics = dict(status=status,constraints_passed=True,target_reached=bool(reached),
        requested_sd_multiplier=float(multiplier),requested_variance_multiplier=float(multiplier)**2,
        achieved_sd_multiplier=achieved,spread_rtol=float(spread_rtol),
        aggregate_sd_definition='sqrt(sum of per-gene sample CPM variances over diagnostic_gene_mask)',
        n_cells=x.shape[0],n_genes=x.shape[1],diagnostic_genes=int(mask.sum()),seed=int(seed),
        cpm_rtol=float(cpm_rtol),cpm_atol=float(cpm_atol),row_depths_exact=row_exact,gene_totals_exact=column_exact,
        arithmetic_cpm_within_tolerance=cpm_pass,max_abs_cpm_change=float(np.abs(delta_mean).max()),
        max_relative_cpm_change=float(np.max(np.divide(np.abs(delta_mean),mean0,out=np.zeros_like(mean0),where=mean0>0))),
        cpm_constraint_violations=int((np.abs(delta_mean)>tolerance).sum()),
        changed_entries=int(np.count_nonzero(x!=original_int)),moved_count_units=changes,
        attempted_pairs=attempted,accepted_pairs=accepted,skipped_depth_gap_pairs=skipped_gap,passes=completed_passes,
        max_relative_depth_gap=float(max_relative_depth_gap),pair_window=pair_window,
        nnz_before=nnz0,nnz_after=nnz1,density_before=nnz0/x.size,density_after=nnz1/x.size,
        duplicate_cells_before=dup0,duplicate_cells_after=dup1,
        duplicate_cells_increased=dup1['duplicate_cells']>dup0['duplicate_cells'],
        per_gene_sd_multiplier_quantiles=quantiles,genes_variance_increased=int((var1>var0+1e-12).sum()),
        genes_variance_decreased=int((var1<var0-1e-12).sum()),
        zero_variance_genes_before=int((var0==0).sum()),zero_total_genes=int((original_int.sum(0)==0).sum()),
        genes_zero_to_positive_variance=int(((var0==0)&(var1>0)).sum()),
        baseline_variance_trace=base_trace,output_variance_trace=float(var1[mask].sum()),
        feasibility_scope='configured near-depth integer pair search; not a global infeasibility proof',
        no_expression_jitter=True,elapsed_seconds=time.monotonic()-start,
        limitation='Exact count margins and bounded CPM do not guarantee unchanged MSE, reach, or measured NMAE; score independently.')
    stats=dict(mean_cpm_before=mean0,mean_cpm_after=mean1,mean_cpm_change=delta_mean,
        mean_cpm_tolerance=tolerance,variance_cpm_before=var0,variance_cpm_after=var1,
        variance_ratio=ratio,sd_multiplier=sd_ratio,diagnostic_gene_mask=mask.copy())
    # Identity preserves dtype and bytes. Other arms return int64 so an input
    # int8/float16 dtype cannot overflow after concentrating counts in one cell.
    return DispersionResult(original.copy() if multiplier==1 else x,diagnostics,stats)
