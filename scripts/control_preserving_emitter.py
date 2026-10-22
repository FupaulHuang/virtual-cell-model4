"""Population-calibrated, identity-preserving stochastic count transport.

Moment fitting adapts the numerical iteration in external/atlasshift/model.py
(MIT; see external/atlasshift/LICENSE). Unlike that emitter, calibration uses
unpooled population controls, the support prior vanishes at zero response,
and emission never forces empirical gene totals. Stochastic thinning coupled
to multinomial births preserves each donor's depth and has the fitted counts
as its conditional expectation.
"""
import numpy as np


def validate_counts(counts, *, max_depth=1_000_000, max_nnz=4_750_000_000):
    x = np.asarray(counts)
    if x.ndim != 2 or not all(x.shape):
        raise ValueError('Counts must be a nonempty cell-by-gene matrix')
    if not np.isfinite(x).all() or (x < 0).any() or (x != np.floor(x)).any():
        raise ValueError('Counts must be finite nonnegative integers')
    depth = x.sum(1, dtype=np.float64)
    if (depth < 1).any() or (depth > max_depth).any():
        raise ValueError('Cell depth outside allowed bounds')
    if np.count_nonzero(x) > max_nnz:
        raise ValueError('Sparse entry cap exceeded')
    return depth.astype(np.int64)


def population_profiles(counts):
    depth = validate_counts(counts)
    probability = np.asarray(counts, dtype=np.float64) / depth[:, None]
    return probability.mean(0), np.asarray(counts).sum(0, dtype=np.float64) / depth.sum()


def project_bulk(probability, bulk, z):
    """Project only to necessary nonnegative depth bounds; identity is feasible."""
    p, b = np.asarray(probability), np.asarray(bulk)
    lower, upper = z.min() * p, z.max() * p
    if np.ptp(z) == 0:
        return p.copy()
    if np.all(b >= lower) and np.all(b <= upper):
        return b.copy()
    lo, hi = 0., 1.
    for _ in range(1024):
        if np.clip(b * hi, lower, upper).sum() >= 1 - 1e-14:
            break
        hi *= 2
    else:
        raise ValueError('Cannot bracket bulk projection')
    for _ in range(80):
        mid = (lo + hi) / 2
        if np.clip(b * mid, lower, upper).sum() < 1:
            lo = mid
        else:
            hi = mid
    return np.clip(b * ((lo + hi) / 2), lower, upper)


def fit_population(counts, requested_mean, requested_bulk, *, iterations=150,
                   tolerance=0.0002, max_fit_l1=0.001, max_projection_l1=0.03):
    """Return per-control expected probabilities and an auditable fit receipt.

Both moments are defined on the entire inference control population. Sampling
400 donors afterwards does not impose exact 400-cell mean or bulk totals.
The returned bulk moment is a ratio of expected sums; the expected arithmetic
profile is exact because stochastic transport preserves each row's depth.
"""
    depth = validate_counts(counts)
    x = np.asarray(counts, dtype=np.float64) / depth[:, None]
    p, requested = [np.asarray(a, dtype=np.float64) for a in (requested_mean, requested_bulk)]
    for a in (p, requested):
        if a.shape != (x.shape[1],) or not np.isfinite(a).all() or (a < 0).any():
            raise ValueError('Invalid requested profile')
        if not np.isclose(a.sum(), 1, atol=1e-8, rtol=0):
            raise ValueError('Requested probabilities must sum to one')
    if (type(iterations) is not int or iterations < 1
            or not np.isfinite([tolerance, max_fit_l1, max_projection_l1]).all()
            or not (0 < tolerance <= max_fit_l1) or max_projection_l1 < 0):
        raise ValueError('Invalid fitting limits')
    z = depth / depth.mean()
    baseline_mean, baseline_bulk = x.mean(0), z @ x / len(x)
    effect_size = max(float(np.abs(p - baseline_mean).sum()),
                      float(np.abs(requested - baseline_bulk).sum()))
    # Floating arithmetic differences from equivalent sparse reductions only.
    identity = (np.allclose(p, baseline_mean, atol=1e-15, rtol=0)
                and np.allclose(requested, baseline_bulk, atol=1e-15, rtol=0))
    if identity:
        return x, dict(identity=True, iterations=0, prior_fraction=0.,
                       effect_l1=effect_size, projection_l1=0.,
                       fitted_mean_l1=0., fitted_bulk_l1=0.)
    bulk = project_bulk(p, requested, z)
    projection = float(np.abs(bulk - requested).sum())
    if projection > max_projection_l1:
        raise ValueError(f'Bulk projection exceeds fixed limit: {projection}')
    # Response-dependent support creates donor-zero activation without any
    # nonzero smoothing floor at the identity limit.
    alpha = min(0.01, 0.01 * effect_size)
    x *= 1 - alpha
    x += alpha * p[None, :]
    desired = len(x) * p
    ratio = np.divide(bulk, p, out=np.ones_like(p), where=p > 0)
    for step in range(iterations):
        col = x.sum(0)
        x *= np.divide(desired, col, out=np.zeros_like(col), where=col > 0)
        first = z @ x
        mean = np.divide(first, desired, out=np.ones_like(first), where=desired > 0)
        second = (z * z) @ x
        var = np.maximum(np.divide(second, desired, out=np.zeros_like(second), where=desired > 0) - mean**2, 0)
        tilt = np.divide(ratio - mean, var, out=np.zeros_like(mean), where=var > 1e-12)
        tilt = np.clip(tilt, -0.95 / np.maximum(z.max() - mean, 1e-12),
                       0.95 / np.maximum(mean - z.min(), 1e-12))
        x *= 1 + tilt[None, :] * (z[:, None] - mean[None, :])
        x /= np.maximum(x.sum(1, keepdims=True), 1e-30)
        if step % 5 == 4:
            error_mean = float(np.abs(x.mean(0) - p).sum())
            error_bulk = float(np.abs(z @ x / len(x) - bulk).sum())
            if max(error_mean, error_bulk) < tolerance:
                break
    error_mean = float(np.abs(x.mean(0) - p).sum())
    error_bulk = float(np.abs(z @ x / len(x) - bulk).sum())
    if max(error_mean, error_bulk) > max_fit_l1:
        raise ValueError(f'Population fit failed: mean={error_mean}, bulk={error_bulk}')
    if not np.isfinite(x).all() or (x < 0).any() or not np.allclose(x.sum(1), 1, atol=1e-12, rtol=0):
        raise ValueError('Invalid fitted population probabilities')
    return x, dict(identity=False, iterations=step + 1, prior_fraction=alpha,
                   effect_l1=effect_size, projection_l1=projection,
                   fitted_mean_l1=error_mean, fitted_bulk_l1=error_bulk)


def emit_counts(donor_counts, expected_probability, seed, *, max_depth=1_000_000,
                max_nnz=4_750_000_000):
    """Coupled thinning/births: E[Y|donor]=depth*p and row sums are exact.

For deficit genes, thin independently with survival lambda/x. Redistribute
the total number removed among surplus genes with probabilities proportional
to lambda-x. Since deficits and surpluses sum equally, each birth expectation
is exactly its surplus. No independent Poisson change to library size occurs.
"""
    depth = validate_counts(donor_counts, max_depth=max_depth, max_nnz=max_nnz)
    raw = np.asarray(donor_counts, dtype=np.int64)
    p = np.asarray(expected_probability, dtype=np.float64)
    if (p.shape != raw.shape or not np.isfinite(p).all() or (p < 0).any()
            or not np.allclose(p.sum(1), 1, atol=1e-10, rtol=0)):
        raise ValueError('Expected probabilities must have matching axes and unit row sums')
    result = raw.copy()
    rng = np.random.default_rng(seed)
    for i, d in enumerate(depth):
        if np.array_equal(p[i], raw[i] / d):
            continue
        expected = p[i] * d
        delta = expected - raw[i]
        # This detects exact mathematical identity after floating division.
        # It is not a biological effect threshold or a variance adjustment.
        if np.max(np.abs(delta)) <= 1e-10:
            continue
        down, up = delta < 0, delta > 0
        lost = rng.binomial(raw[i, down], np.clip(-delta[down] / raw[i, down], 0, 1))
        moved = int(lost.sum())
        result[i, down] -= lost
        if moved:
            if not up.any() or delta[up].sum() <= 0:
                raise ValueError('Removed counts have no birth support')
            result[i, up] += rng.multinomial(moved, delta[up] / delta[up].sum())
    if not np.array_equal(result.sum(1), depth):
        raise AssertionError('Stochastic transport changed a donor library size')
    validate_counts(result, max_depth=max_depth, max_nnz=max_nnz)
    return result.astype(np.int32)
