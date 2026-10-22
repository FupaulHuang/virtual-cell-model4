"""Fixed-size raw-donor sampling with known equal inclusion probabilities."""
import numpy as np


def depth_blocks(depths, cells=400):
    d = np.asarray(depths)
    if (d.ndim != 1 or not len(d) or not np.isfinite(d).all() or (d <= 0).any()
            or (d != np.floor(d)).any() or type(cells) is not int or cells < 1
            or len(d) % cells or len(d) < cells):
        raise ValueError('Positive integer depths and a population divisible by sample size required')
    # Ties retain the original row order of the frozen control file.
    return np.argsort(d, kind='stable').reshape(cells, len(d)//cells)


def select_donors(depths, seed, policy, cells=400):
    blocks = depth_blocks(depths, cells)
    rng = np.random.default_rng(seed)
    if policy == 'uniform':
        return rng.choice(blocks.size, cells, replace=False)
    if policy == 'depth_stratified':
        return blocks[np.arange(cells), rng.integers(0, blocks.shape[1], cells)]
    raise ValueError('Unknown donor policy')


def exact_mean_sampling_variance(values, blocks):
    """Per-feature variance of the unweighted sample mean under both designs.

    Population groups have equal size, one independent uniform draw per group.
    This applies to linear estimators (CPM means/count means), not bulk ratios.
    """
    x = np.asarray(values, dtype=np.float64)
    groups = np.asarray(blocks)
    if (x.ndim != 2 or not np.isfinite(x).all() or groups.ndim != 2
            or not np.issubdtype(groups.dtype, np.integer)
            or not np.array_equal(np.sort(groups.ravel()), np.arange(len(x)))
            or len(x) < 2):
        raise ValueError('Finite feature matrix and a complete equal-sized partition required')
    n, width = groups.shape
    uniform = (1 - n/len(x)) * x.var(0, ddof=1) / n
    # A single draw in each group has the finite population ddof=0 variance.
    stratified = np.zeros(x.shape[1], dtype=np.float64)
    for block in groups:
        stratified += x[block].var(0, ddof=0)
    return uniform, stratified / n**2
