import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import numpy as np
import pytest
from control_preserving_emitter import emit_counts, fit_population, population_profiles, project_bulk


def test_identity_is_exact_and_seed_independent():
    raw = np.array([[20, 0, 3], [0, 7, 0], [11, 2, 5]], dtype=np.int32)
    p, b = population_profiles(raw)
    fit, audit = fit_population(raw, p, b)
    assert audit['identity'] and audit['prior_fraction'] == 0
    for seed in [0, 1, 932]:
        assert np.array_equal(emit_counts(raw, fit, seed), raw)


def test_conditional_expectation_activation_depth_and_reproducibility():
    # A zero donor entry must activate. Across 30,000 independent cells the
    # conditional expected counts are known, rather than fitted to outcomes.
    raw = np.tile([30, 20, 0, 50], (30000, 1))
    p = np.tile([.20, .10, .15, .55], (len(raw), 1))
    emitted = emit_counts(raw, p, 1)
    assert np.array_equal(emitted, emit_counts(raw, p, 1))
    assert np.all(emitted.sum(1) == 100)
    assert np.issubdtype(emitted.dtype, np.integer) and emitted.min() >= 0
    assert emitted[:, 2].mean() > 14.8
    np.testing.assert_allclose(emitted.mean(0), p[0]*100, atol=.12, rtol=0)
    assert not np.array_equal(emitted.sum(0), (p[0]*100*len(raw)).astype(int))


def test_nonzero_population_fit_and_zero_limit():
    rng = np.random.default_rng(9)
    raw = rng.poisson(np.linspace(2, 20, 12), (100, 12))
    base = raw / raw.sum(1)[:, None]
    target = base * np.linspace(.9, 1.1, 12)
    target /= target.sum(1)[:, None]
    p, b = target.mean(0), raw.sum(1) @ target / raw.sum()
    fit, audit = fit_population(raw, p, b, tolerance=1e-7, max_fit_l1=1e-5)
    assert not audit['identity']
    np.testing.assert_allclose(fit.mean(0), p, atol=1e-5, rtol=0)
    np.testing.assert_allclose(raw.sum(1) @ fit / raw.sum(), b, atol=1e-5, rtol=0)
    p0, b0 = population_profiles(raw)
    near, _ = fit_population(raw, .999999*p0+.000001*p, .999999*b0+.000001*b)
    assert np.abs(near-base).sum() < 1e-3


@pytest.mark.parametrize('raw,p,kwargs', [
    ([[1, -1]], [[.5, .5]], {}),
    ([[1.5, 1]], [[.5, .5]], {}),
    ([[0, 0]], [[.5, .5]], {}),
    ([[1, 1]], [[.6, .6]], {}),
    ([[3, 1]], [[.5, .5]], {'max_depth': 3}),
    ([[1, 1]], [[.5, .5]], {'max_nnz': 1}),
])
def test_invalid_inputs_and_caps(raw, p, kwargs):
    with pytest.raises(ValueError):
        emit_counts(raw, p, 1, **kwargs)


def test_output_density_cap_is_checked_after_activation():
    with pytest.raises(ValueError, match='Sparse entry'):
        emit_counts([[100, 0, 0]], [[.2, .4, .4]], 2, max_nnz=2)


def test_transport_analytic_variance_and_covariance():
    # Two independent loss processes have mean20 and variance35/3.
    # Redistribution splits those losses .75/.25 between two gain genes.
    n = 24000
    raw = np.tile([30, 20, 0, 50], (n, 1))
    probability = np.tile([.20, .10, .15, .55], (n, 1))
    observed = emit_counts(raw, probability, 20261021).astype(float)
    expected_mean = np.array([20., 10., 15., 55.])
    loss_mean, loss_variance = 20., 35./3
    r = np.array([.75, .25])
    expected_variance = np.r_[20./3, 5., r*(1-r)*loss_mean + r*r*loss_variance]
    assert np.all(abs(observed.mean(0)-expected_mean) < 6*np.sqrt(expected_variance/n))
    centered = observed - expected_mean
    # Six standard errors of the variance estimate, estimated from the fourth
    # centered moment; this is a stochastic-law check, not an arbitrary tolerance.
    variance_se = np.sqrt(np.maximum((centered**4).mean(0)-expected_variance**2, 0)/n)
    assert np.all(abs(observed.var(0)-expected_variance) < 6*variance_se)
    products = centered[:, 2]*centered[:, 3]
    expected_covariance = r.prod()*(loss_variance-loss_mean)
    assert abs(products.mean()-expected_covariance) < 6*products.std(ddof=1)/np.sqrt(n)
    assert np.all(observed.sum(1) == 100)


def test_exact_integer_destination_is_deterministic_at_boundary():
    raw = np.tile([30, 0, 70], (50, 1))
    for seed in (0, 1, 78):
        output = emit_counts(raw, np.tile([0., 1., 0.], (50, 1)), seed)
        np.testing.assert_array_equal(output, np.tile([0, 100, 0], (50, 1)))


def test_nonzero_transport_converges_in_probability_to_identity():
    n = 20000
    raw = np.tile([100, 0], (n, 1))
    # The chance of any edit is exactly 1-(1-epsilon)^100.
    for epsilon in (.0001, .001):
        output = emit_counts(raw, np.tile([1-epsilon, epsilon], (n, 1)), 901)
        changed = (output != raw).any(1).mean()
        expected = 1-(1-epsilon)**100
        se = np.sqrt(expected*(1-expected)/n)
        assert abs(changed-expected) < 6*se
        assert changed < 100*epsilon + 6*se  # Markov bound P(edit)<=E[removed]


def test_equal_depth_population_projects_bulk_to_arithmetic_mean():
    raw = np.array([[40, 60], [60, 40], [20, 80], [80, 20]])
    desired_mean, desired_bulk = np.array([.45, .55]), np.array([.46, .54])
    fit, audit = fit_population(raw, desired_mean, desired_bulk)
    assert audit['projection_l1'] == pytest.approx(.02)
    np.testing.assert_allclose(fit.mean(0), desired_mean, atol=.0002, rtol=0)
    np.testing.assert_allclose(raw.sum(1) @ fit / raw.sum(), fit.mean(0), atol=1e-15, rtol=0)
    with pytest.raises(ValueError, match='projection'):
        fit_population(raw, desired_mean, [.65, .35])


def test_population_support_prior_activates_previously_absent_gene():
    rng = np.random.default_rng(4)
    raw = np.c_[rng.poisson([10, 30, 60], (80, 3)), np.zeros(80, dtype=int)]
    depth = raw.sum(1)
    desired = .98*(raw/depth[:, None])
    desired[:, 3] = .02
    p, b = desired.mean(0), depth @ desired / depth.sum()
    fit, audit = fit_population(raw, p, b)
    assert audit['prior_fraction'] > 0
    assert np.all(fit[:, 3] > 0)
    assert abs(fit[:, 3].mean()-.02) < .0002
    output = emit_counts(raw, fit, 72)
    assert np.count_nonzero(output[:, 3]) > 0
    np.testing.assert_array_equal(output.sum(1), depth)


def test_zero_requested_gene_has_zero_fitted_and_emitted_mass():
    raw = np.array([[10, 10, 10], [10, 30, 20], [20, 10, 10]])
    desired = raw[:, :2]/raw[:, :2].sum(1)[:, None]
    desired = np.c_[desired, np.zeros(len(raw))]
    depth = raw.sum(1)
    fit, _ = fit_population(raw, desired.mean(0), depth @ desired / depth.sum())
    np.testing.assert_array_equal(fit[:, 2], 0)
    np.testing.assert_array_equal(emit_counts(raw, fit, 19)[:, 2], 0)


def test_boundary_control_moments_are_not_shrunk_inward():
    raw = np.array([[100, 0], [0, 200]])
    p, b = population_profiles(raw)
    projected = project_bulk(p, b, raw.sum(1)/raw.sum(1).mean())
    np.testing.assert_array_equal(projected, b)
    fit, audit = fit_population(raw, p, b)
    assert audit['identity'] and audit['projection_l1'] == 0
    np.testing.assert_array_equal(emit_counts(raw, fit, 13), raw)


def test_maximum_depth_identity_and_single_gene_population():
    raw = np.array([[500001, 499999, 0], [1, 1, 999998]], dtype=np.int32)
    p, b = population_profiles(raw)
    fit, audit = fit_population(raw, p, b)
    assert audit['identity']
    np.testing.assert_array_equal(emit_counts(raw, fit, 99), raw)
    single = np.array([[1], [17], [1000000]])
    fit, audit = fit_population(single, [1.], [1.])
    assert audit['identity']
    np.testing.assert_array_equal(emit_counts(single, fit, 10), single)


@pytest.mark.parametrize('limit', [float('nan'), float('inf'), -1.])
def test_nonfinite_or_negative_projection_limit_is_rejected(limit):
    raw = np.array([[40, 60], [60, 40]])
    with pytest.raises(ValueError):
        fit_population(raw, [.45, .55], [.65, .35], max_projection_l1=limit)
