import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import itertools
import numpy as np
import pytest
from stratified_donor_sampling import depth_blocks, select_donors, exact_mean_sampling_variance


def test_exact_partition_stable_ties_and_reproducibility():
    depth = np.array([4, 1, 2, 1, 3, 4, 2, 3])
    blocks = depth_blocks(depth, 4)
    np.testing.assert_array_equal(blocks, [[1,3],[2,6],[4,7],[0,5]])
    for seed in [1, 2, 303]:
        rows = select_donors(depth, seed, 'depth_stratified', 4)
        assert len(set(rows)) == 4
        assert all(row in block for row,block in zip(rows,blocks))
        np.testing.assert_array_equal(rows, select_donors(depth,seed,'depth_stratified',4))


def test_exact_enumeration_unbiasedness_and_variance():
    # Enumerate every possible outcome rather than rely on a loose MC tolerance.
    x = np.array([[1.,10.],[2.,0.],[8.,12.],[11.,-3.],[20.,4.],[23.,1.]])
    blocks = depth_blocks(np.arange(1,7), 3)
    unif = np.array([x[list(rows)].mean(0) for rows in itertools.combinations(range(6),3)])
    strat = np.array([x[list(rows)].mean(0) for rows in itertools.product(*blocks)])
    uvar,svar = exact_mean_sampling_variance(x,blocks)
    np.testing.assert_allclose(unif.mean(0),x.mean(0))
    np.testing.assert_allclose(strat.mean(0),x.mean(0))
    np.testing.assert_allclose(unif.var(0),uvar)
    np.testing.assert_allclose(strat.var(0),svar)
    assert svar[0] < uvar[0]
    # Inclusion counts are exactly equal across all possible stratified samples.
    selections = np.array(list(itertools.product(*blocks)))
    assert np.all(np.bincount(selections.ravel(),minlength=6)==len(selections)/2)


def test_uniform_matches_existing_generator_draw():
    d=np.arange(1,4001)
    np.testing.assert_array_equal(select_donors(d,20261021,'uniform'),
                                 np.random.default_rng(20261021).choice(4000,400,replace=False))


@pytest.mark.parametrize('depth,cells', [([1,2,3],2),([1,0],1),([1,np.nan],1),([1,2.5],1),([[1,2]],1),([1,2],0),([1,2],3)])
def test_bad_design_rejected(depth,cells):
    with pytest.raises(ValueError):depth_blocks(depth,cells)


def test_unknown_policy_rejected():
    with pytest.raises(ValueError):select_donors([1,2],1,'adaptive',1)


def test_bad_partition_rejected():
    with pytest.raises(ValueError):exact_mean_sampling_variance(np.ones((4,2)),np.array([[0,1],[1,2]]))
