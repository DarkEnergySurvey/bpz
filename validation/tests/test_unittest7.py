import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval

""" =============================
tests on metric functions
=================================

"""

def test_delta_sigma_crit_1():
    """Test weak lensing test_delta_sigma_crit metric I"""
    z1 = np.random.normal(size=int(1e5))*0.1 + 1
    z2 = np.random.normal(size=int(1e5))*0.1 + 2
    weights = np.ones(len(z1))
    res = pval.delta_sigma_crit(z1, z2, weights=weights, z_lens=0.5)
    print res
    np.testing.assert_almost_equal(res, 0.9, 2)


def test_delta_sigma_crit_2():
    """Test weak lensing test_delta_sigma_crit metric II"""
    z1 = np.random.normal(size=int(1e5))*0.1 + 0.6
    z2 = np.random.normal(size=int(1e5))*0.1 + 0.6
    weights = np.ones(len(z1))
    res = pval.delta_sigma_crit(z1, z2, weights=weights, z_lens=0.5)
    print res
    np.testing.assert_almost_equal(res, 0, 2)



def test_delta_sigma_crit_3():
    """Test weak lensing test_delta_sigma_crit metric III"""
    z1 = np.random.normal(size=int(1e5))*0.1 + 0.55
    z2 = np.random.normal(size=int(1e5))*0.1 + 0.6
    weights = 1 + z1
    res = pval.delta_sigma_crit(z1, z2, weights=weights, z_lens=0.5)
    print res
    np.testing.assert_almost_equal(res, 0.05, 2)

