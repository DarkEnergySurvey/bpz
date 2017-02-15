import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval
import vlfn 

""" =============================
tests on metric functions
=================================

"""

def test_delta_sigma_crit_1():
    """Test weak lensing test_delta_sigma_crit metric I"""
    z1 = np.random.normal(size=int(3e5))*0.1 + 1
    z2 = np.random.normal(size=int(3e5))*0.1 + 2
    weights = np.ones(len(z2))
    weights = weights / np.sum(weights)
    res = vlfn.process_function(pval.delta_sigma_crit, z1, z2, weights=weights, extra_params=0.5)

    print res, 0.423893/0.63487
    np.testing.assert_almost_equal(res, 0.423893/0.63487, 2)


def test_delta_sigma_crit_2():
    """Test weak lensing test_delta_sigma_crit metric II"""
    z1 = np.random.normal(size=int(3e5))*0.1 + 0.6
    z2 = np.random.normal(size=int(3e5))*0.1 + 0.6
    weights = np.ones(len(z2))
    weights = weights / np.sum(weights)
    res = vlfn.process_function(pval.delta_sigma_crit, z1, z2, weights=weights, extra_params=0.5)
    print res
    np.testing.assert_almost_equal(res, 1.0, 1)


def test_delta_sigma_crit_3():
    """Test weak lensing test_delta_sigma_crit metric III"""
    z1 = np.random.normal(size=int(3e5))*0.1 + 0.55
    z2 = np.random.normal(size=int(3e5))*0.1 + 0.6
    weights = 1 + z1
    res = vlfn.process_function(pval.delta_sigma_crit, z1, z2, weights=weights, extra_params=0.5)
    np.testing.assert_almost_equal(res, 0.09915666037/0.1399354, 2)

