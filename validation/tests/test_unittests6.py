import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval

""" =============================
tests on metric functions
=================================

"""


def test_outFrac_2sigma68():
    """Test outlier fraction rates"""

    arr = np.random.normal(size=int(4e6))
    res = pval.outFrac_2sigma68(arr)
    val = np.sum(np.abs(arr) > 2 * pval.sigma_68(arr)) * 1.0 / len(arr)
    np.testing.assert_almost_equal(res, val, 4)


def test_outFrac_2sigma68v1():
    """Test outlier fraction rates"""

    arr = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 5])
    res = pval.outFrac_2sigma68(arr)
    val = np.sum(np.abs(arr) > 2 * pval.sigma_68(arr)) * 1.0 / len(arr)
    np.testing.assert_almost_equal(res, val, 4)


def test_outFrac_3sigma68():
    """Test outlier fraction rates"""

    arr = np.random.normal(size=int(4e6))
    res = pval.outFrac_3sigma68(arr)
    val = np.sum(np.abs(arr) > 3 * pval.sigma_68(arr)) * 1.0 / len(arr)
    np.testing.assert_almost_equal(res, val, 4)


def test_outFrac_3sigma68v1():
    """Test outlier fraction rates"""

    arr = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 5])
    res = pval.outFrac_3sigma68(arr)
    val = np.sum(np.abs(arr) > 3 * pval.sigma_68(arr)) * 1.0 / len(arr)
    np.testing.assert_almost_equal(res, val, 4)


def test_wl_metric():
    """Test weak lensing metric"""
    z1 = np.random.normal(size=int(1e5)) + 1
    z2 = np.random.normal(size=int(1e5)) + 2
    res = pval.wl_metric(z1, z2)
    np.testing.assert_almost_equal(res, 1, 2)