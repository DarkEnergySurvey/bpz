import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval


"""To do
 add unit tests to pval.binned_statistic_dist1_dist2
 add unit test pval.binned_pdf_point_stats():
 """

"""==============================================================
Tests on distributions ==
=================================================================
"""


def test_xval_cumaltive_at_ypoint():
    """check we can correctly identify the x-axis values at a y-axis point on nD-cdf 1"""

    ngals = 7
    sigs = np.random.uniform(size=ngals) * 0.1 + 0.3
    xcentr = np.linspace(0, 3, 500)
    arr = np.zeros((ngals, len(xcentr)))
    for i in range(ngals):
        arr[i, :] = pval.dist_pdf_weights(np.random.normal(size=100000) * sigs[i] + 5 * sigs[i], xcentr)

    median = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.5)
    print 'median', median

    s2 = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.84075)
    s1 = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.15825)
    s68 = (s2 - s1) / 2.0

    for i in range(ngals):
        print 'values here', median[i]*0.2, s68[i], sigs[i]
        np.testing.assert_almost_equal(median[i]*0.2, sigs[i], decimal=2)
        np.testing.assert_almost_equal(s68[i], sigs[i], decimal=2)


def test_xval_cumaltive_at_ypoint1():
    """check we can correctly identify the x-axis values at a y-axis point on 1-d cdf 2"""

    ngals = 1
    sigs = np.random.uniform() * 0.1 + 0.1
    xcentr = np.linspace(0, 3, 400)
    arr = pval.dist_pdf_weights(np.random.normal(size=100000) * sigs + 4 * sigs, xcentr)

    s1 = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.15825)
    median = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.5)
    s2 = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.84075)

    print 'median', median
    #s1 = pval.xval_cumaltive_at_ypoint(arr, xcentr, )
    s68 = (s2 - s1) / 2.0

    print 'values here', median*0.25, s68, sigs
    np.testing.assert_almost_equal(median*0.25, sigs, decimal=2)
    np.testing.assert_almost_equal(s68, sigs, decimal=2)


def test_npoisson1():
    """test *un* normalised dndz in redshfit bins"""
    arr1 = np.array([25, 25, 25])
    arr2 = np.array([20, 20, 20])
    res = pval.npoisson(arr1, arr2)
    np.testing.assert_equal(res, 1)


def test_interpolate_dist():
    """test correctly coded interpolator"""
    x = np.arange(100)/50.0
    y = np.power(x, 2)
    xnew = np.array([0.1, 1, 1.5])
    xnew2 = xnew * xnew
    res = pval.interpolate_dist(y, x, xnew)
    for i in range(len(res)):
        np.testing.assert_almost_equal(res[i], xnew2[i], decimal=4)
