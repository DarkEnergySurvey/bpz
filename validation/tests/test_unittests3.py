#other test file takes 35, and 30 secs to run. Start a new file
import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval
import time

"""==============================================================
test we can catch errors in how the test files are constructed ==
=================================================================
"""


def test_valid_tests1():
    """test valid test input"""
    t1 = [{'pdf': {'individual': {'metrics': ['bh_photo_z_validation.eval_pdf_point'], 'truths': 'Z_SPEC', 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'tolerance': [0.7, 20], 'weights': 'WEIGHTS'}, 'stacks': {'metrics': ['bh_photo_z_validation.kstest', 'bh_photo_z_validation.npoisson', 'bh_photo_z_validation.log_loss'], 'tolerance': [0.7, 20], 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'weights': 'WEIGHTS', 'truth_bins': [{'Z_SPEC': 'numpy.linspace(0, 2, 4)'}]}}, 'test_name': 'photoz-wg', 'point': {'metrics': ['numpy.std', 'numpy.median', 'bh_photo_z_validation.sigma_68', 'bh_photo_z_validation.outlier_fraction'], 'weights': 'WEIGHTS', 'error_function': ['bh_photo_z_validation.bootstrap_mean_error'], 'tolerance': [0.4, 0.001, 0.02, 5], 'truths': 'Z_SPEC', 'predictions': ['MODE_Z', 'MEAN_Z', 'Z_MC'], 'bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}, {'MODE_Z': 'numpy.linspace(0, 2, 4)'}]}}]
    try:
        r = pval.valid_tests(t1)
        np.testing.assert_equal(True, r)
    except:
        np.testing.assert_equal(True, False)


def test_valid_tests2():
    """test valid test input; a bit ugly unit test. Should die and therefore pass the test"""

    #error here: numpy.linspace(0)
    t1 = [{'point': {'metrics': ['numpy.std', 'numpy.median', 'bh_photo_z_validation.sigma_68', 'bh_photo_z_validation.outlier_fraction'], 'weights': 'WEIGHTS', 'error_function': ['bh_photo_z_validation.bootstrap_mean_error'], 'tolerance': [0.4, 0.001, 0.02, 5], 'truths': 'Z_SPEC', 'predictions': ['MODE_Z', 'MEAN_Z', 'Z_MC'], 'bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}, {'MODE_Z': 'numpy.linspace(0)'}]}}]
    try:
        pval.valid_tests(t1)
        np.testing.assert_equal(False, True)
    except:
        np.testing.assert_equal(True, True)


def test_valid_tests3():
    """test valid test input; a bit ugly unit test. Should die and therefore pass the test"""

    #error here: '[ 17.5, ,19, 22, 25]'
    t1 = [{'pdf': {'individual': {'metrics': ['bh_photo_z_validation.eval_pdf_point'], 'truths': 'Z_SPEC', 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, ,19, 22, 25]'}], 'tolerance': [0.7, 20], 'weights': 'WEIGHTS'}, 'stacks': {'metrics': ['bh_photo_z_validation.kstest', 'bh_photo_z_validation.npoisson', 'bh_photo_z_validation.log_loss'], 'tolerance': [0.7, 20], 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'weights': 'WEIGHTS', 'truth_bins': [{'Z_SPEC': 'numpy.linspace(0, 2, 4)'}]}}, 'test_name': 'photoz-wg', 'point': {'metrics': ['numpy.std', 'numpy.median', 'bh_photo_z_validation.sigma_68', 'bh_photo_z_validation.outlier_fraction'], 'weights': 'WEIGHTS', 'error_function': ['bh_photo_z_validation.bootstrap_mean_error'], 'tolerance': [0.4, 0.001, 0.02, 5], 'truths': 'Z_SPEC', 'predictions': ['MODE_Z', 'MEAN_Z', 'Z_MC'], 'bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}, {'MODE_Z': 'numpy.linspace(0, 2, 4)'}]}}]
    try:
        pval.valid_tests(t1)
        np.testing.assert_equal(False, True)
    except:
        np.testing.assert_equal(True, True)


def test_binned_pdf_point_stats():
    """to do, write this test"""
    np.testing.assert_equal(True, True)


def test_xval_cumaltive_at_ypoint():
    """check we can correctly identify the x-axis values at a y-axis point on cdf"""

    ngals = 7
    sigs = np.random.uniform(size=ngals) * 0.1 + 0.02
    x = np.arange(401)/400.0 * 2
    xcentr = x[1:] - (x[1] - x[0]) / 2.0
    arr = np.zeros((ngals, len(xcentr)))

    for i in range(ngals):
        h = np.histogram(np.random.normal(size=1e5) * sigs[i] + 2 * sigs[i], bins=x)
        arr[i, :] = h[0] * 1.0

    median = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.5)
    s2 = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.84075)
    s1 = pval.xval_cumaltive_at_ypoint(arr, xcentr, 0.15825)
    s68 = (s2 - s1) / 2.0

    for i in range(ngals):
        print median[i]*0.5, s68[i], sigs[i]
        np.testing.assert_almost_equal(median[i]*0.5, sigs[i], decimal=2)
        np.testing.assert_almost_equal(s68[i], sigs[i], decimal=2)


def test_npoisson1():
    """test *un normalised dndz in redshfit bins"""
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
    print res
    print xnew2
    for i in range(len(res)):
        np.testing.assert_almost_equal(res[i], xnew2[i], decimal=4)

