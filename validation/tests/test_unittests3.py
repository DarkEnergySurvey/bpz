#other test file takes 35, and 30 secs to run. Start a new file
import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval
import time

from weighted_kde import gaussian_kde as gss_kde

"""==============================================================
test we can catch errors in how the test files are constructed ==
=================================================================
"""


def test_valid_tests1():
    """test valid test input"""
    t1 = [{'individual': {'metrics': ['bh_photo_z_validation.eval_pdf_point'], 'truths': 'Z_SPEC', 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'tolerance': [0.7, 20], 'weights': 'WEIGHTS'}, 'stacks': {'metrics': ['bh_photo_z_validation.kstest', 'bh_photo_z_validation.npoisson', 'bh_photo_z_validation.log_loss'], 'tolerance': [0.7, 20], 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'weights': 'WEIGHTS', 'truth_bins': [{'Z_SPEC': 'numpy.linspace(0, 2, 4)'}]}}]
    try:
        r = pval.valid_tests(t1)
        np.testing.assert_equal(True, r)
    except:
        np.testing.assert_equal(True, False)


def test_valid_tests2():
    """test valid test input; a bit ugly unit test. Should die and therefore pass the test"""

    #error here: numpy.linspace(0)
    t1 = [{'metrics': ['numpy.std', 'numpy.median', 'bh_photo_z_validation.sigma_68', 'bh_photo_z_validation.outlier_fraction'], 'weights': 'WEIGHTS', 'error_function': ['bh_photo_z_validation.bootstrap_mean_error'], 'tolerance': [0.4, 0.001, 0.02, 5], 'truths': 'Z_SPEC', 'predictions': ['MODE_Z', 'MEAN_Z', 'Z_MC'], 'bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}, {'MODE_Z': 'numpy.linspace(0)'}]}]
    try:
        pval.valid_tests(t1)
        np.testing.assert_equal(False, True)
    except:
        np.testing.assert_equal(True, True)


def test_valid_tests3():
    """test valid test input; a bit ugly unit test. Should die and therefore pass the test"""

    #error here: '[ 17.5, ,19, 22, 25]'
    t1 = [{'individual': {'metrics': ['bh_photo_z_validation.eval_pdf_point'], 'truths': 'Z_SPEC', 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, ,19, 22, 25]'}], 'tolerance': [0.7, 20], 'weights': 'WEIGHTS'}, 'stacks': {'metrics': ['bh_photo_z_validation.kstest', 'bh_photo_z_validation.npoisson', 'bh_photo_z_validation.log_loss'], 'tolerance': [0.7, 20], 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'weights': 'WEIGHTS', 'truth_bins': [{'Z_SPEC': 'numpy.linspace(0, 2, 4)'}]}}]
    try:
        pval.valid_tests(t1)
        np.testing.assert_equal(False, True)
    except:
        np.testing.assert_equal(True, True)


def test_valid_tests4():
    """test valid test input; a bit ugly unit test. Should die and therefore pass the test"""

    #error here: '[ 17.5, ,19, 22, 25]'
    t1 = [{'metrics': ['numpy.std', 'numpy.median', 'bh_photo_z_validation.sigma_68', 'bh_photo_z_validation.outlier_fraction'], 'weights': 'WEIGHTS', 'error_function': ['bh_photo_z_validation.bootstrap_mean_error'], 'tolerance': [0.4, 0.001, 0.02, 5], 'truths': 'Z_SPEC', 'predictions': ['MODE_Z', 'MEAN_Z', 'Z_MC'], 'bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}, {'MODE_Z': 'numpy.linspace(0, 2, 4)'}]}]
    try:
        pval.valid_tests(t1)
        np.testing.assert_equal(False, True)
    except:
        np.testing.assert_equal(True, True)

