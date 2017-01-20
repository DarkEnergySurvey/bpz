import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval


""" =============================
tests on file checking tools  ======
Are the files formatted as we expected?
=================================
To do. Support N-d arrays on pval.mad
"""


def test_required_cols1():
    """extracts the required columns from a test file 1"""

    t = [{'metrics': ['bh_photo_z_validation.sigma_68'], 'truths': 'Z_SPEC', 'weights': 'WEIGHTS', 'predictions': ['Z_MC'], 'bins': [{'Z_MC': '[0, 0.5, 1.0, 2.0]'}]}]

    cols = pval.required_cols(t, 'point')
    corr_cols = ['COADD_OBJECTS_ID', 'Z_SPEC', 'Z_MC', 'WEIGHTS']
    print cols, corr_cols
    for c in corr_cols:

        np.testing.assert_equal(c in cols, True)


def test_required_cols2():
    """extracts the required columns from a test file 2"""

    t = [{'individual': {'metrics': ['bh_photo_z_validation.eval_pdf_point'], 'truths': 'Z_SPEC', 'tolerance': [0.7, 20], 'weights': 'WEIGHTS', 'bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}]}}]

    cols = pval.required_cols(t, 'pdf')

    corr_cols = ['COADD_OBJECTS_ID', 'Z_SPEC', 'MAG_DETMODEL_I', 'WEIGHTS']
    print cols, corr_cols
    for c in corr_cols:
        np.testing.assert_equal(c in cols, True)


def test_required_cols3():
    """extracts the required columns from a test file 3"""

    t = [{'stacks': {'metrics': ['bh_photo_z_validation.kstest', 'bh_photo_z_validation.npoisson', 'bh_photo_z_validation.log_loss'], 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'weights': 'WEIGHTS_LSS', 'truth_bins': [{'Z_SPEC': 'numpy.linspace(0, 2, 4)'}], 'tolerance': [0.7, 20], 'truths': 'Z_SPEC'}}]

    cols = pval.required_cols(t, 'pdf')

    corr_cols = ['COADD_OBJECTS_ID', 'Z_SPEC', 'MAG_DETMODEL_I', 'WEIGHTS_LSS']
    print cols, corr_cols
    for c in corr_cols:
        np.testing.assert_equal(c in cols, True)


def test_required_cols4():
    """extracts the required columns from a test file 4"""

    t = [{'stacks': {'metrics': ['bh_photo_z_validation.kstest', 'bh_photo_z_validation.npoisson', 'bh_photo_z_validation.log_loss'], 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'weights': 'WEIGHTS_LSS', 'truth_bins': [{'Z_SPEC': 'numpy.linspace(0, 2, 4)'}], 'tolerance': [0.7, 20], 'truths': 'Z_SPEC'}}, {'stacks': {'metrics': ['bh_photo_z_validation.kstest', 'bh_photo_z_validation.npoisson', 'bh_photo_z_validation.log_loss'], 'metric_bins': [{'MAG_DETMODEL_I': '[ 17.5, 19, 22, 25]'}], 'weights': 'WEIGHTS_LSS', 'truth_bins': [{'Z_SPEC': 'numpy.linspace(0, 2, 4)'}], 'tolerance': [0.7, 20], 'truths': 'Z_SPEC'}}]

    cols = pval.required_cols(t, 'pdf')

    corr_cols = ['COADD_OBJECTS_ID', 'Z_SPEC', 'MAG_DETMODEL_I', 'WEIGHTS_LSS']
    print cols, corr_cols

    for c in corr_cols:
        np.testing.assert_equal(c in cols, True)


def test_nmad1():
    #http://mathbits.com/MathBits/TISection/Statistics1/MAD.html
    arr = np.array([92, 83, 88, 94, 91, 85, 89, 90])
    res = pval.mean_mad(arr)
    expected = 2.75
    np.testing.assert_equal(res, expected)


def test_mmad2():
    """N-d arrays are not supported"""
    #https://en.wikipedia.org/wiki/Median_absolute_deviation
    arr = np.array([1, 1, 2, 2, 4, 6, 9])
    res = pval.mad(arr)
    expected = 1
    np.testing.assert_equal(res, expected)
