import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval


""" =============================
tests on hdf5 file format  ======
=================================
"""


def test_hdf1():
    """test get correct error with non-existent hdf5 file """
    filename = 'data/__nonValidHDF__.hdf5'
    cols = ['COADD_OBJECTS_ID', 'Z_SPEC', 'pdf_0']
    err, mess = pval.valid_file(filename, cols)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'file does not exist')


def test_hdf2():
    """test get correct error with non valid hdf5 file """
    filename = 'data/invalidHDF.hdf5'
    cols = ['COADD_OBJECTS_ID', 'Z_SPEC', 'pdf_0']
    err, mess = pval.valid_file(filename, cols)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'missing column COADD_OBJECTS_ID')


def test_hdf4():
    """test can load valid hdf5 file """
    filename = 'data/validHDF.hdf5'
    cols = ['COADD_OBJECTS_ID', 'Z_SPEC']
    err, mess = pval.valid_hdf(filename, cols)
    np.testing.assert_equal(err, True)

    #To do, fix the data type comparison
    #np.testing.assert_equal(type(mess), pd.core.frame.DataFrame)


""" =============================
tests on fits file format  ======
=================================
"""


def test_fits1():
    """test get correct error with non-existent fits file """
    filename = 'data/NotHerevalidPointPrediction.fits'
    cols = ['COADDED_OBJECTS_ID', 'MEAN_Z', 'Z_SPEC']
    err, mess = pval.valid_file(filename, cols)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'file does not exist')


def test_fits2():
    """test get correct error with invalid fits file """
    filename = 'data/invalidPointPrediction.fits'
    cols = ['COADDED_OBJECTS_ID', 'MEAN_Z', 'Z_SPEC']
    err, mess = pval.valid_file(filename, cols)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'missing column MEAN_Z of ' + filename)


def test_fits3():
    """test can load a valid fits file """
    filename = 'data/validPointPrediction.fits'
    cols = ['COADD_OBJECTS_ID', 'MEAN_Z', 'Z_SPEC']
    err, mess = pval.valid_file(filename, cols)
    if err is False:
        print mess
    np.testing.assert_equal(err, True)


""" =============================
tests on point predictions ======
=================================
"""


def test_sigma68_1():
    arr1 = np.arange(1000)
    np.random.shuffle(arr1)
    print 2 * pval.sigma_68(arr1)
    np.testing.assert_equal(np.abs(2 * pval.sigma_68(arr1) - 681.8) < 0.1, True)


def test_sigma68_2():
    arr2 = np.arange(1e4)*1.0
    np.random.shuffle(arr2)
    print 2 * pval.sigma_68(arr2)
    np.testing.assert_equal(np.abs(2 * pval.sigma_68(arr2) - 6824) < 1, True)


def test_sigma68_3():
    arr1 = np.random.normal(size=1e7)
    np.testing.assert_equal(np.abs(pval.sigma_68(arr1) - np.std(arr1)) < 0.01, True)


def test_sigma95_1():
    """test sigma95 to high level of accuracy"""
    arr1 = np.random.normal(size=1e7)
    np.testing.assert_equal(np.abs(pval.sigma_95(arr1) - 2 * np.std(arr1)) < 0.01, True)


def test_sigma95_2():
    """test sigma95 to high intermediate level of accuracy"""
    arr1 = np.arange(1000)
    np.random.shuffle(arr1)
    print 2 * pval.sigma_95(arr1)
    np.testing.assert_equal(np.abs(2 * pval.sigma_95(arr1) - 953) < 1, True)


def test_sigma95_3():
    """test sigma95 to high level of accuracy"""
    arr2 = np.arange(1e4)*1.0
    np.random.shuffle(arr2)
    print 2 * pval.sigma_95(arr2)
    np.testing.assert_equal(np.abs(2 * pval.sigma_95(arr2) - 9540) < 1, True)


def test_delta_z():
    """can we subtract two arrays!"""
    for N in [1, 100, 1000]:
        deltaz = pval.delta_z(np.arange(N), np.arange(N))
        np.testing.assert_equal(np.sum(np.abs(deltaz)), 0)


def test_outlierRate():
    """test outlier rate usig gaussian"""
    #hard coded outlier rate
    arr1 = np.random.normal(size=1e6) * 0.15
    outRate = pval.outlier_rate(arr1)
    np.testing.assert_equal(outRate + 0.68 > 0.99, True)


def test_outlierFraction():
    """test outlier fraction, using gaussian with 0.15"""
    #hard coded outlier rate
    arr1 = np.random.normal(size=1e6) * 0.15
    outRate = pval.outlier_fraction(arr1)
    np.testing.assert_equal(outRate + 68.1 > 99.7, True)


def test_sigma_68_axis1():
    """sigma_68 = std dev for gauss for N-d array along 1 axis"""
    arr = np.zeros((5e5, 10))
    sigs = np.arange(10)*0.1 + 0.5
    for i in range(10):
        arr[:, i] = np.random.normal(size=5e5) * sigs[i] + sigs[i]
    sig68 = pval.sigma_68(arr, axis=0)
    print sig68
    print sigs
    np.testing.assert_almost_equal(sigs, sig68, decimal=2)


def test_sigma_68_axis2():
    """ test sigma_68 = std dev for gauss for N-d array along diff axis"""
    #replicating a 'pdf'
    arr = np.zeros((10, 5e5))
    sigs = np.arange(10)*0.1 + 0.5
    for i in range(10):
        arr[i, :] = np.random.normal(size=5e5) * sigs[i] + sigs[i]
    sig68 = pval.sigma_68(arr, axis=1)
    print sig68
    print sigs
    np.testing.assert_almost_equal(sigs, sig68, decimal=2)


def test_sigma_95_axis():
    """ test sigma_95 = 2xstd dev for gauss for N-d array"""
    arr = np.zeros((5e5, 10))
    sigs = np.arange(10)*0.1 + 0.5
    for i in range(10):
        arr[:, i] = np.random.normal(size=5e5) * sigs[i] + sigs[i]
    sig95 = pval.sigma_95(arr, axis=0)
    print sig95
    print sigs * 2.0
    np.testing.assert_almost_equal(sigs * 2.0, sig95, decimal=2)


""" ===========================
tests on error functions ======
===============================
"""


def test_bootstrap_mean_error1():
    """test error on sigma and sigma"""
    for mean, err in [[1, 1.1], [0.4, 1e-3], [79, 33]]:
        arr = np.random.normal(size=1e5) * err + mean
        weight = np.array([err] * len(arr))
        res = pval.bootstrap_mean_error(arr, weight, np.std)

        #test std and error on std
        np.testing.assert_approx_equal(res['mean'], err, 2)

        errorONerror = err / np.sqrt(2 * len(arr))
        np.testing.assert_approx_equal(res['sigma'], errorONerror, 1)


def test_bootstrap_mean_error2():
    """ test error on mean, and error on mean"""
    for mean, err in [[1, 1.1], [0.4, 1e-3], [79, 33]]:
        arr = np.random.normal(size=1e5) * err + mean
        weight = np.array([err] * len(arr))

        #test for mean
        res = pval.bootstrap_mean_error(arr, weight, np.mean)
        #test mean and std
        errorOnmean = err/np.sqrt(len(arr))
        np.testing.assert_approx_equal(res['mean'], mean, 2)
        np.testing.assert_approx_equal(np.abs(res['sigma'] - errorOnmean) / errorOnmean < 0.1, True)


def test_bootstrap_mean_error_binned():
    """test on binned sigma, mean"""
    bins1 = [0, 2, 10, 30]
    means = [1.1, 5, 20]
    sigs = [1e-3, 0.01, 0.4]
    x = np.array(0)
    for i in range(len(means)):
        x = np.append(x, np.random.normal(size=5e5) * sigs[i] + means[i])
    h = np.histogram(x, bins=500)
    hb = h[1][1:] - (h[1][1] - h[1][1])/2.0

    del x
    weight = h[0]

    bsS = pval.bootstrap_mean_error_binned(hb, hb, weight, bins1, np.std)
    bsM = pval.bootstrap_mean_error_binned(hb, hb, weight, bins1, np.mean)

    print bsS
    print bsM
    np.testing.assert_array_almost_equal(bsM['mean'], means, decimal=1)
    np.testing.assert_array_almost_equal(bsS['mean'], sigs, decimal=2)


""" ==========================
tests two distributions ======
==============================
"""


def test_normalize_pdf():
    """ test can normalise a pdf"""
    N = 1e4
    pdf = np.random.dirichlet(np.arange(N))
    z = np.arange(N, dtype=float) / N
    normPdf = pval.normalize_pdf(pdf, z)
    integ = np.trapz(normPdf, x=z)
    np.testing.assert_equal(np.abs(integ - 1) < 0.001, True)


def test_ks1():
    """ for truth check out http://scistatcalc.blogspot.de/2013/11/kolmogorov-smirnov-test-calculator.html
    """
    arr1 = np.arange(12)
    arr2 = np.arange(12) + 1
    truth = 0.083333
    Pval = 1.0

    ks = pval.ks_test(arr1, arr2)
    np.testing.assert_equal(np.abs(ks - truth) < 0.01, True)

    prob = pval.ks_test_prob(arr1, arr2)
    np.testing.assert_equal(np.abs(prob - Pval) < 0.01, True)


def test_ks2():
    """ for truth check out http://scistatcalc.blogspot.de/2013/11/kolmogorov-smirnov-test-calculator.html
    """
    arr1 = np.arange(300)
    arr2 = np.arange(300)*1.6
    truth = 0.37666
    Pval = 0

    ks = pval.ks_test(arr1, arr2)
    np.testing.assert_equal(np.abs(ks - truth) < 0.01, True)

    prob = pval.ks_test_prob(arr1, arr2)
    np.testing.assert_equal(np.abs(prob - Pval) < 0.01, True)


def test_kulbachLeiber_bins1():
    """KL divergence expects 0"""
    P = np.ones(50, dtype=float)
    Q = np.ones(50, dtype=float)
    KPQ = pval.kulbachLeiber_bins(P, Q)
    expected = 0
    np.testing.assert_equal(expected, KPQ)


def test_kulbachLeiber_bins2():
    """KL divergence expects value """
    P = np.ones(50, dtype=float)
    Q = np.array([np.exp(1)] * len(P))
    KPQ = pval.kulbachLeiber_bins(P, Q)
    expected = -1 * len(P)
    np.testing.assert_equal(expected, KPQ)


def test_kulbachLeiber_bins3():
    """KL tests expects NAN"""
    P = np.ones(50, dtype=float)
    Q = np.array([np.exp(1)] * len(P))
    Q[5] = 0
    KPQ = pval.kulbachLeiber_bins(P, Q)
    expected = np.nan
    np.testing.assert_equal(expected, KPQ)

""" ==================================
tests point distributions  in pdf ====
======================================
"""

def test_eval_pdf_point():
    #make a histogram
    h = np.histogram(np.random.normal(size=5e5), bins=500)
    hb = h[1][1:] - (h[1][1]-h[1][0])/2.0
    h = h[0]
    
    #turn into a pdf
    
    
""" ===========================
tests on useful functions =====
===============================
"""


def test_random_choice1():
    """test can shuffle n-d array alone a dimension"""
    arr = np.arange(55).reshape((11, 5))
    arr2 = pval.random_choice(arr)
    for i in range(5):
        np.testing.assert_equal(np.all(arr[:, i] == arr2), False)


def test_random_choice2():
    """test can shuffle 1-d array"""
    arr = np.arange(5050)
    arr2 = pval.random_choice(arr)
    for i in range(5):
        np.testing.assert_equal(np.all(arr[i] == arr2), False)
