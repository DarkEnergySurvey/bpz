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
import matplotlib.mlab as mlab

def test_cumaltive_to_point02():
    """check we can we draw a set of z_mc from a pdf, that are flat in histogram heights <h>
    e.g. Bordoloi test for sets of pdfs with 3 bins values """
    from scipy.stats import entropy

    x = np.arange(0.01, 3.01, 0.01)
    dx = (x[1]-x[0]) / 2.0

    for i in range(len(x)-3):
        pdf = np.zeros_like(x)
        pdf[i:i+3] = [0.1, 0.2, 0.1]
        z_mc = pval.get_mc(pdf, x, N=1000000)
        c = pval.cumaltive_to_point(pdf, x, z_mc)
     
        h = np.histogram(c, bins=np.arange(0, 1.05, 0.05))
        print h
        print len(h), 'len(h)'
        h = h[0]
        res = entropy(h, [np.mean(h)]*len(h))
        print 'res', res
        np.testing.assert_array_less(res, 0.005)


def test_cumaltive_to_point01():
    """check we can we draw a set of z_mc from a pdf, that are flat in histogram heights <h>
    e.g. Bordoloi test for sets of pdfs with 2 bins values """
    from scipy.stats import entropy

    x = np.arange(0.01, 3.01, 0.01)
    dx = (x[1]-x[0]) / 2.0

    for i in range(len(x)-2):
        pdf = np.zeros_like(x)
        pdf[i:i+2] = [0.1, 0.2]
        z_mc = pval.get_mc(pdf, x, N=100000)
        c = pval.cumaltive_to_point(pdf, x, z_mc)
     
        h = np.histogram(c, bins=np.arange(0, 1.05, 0.05))
        print h
        print len(h), 'len(h)'
        h = h[0]
        res = entropy(h, [np.mean(h)]*len(h))
        print 'res', res
        np.testing.assert_array_less(res, 0.005)


def test_cumaltive_to_point0():
    """check we can we draw a set of z_mc from a pdf, that are flat in histogram heights <h>
    e.g. Bordoloi test"""
    from scipy.stats import entropy
    ngals = 20
    sigs = np.random.uniform(size=ngals) * 0.5 + 0.1
    x = np.arange(0.01, 3.01, 0.01)
    dx = (x[1]-x[0]) / 2.0
    print 'sigs', sigs
    
    for i in range(ngals):
        pdf = mlab.normpdf(x, sigs[i],  0.02)
        z_mc = pval.get_mc(pdf, x, N=100000)
        c=pval.cumaltive_to_point(pdf, x, z_mc)
        h = np.histogram(c, bins=np.arange(0, 1.05, 0.05))
        print h
        print len(h), 'len(h)'
        h = h[0]
        res = entropy(h, [np.mean(h)]*len(h))
        print 'res', res
        np.testing.assert_array_less(res, 0.005)


def test_xval_cumaltive_at_ypoint():
    """check we can correctly identify the x-axis values at a y-axis point on nD-cdf 1"""

    ngals = 1
    sigs = np.random.uniform(size=ngals) * 0.1 + 0.4
    xcentr = np.arange(0.01, 3.01, 0.01)
    median = np.zeros(ngals)
    s68 = np.zeros(ngals)
    dx = (xcentr[1]-xcentr[0])/2.0/2.0
    for i in range(ngals):
        pdf = mlab.normpdf(xcentr, sigs[i],  sigs[i]*0.02)
        median[i] = pval.xval_cumaltive_at_ypoint(pdf, xcentr, 0.5)
        s68[i] = pval.get_sig68(pdf, xcentr)

    for i in range(ngals):
        print 'values here', i, median[i], s68[i]/0.02, sigs[i], dx
        np.testing.assert_almost_equal(median[i], sigs[i], decimal=4)
        np.testing.assert_almost_equal(s68[i]/0.02, sigs[i], decimal=4)



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
