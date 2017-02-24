#other test file takes 35 secs to run. Start a new file
import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval


""" =============================
tests on pdfs tools  ======
=================================
"""


def test_stackpdfs():
    """can we stack a pdf over many galaxies"""
    ngals = 56
    pdfs = np.zeros((ngals, 500))
    for i in range(ngals):
        pdfs[i, :] = i

    res = pval.stackpdfs(pdfs)

    #what is the expected result? (sum of 0,1,..55 in each bin)
    exp = np.array([np.asscalar(np.sum(np.arange(ngals)))] * 500)
    np.testing.assert_array_equal(res, exp)


def test_normalisepdfs1():
    """ can we renomalised a mutli-d pdf along each axis"""
    ngals = 56
    pdfs = np.zeros((ngals, 500))
    for i in range(ngals):
        pdfs[i, :] = i * i + 1

    x = np.arange(500)*0.1
    npdfs = pval.normalisepdfs(pdfs, x)

    for i in range(ngals):
        np.testing.assert_almost_equal(np.trapz(npdfs[i], x), 1, 4)


def test_normalisepdfs2():
    """ can we renomalised a 1-d pdf along one axis"""
    pdf = np.arange(500)

    x = np.arange(500)*0.1
    npdfs = pval.normalisepdfs(pdf, x)

    np.testing.assert_almost_equal(np.trapz(npdfs, x), 1, 4)


def test_integrate_dist_bin1():
    """ can we integrate a 1-d pdf in a specified bin"""

    pdf = np.array([1] * 500)
    x = np.arange(500) * 0.1
    minval = 4
    maxval = 24
    #shoud == maxval-minval
    tot = pval.integrate_dist_bin(pdf, x, minval, maxval)
    np.testing.assert_almost_equal(tot, 20, 4)


def test_integrate_dist_bin2():
    """ can we integrate a n-d pdf in a bin"""

    ngals = 56
    pdfs = np.zeros((ngals, 500))
    for i in range(ngals):
        pdfs[i, :] = i + 1

    x = np.arange(500) * 0.1
    minval = 4
    maxval = 24
    #shoud == (maxval-minval) * value in pdfs[i, 0]
    tot = pval.integrate_dist_bin(pdfs, x, minval, maxval)

    for i in range(ngals):
        np.testing.assert_almost_equal(tot[i], (maxval-minval) * (i + 1), 4)

def test_cumaltive_to_point1():
    """ can we determine cumulative 1-d pdf at a set of points"""

    #generate a fake df, flat across the 500 bins
    pdf = np.ones(500, dtype=float)
    pdf = pdf/np.sum(pdf)
    #calcalate cumulative df of this, up to 0, 1, 2, 3,.. ngals etc
    for i in np.arange(40)+1:
        res = pval.cumaltive_to_point(pdf, np.arange(500)+0.5, i)
        print i, np.sum(pdf[0:i])/500.0, res
        np.testing.assert_almost_equal(res, np.sum(pdf[0:i])/500.0, 3)


def test_cumaltive_to_point2():
    """ can we determine cumulative n-d pdf at a set of points"""

    #generate a fake df, flat across the 500 bins
    ngals = 56
    pdfs = np.zeros((ngals, 500))
    binCenters = np.arange(500)+0.5
    for i in np.arange(ngals):
        pdfs[i, :] = i + 1

    #npfds = pval.normalisepdfs(pdfs, binCenters)

    #calcalate cumulative df of this, up to 0, 1, 2, 3,.. ngals etc
    res = pval.cumaltive_to_point(npfds, binCenters, np.arange(ngals)+0.5)

    for i in np.arange(ngals-1) + 1:
        print i, res[i], np.sum(npfds[i, 0:i+1])
        np.testing.assert_almost_equal(res[i], np.sum(npfds[i, 0:i+1]), 1)


def test_gini_criteria():
    """ test gini codes, using values found
    https://en.wikipedia.org/wiki/Gini_coefficient
    """
    np.testing.assert_equal(pval.gini(np.ones(10)), 0)
    np.testing.assert_almost_equal(pval.gini(np.sqrt(np.arange(100)/101.0)), 0.2, 2)
    np.testing.assert_almost_equal(pval.gini(np.power(np.arange(100)/101.0, 2)), 0.5, 2)
    np.testing.assert_almost_equal(pval.gini(np.power(np.arange(100)/101.0, 3)), 0.6, 2)


def test_dfs_mode1():
    """test we can extract the mode from an N-d distribution"""
    pdfs = np.array([[0, 1, 12, 3, 5], [10, 9, 8, 7, 0], [0, 1, 2, 3, 4]])
    print np.shape(pdfs)
    x = np.arange(5)
    mds = pval.dfs_mode(pdfs, x)
    np.testing.assert_array_equal(mds, [2, 0, 4])


def test_dfs_mode2():
    """test we can extract the mode from a distribution"""
    pdfs = np.array([0, 1, 12, 3, 5])
    x = np.arange(5)
    mds = pval.dfs_mode(pdfs, x)
    np.testing.assert_equal(mds, 2)


def test_Bordoloi_pdf():
    """test the cum pdf distribution upto a point estimate is flat between 0-1 for a heap of samples.
    We use the gini criteria to check for equivilence"""

    #make a normalised pdf
    h = np.histogram(np.random.normal(size=5e6)*0.1 + 2, bins=800)
    dist = h[0]
    bn = h[1][1:] - (h[1][1] - h[1][0]) / 2.0
    pdf = pval.normalisepdfs(dist, bn)

    ngals = 1e5
    pdfs = np.tile(pdf, [ngals, 1])

    specz = np.random.normal(size=ngals)*0.1 + 2

    gini = pval.Bordoloi_pdf_test(pdfs, bn, specz)

    print gini
    np.testing.assert_almost_equal(gini, 0, decimal=1)


