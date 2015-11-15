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


