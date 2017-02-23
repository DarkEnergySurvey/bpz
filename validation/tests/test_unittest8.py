import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval
import matplotlib.mlab as mlab
import math
""" =============================
tests on metric functions
=================================

"""

def test_z_mc1():

    N = 10000
    z_arr = np.array([0.6] * N)#np.random.uniform(size=N)*1.5 + 0.4
    z_sigma = np.array([0.1] * N)#np.random.uniform(size=N)*0.1 + 0.1
    x = np.arange(0.1, 3.01, 0.01)

    z_mc = np.zeros(N)
    for i in range(N):
        pdf = mlab.normpdf(x, z_arr[i], z_sigma[i])
        z_mc[i] = pval.get_mc(pdf, x)

    print np.median(z_mc) - np.median(z_arr), (x[1]-x[0]) / 2.0
    np.testing.assert_almost_equal(np.mean(z_mc), np.mean(z_arr), decimal=4)


def test_z_mc2():
    N = 1000000
    z_arr = np.random.uniform(size=N)*1.5 + 0.4
    z_sigma = np.array([0.1] * N)#np.random.uniform(size=N)*0.1 + 0.1
    x = np.arange(0.1, 3.01, 0.01)

    z_mc = np.zeros(N)
    for i in range(N):
        pdf = mlab.normpdf(x, z_arr[i], z_sigma[i])
        z_mc[i] = pval.get_mc(pdf, x)

    print np.median(z_mc) - np.median(z_arr), (x[1]-x[0]) / 2.0
    np.testing.assert_almost_equal(np.median(z_mc), np.median(z_arr), decimal=4)

def test_z_mc3():
    N = 1000000
    z_arr = np.random.uniform(size=N)*1.5 + 0.4
    z_sigma = np.random.uniform(size=N)*0.1 + 0.1
    x = np.arange(0.1, 3.01, 0.01)

    z_mc = np.zeros(N)
    for i in range(N):
        pdf = mlab.normpdf(x, z_arr[i], z_sigma[i])
        z_mc[i] = pval.get_mc(pdf, x)

    print np.median(z_mc) - np.median(z_arr), (x[1]-x[0]) / 2.0
    np.testing.assert_almost_equal(np.median(z_mc), np.median(z_arr), decimal=3)
