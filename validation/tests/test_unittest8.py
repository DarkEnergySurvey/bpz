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
    """test can recover <z_mc> for fixed gaussian sigma, and fixed mean """
    N = 1000000
    z_arr = np.array([0.6] * N)#np.random.uniform(size=N)*1.5 + 0.4
    z_sigma = np.array([0.1] * N)#np.random.uniform(size=N)*0.1 + 0.1
    x = np.arange(0.1, 3.01, 0.01)

    z_mc = np.zeros(N)
    for i in range(N):
        pdf = mlab.normpdf(x, z_arr[i], z_sigma[i])
        z_mc[i] = pval.get_mc(pdf, x)

    print np.median(z_mc) - np.median(z_arr), (x[1]-x[0]) / 2.0
    np.testing.assert_almost_equal(np.median(z_mc)-np.median(z_arr), 0, decimal=3)


def test_z_mc2():
    """test can get <z_mc> for a range of z_mean, and fixed sigma"""
    N = 100000
    z_arr = np.random.uniform(size=N)*1.5 + 0.4
    z_sigma = np.array([0.1] * N)
    x = np.arange(0.1, 3.01, 0.01)

    z_mc = np.zeros(N)
    for i in range(N):
        pdf = mlab.normpdf(x, z_arr[i], z_sigma[i])
        z_mc[i] = pval.get_mc(pdf, x)

    print np.median(z_mc) - np.average(z_arr, weights=z_sigma), (x[1]-x[0]) / 2.0
    np.testing.assert_almost_equal(np.median(z_mc)- np.average(z_arr, weights=z_sigma), 0, decimal=4)

def test_z_mc3():
    """test get <z_mc> from varying gaussians with random centers and sigmas"""
    N = 1000000
    z_arr = np.random.uniform(size=N)*1.5 + 0.4
    z_sigma = np.random.uniform(size=N)*0.1 + 0.1
    x = np.arange(0.1, 3.01, 0.01)
    dx = (x[1]-x[0])/2.0
    z_mc = np.zeros(N)
    for i in range(N):
        pdf = mlab.normpdf(x, z_arr[i], z_sigma[i])
        z_mc[i] = pval.get_mc(pdf, x)

    print np.mean(z_mc) - np.average(z_arr, weights=z_sigma)
    np.testing.assert_almost_equal(np.median(z_mc) - np.average(z_arr, weights=z_sigma), 0, decimal=3)


def test_z_mc4():
    """test properties of z_mc, is STD correct, is 2*STD correct"""
    N = 3
    z_arr = np.random.uniform(size=N)*1.5 + 0.4
    z_sigma = np.random.uniform(size=N)*0.1 + 0.1
    x = np.arange(0.1, 3.01, 0.01)
    z_mc_std = np.zeros(N)
    z_mc_mean = np.zeros(N)
    z_mc_95 = np.zeros(N)
    for i in range(N):
        pdf = mlab.normpdf(x, z_arr[i], z_sigma[i])
        z_mc = np.array([pval.get_mc(pdf, x) for j in range(100000)])
        z_mc_std[i] = np.std(z_mc)
        z_mc_mean[i] = np.mean(z_mc)
        z_mc_95[i] = pval.sigma_95(z_mc)
    dx = (x[1]-x[0])/2.0
    print 'diff sig', z_sigma - z_mc_std, dx
    print 'diff mean', z_arr - z_mc_mean, dx
    print 'diff 2xsig', z_mc_95 - z_mc_std*2, dx
    for i in range(N):
        np.testing.assert_almost_equal(z_sigma[i]- z_mc_std[i], 0, decimal=3)
        np.testing.assert_almost_equal(z_arr[i] - z_mc_mean[i], 0, decimal=3)
        np.testing.assert_almost_equal(2*z_sigma[i] -z_mc_95[i], 0, decimal=2)
