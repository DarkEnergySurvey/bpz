import numpy as np
import matplotlib.pyplot as plt
import math


def wl_metric(z1, z2, weights=None):
    """Determine the WL metric of choice
    |<z1> - <z2>|
    """
    w_ = np.ones(len(z1), dtype=float)
    if weights is not None:
        w_ = weights
    w_ = w_ / np.sum(w_)
    ind = np.random.choice(np.arange(len(w_), dtype=int), size=len(w_), replace=True, p=w_)

    return np.abs(np.mean(z1[ind]) - np.mean(z2[ind]))


def delta_z(z_spec, z_phot):
    return z_spec - z_phot


def delta_z_1pz(z_spec, z_phot):
    return delta_z(z_spec, z_phot) / (1 + z_spec)


def sigma_68(arr, axis=None):
    """Input: an (multi-dimensional) array
    Optional input: the axis along which to calculate the metric
    Outputs: the 68% spread of data about the median value of the array
    """
    upper, lower = np.percentile(arr, [84.075, 15.825], axis=axis)
    return (upper - lower) / 2.0

def sigma_95(arr, axis=None):
    upper, lower = np.percentile(arr, [97.7, 2.3], axis=axis)
    return (upper - lower) / 2.0


def outlier_rate(arr, outR=None):
    """assumes frac outliers >0.15
    """
    if outR is None:
        outR = 0.15
    return np.sum(np.abs(arr) > outR)*1.0/len(arr)


def outlier_fraction(arr):
    return outlier_rate(arr)*1e2


zbins = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
for i in range(len(zbins)-1):
    ind = (z_phot_patch1 > zbins[i]) * (z_phot_patch1 < zbins[i + 1])
    dz = delta_z(z_sim_patch1[ind], z_phot_patch1[ind])

    mu, sig86, sig95, outF = np.median(dz), sigma_68(dz), sigma_95(dz), outlier_fraction(dz)

    wlmetric = wl_metric(np.mean(z_sim_patch1[ind]), np.mean(z_phot_patch1[ind]))

    #save mu, sig86, sig95, outF, wlmetric for each zbin
    #and then compare this with the metrics as measured across the entire footprint
    ind = (z_phot > zbins[i]) * (z_phot < zbins[i + 1])
    dz = delta_z(z_sim[ind], z_phot[ind])
    wlmetric = wl_metric(np.mean(z_sim[ind]), np.mean(z_phot[ind]))

