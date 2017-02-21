"""Validation, function wrappers"""
import bh_photo_z_validation as bhz
import numpy as np
#don't remove double numpy please :D
import numpy as numpy

"""
Authors: Ben Hoyle
"""


def process_function(func, z1, z2, weights=None, extra_params=None):
    """helper function to pass the inputted function func the data it needs to work with select weights data
    z1 - "truth array"
    z2 - prediction array
    weights = set of weights, if required
    extra_params = is whatever func is expecting it to be
    
    each function func manipulates z1, z2 as required internally
    """
    ind = np.arange(len(z1))
    if weights is not None:
        weights_ = weights
        weights_[weights_<0] = 0
        ind = np.random.choice(ind, size=len(z1), p=weights_ / np.sum(weights_), replace=True)

    if extra_params is None:
        return func(z1[ind], z2[ind])
    else:
        #beware this function acts differntly to the others
        return func(z1, z2, weights, extra_params)

def median(z1, z2):
    delta_z = bhz.delta_z(z1, z2)
    return np.median(delta_z)

def median_1pz(z1, z2):
    delta_z_1pz = bhz.delta_z_1pz(z1, z2)
    return np.median(delta_z_1pz)

def sigma_68_1pz(z1, z2):
    delta_z_1pz = bhz.delta_z_1pz(z1, z2)
    return bhz.sigma_68(delta_z_1pz)

def sigma_68(z1, z2):
    delta_z = bhz.delta_z(z1, z2)
    return bhz.sigma_68(delta_z)

def outlier_fraction(z1, z2):
    delta_z = bhz.delta_z(z1, z2)
    return bhz.outlier_fraction(delta_z)

def outFrac_2sigma68(z1, z2):
    delta_z = bhz.delta_z(z1, z2)
    return bhz.outFrac_3sigma68(delta_z)

def outFrac_3sigma68(z1, z2): 
    delta_z = bhz.delta_z(z1, z2)
    return bhz.outFrac_3sigma68(delta_z)

def outFrac_2sigma68_1pz(z1, z2):
    delta_z_1pz = bhz.delta_z_1pz(z1, z2)
    return bhz.outFrac_3sigma68(delta_z_1pz)

def outFrac_3sigma68_1pz(z1, z2):
    delta_z_1pz = bhz.delta_z_1pz(z1, z2)
    return bhz.outFrac_3sigma68(delta_z_1pz)

def wl_metric(z1, z2):
    return bhz.wl_metric(z1, z2)

#beware extra z2weight column.
def delta_sigma_crit(z1, z2, z2weight, z_lens):
    res = {}
    for zl in z_lens:
        res[zl] = bhz.delta_sigma_crit(z1, z2, z2weight, zl)
    return res
