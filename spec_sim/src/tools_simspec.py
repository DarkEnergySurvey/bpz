import astropy.io.fits as pf
import numpy as np
import pylab as mplot
from scipy import interpolate
import sys
import dicts
import random
from tools_simspec import *


h = 6.626 * 1e-34            # m^2 kg/s
c = 3 * 1e8                  # m/s

# DES filters, this is only used for normalising the spectra
L_mean = np.array([ 4735.0, 6385.0, 7755.0, 9225.0, 9950.0 ])      # A
delta_L = np.array([ 1470.0, 1410.0, 1470.0, 1470.0, 500.0 ])      # A


def multiply_coeffs(template, coeffs):

    spec = template[0] * 0.0
    tspec = []
    for i in range(len(coeffs)):
        spec += template[i] * coeffs[i]
        tspec.append(template[i] * coeffs[i])
    return spec, tspec

def redshift_spectrum(lam, spec, z):

    # take spectrum and redshift it
    # spectrum has units originally (erg/s/cm^2/A)

    nu = c / (lam * 1e-10 )         # 1/s
    d_lam = lam[1] - lam[0]         # size of bin in A
    lam_min = lam - 0.5 * d_lam     # edges of bin in A
    lam_max = lam + 0.5 * d_lam     # edges of bin in A

    lam_min_shift = lam_min * (1 + z)                   # redshifted edges of bin in A
    lam_max_shift = lam_max * (1 + z)                   # redshifted edges of bin in A
    d_lam_shift = lam_max_shift - lam_min_shift         # redshifted size of bin in A
    lam_shift = (lam_min_shift + lam_max_shift) / 2.0   # redshifted mean of bin in A
    nu_shift = c / (lam_shift * 1e-10)                  # redshifted nu
    
    # shifted spectra, main point is the number of photons in each bin is conserved
    spec_shift = spec / ( 1 + z )**1 
    #(spec * d_lam * 1e-3 / h / nu) / d_lam_shift * h * nu_shift * 1e3

    return lam_shift, spec_shift

def calculate_magnitude(lam, spec, Nfilt):

    # given filter specification, total throughput and a spectrum, calcualte magnitude
    # throughput not included - assumed constant
    # nu lam = c, nu = c / lam, d_nu = c / lam^2 d_lam
    l_mean = L_mean[Nfilt]
    delta_l = delta_L[Nfilt]
    l_min = l_mean - 0.5 * delta_l
    l_max = l_mean + 0.5 * delta_l

    mask = (lam > l_min)*(lam < l_max)
    lam = lam[mask]
    d_lam = lam[1:] - lam[:-1]          # A
    lam = (lam[1:] + lam[:-1]) / 2      # A
    spec = spec[mask]
    spec = (spec[1:] + spec[:-1]) / 2                                                 # erg/s/cm^2/A
    nu = c / (lam * 1e-10 )                                                           # 1/s
    int1 = np.sum(spec * lam**2 / (c*1e10) * d_lam )                                 
    int2 = np.sum(d_lam)  # 1/s=Hz
    # write this again
    mag = -2.5 * np.log10( int1 / int2 ) - 48.6 

    return mag

def calculate_photons(lam, spec, expt, AA):

    d_lam = lam[1:] - lam[:-1]          # A
    lam = (lam[1:] + lam[:-1]) / 2      # A
    spec = (spec[1:] + spec[:-1]) / 2                                                 # erg/s/cm^2/A
    nu = c / (lam * 1e-10 )                                                           
    photons = spec * 10**-7 / h / nu * expt * AA * d_lam

    return photons

def calculate_spec(lam, photon, expt, AA):

    d_lam = lam[1:] - lam[:-1]          # A
    lam = (lam[1:] + lam[:-1]) / 2      # A
    photon = (photon[1:] + photon[:-1]) / 2                                                 # erg/s/cm^2/A
    nu = c / (lam * 1e-10 )                                                           
    spec = photon/10**-7* h * nu / expt / AA / d_lam

    return spec

