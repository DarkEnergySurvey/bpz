#! /usr/bin/env python

"""make pixelised maps of column properties, for each file at a user specified resolution"""
"""Author: Ben Hoyle """
"""v0.1 date: 13th May 2016"""
"""Requires, scipy, numpy polspice 'spice'"""

import math
import numpy as np
import sys
#leave double import please!
import numpy as numpy
import healpy as hp
from astropy.io import fits
from scipy.stats import binned_statistic
import os
import shutil
from spice import spice


def poissonNumberToDensity(numberMap, mask):
    #convert discreetized map back to density map
    meanInMask = np.mean(numberMap[mask])
    densityMap = (numberMap - meanInMask) / (meanInMask)
    return densityMap, meanInMask


#clnoise is an large object with many poisson samples of sky
def theoryErrorsNbar(cl, clnoise, fsky):
    l = np.arange(len(cl))
    #deltaCl' = delta(Cl+1/n) = sqrt(2/(2l+1)*1/fsky)* Cl'
    #error propagation formulea
    deltaCl = np.sqrt(2.0 / ((2 * l + 1) * fsky)) * cl
    deltaCl2 = deltaCl**2
    #DeltaN
    deltaClN = [np.sqrt(2.0/((2 * l + 1) * fsky)) * cli for cli in clnoise]
    Nbar = np.mean(deltaClN, axis=0) + np.std(deltaClN, axis=0)
    #print (NbarErr/Nbar)**2
    #x,y     = np.shape(deltaClN)
    #sig681 = [sig68(np.array([deltaClN[j][i] for j in np.arange(x)])) for i in np.arange(y)]
    #print np.ptp(deltaClN,axis=0)/Nbar
    print deltaCl / Nbar
    #print Nbar[0:10]
    #print (NbarErr**2/Nbar**4)[0:10]
    #deltaN2 = 1.0/(Nbar**4) * (NbarErr**2)
    #print deltaN2
    #get combined error for each possion sampling
    ClNerr = np.sqrt(Nbar**2 + deltaCl2)
    return ClNerr


#get weighted mean and error for binned datapoints
#lbin correspond to min and max of bins
def wbinCls(l, cl, err, lbin):
    #binned_statistics
    Clbin = np.zeros(len(lbin) - 1)
    Clberr = np.zeros(len(lbin) - 1)
    lmean = np.zeros(len(lbin) - 1)
    for i in range(len(lbin) - 1):
        ind = np.where((l >= lbin[i]) * (l < lbin[i + 1]))[0]
        if len(ind) > 0:
            w = (1.0 / err[ind]) ** 2
            Clbin[i] = np.sum(cl[ind] * w) / np.sum(w)
            Clberr[i] = np.sqrt(1.0 / np.sum(w))
            lmean[i] = np.mean(l[ind])
    return Clbin, Clberr, lmean


def spice_cls(fileName1, fileName2, maskfile, lmax, thetamax=10):
    use_weights = True
    beam = 'NO'
    pixelfile = 'YES'
    weightfile = 'NO'
    apodizesigma = thetamax
    covfileout = "YES"
    print fileName1, fileName2, maskfile, lmax, thetamax
    cl_spice = spice(bin=False, norm=False, mapfile=fileName1, mapfile2=fileName2, nlmax=lmax,
                     polarization="NO", beam_file="NO", beam_file2="NO", pixelfile=pixelfile,
                     weightfile=weightfile, apodizesigma=apodizesigma, thetamax=thetamax,
                     covfileout=covfileout, maskfile2=maskfile, maskfile=maskfile
                     )

    clc5_s = (np.array(cl_spice['TT']))[0: lmax]

    covF = pyfits.open('spice.covariance.fits')
    cov = covF[0].data[0]
    di = np.diag_indices(len(cov))
    err_s = np.sqrt(cov[di])[1: lmax + 1]

    return clc5_s, err_s


def err_message(args):
    print "two_d_correlate.py pathTohealPix/FitFiles*.fits [, mask=maskFile, thetamax=thetamax, output_file=fileName]"


if __name__ == "__main__":
    import glob

    args = sys.argv[1:]
    if 'help' in args or len(args) < 1:
        err_message(args)
        sys.exit()

    inArgs = {}
    for i in args:
        if '='in i:
            k, v = i.split('=')
            inArgs[k] = v

    fils = [i for i in args if '=' not in i]

    #copy mask file to a location + name that PolSpice likes
    mask_file = 'total_mask.fits'
    if 'mask' in inArgs:
        shutil.copy(inArgs['mask'], 'total_mask.fits')
    
    thetamax = 10
    if 'thetamax' in inArgs:
        thetamax = float(inArgs['thetamax'])

    if 'output_file' in inArgs:
        output_file = inArgs['output_file']

    mask = hp.read_map(mask_file)

    #lmax is 1.5x Nside
    n_side = int(np.sqrt(len(mask) / 12.0))
    lmax = int(n_side * 1.5)

    output_file = 'outFile.Nside_' + str(n_side) + '.p'
    if 'output_file' in inArgs:
        output_file = inArgs['output_file']

    corr_res = {
        'files': fils,
        'L': np.arange(lmax),
        'mask_file': mask_file,
        'Nside': n_side
        }

    for fil1 in fils:
        corr_res[fil1] = {}
        for fil2 in fils:
            corr_res[fil1][fil2] = {}

            #call polspice on written map #spiceCls
            clT1, clT1err = spice_cls(fil1, fil2, mask_file, lmax, thetamax=thetamax)
            corr_res[fil1][fil2]['Cl'] = clT1
            corr_res[fil1][fil2]['Cl_err'] = clT1err
            pickle.dump({'correlation_results': corr_res, 'files1': fils, 'files2': fils}, open(output_file, 'w') )

    print output_file
