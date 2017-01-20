#! /usr/bin/env python
import sys
import numpy as np
import astropy.io.fits as pyfits
import random as rdm
import copy
from joblib import Parallel, delayed
import bh_photo_z_validation as pval

class Parrallelise:
    """ Creates a generic method to perform
    trivially Parrallel operations
    loop is a tuple of things that one wants to use in the declared function
    call like:
    from bh_parallelise import Parrallelise
    p = Parrallelise(n_jobs=2, method=myFunction, loop=[Tuple,Tuple,..]).run()
    Tuple is passed to myFunction and containes all the info required by the function

    Exampe call:

    #load all the files
    files = glob.glob('/Volumes/Untitled/DES/Y1_GOLD_V1_0_1/*.csv')

    arr = []
    for n, f in enumerate(files):
        arr.append(['ALL_DES_SPECTRA.fits', f])

    res1 = Parrallelise(n_jobs=5, method=stiltsMatch, loop=arr).run()

    """

    def __init__(self, n_jobs=-1, method=None, loop=None):
        self._n_jobs = n_jobs
        self._method = method
        self._loop = loop

    def run(self):
        results = Parallel(n_jobs=self._n_jobs)(delayed(self._method)(l) for l in self._loop)
        return results


def get_mc(pdf, zarr):
    # renorm incase there is probability at higher-z that we've cut, or some error.
    if np.sum(pdf) > 0:
        targ_prob = rdm.random()
        return pval.xval_cumaltive_at_ypoint(pdf, zarr, targ_prob)
    else:
        return -1.

def get_mean_and_sig(pdf, zarr):
    if np.sum(pdf) > 0:
        zm = np.average(zarr, weights=pdf)
        sig = np.sqrt(np.average((zarr-zm)*(zarr-zm), weights=pdf))
        return zm, sig
    else:
        return -1.

def get_median(pdf, zarr):
    if np.sum(pdf) > 0:
        return pval.xval_cumaltive_at_ypoint(pdf, zarr, 0.5)
    else:
        return -1.

def get_sig68(pdf, zarr):
    s2 = pval.xval_cumaltive_at_ypoint(pdf, zarr, 0.84075)
    s1 = pval.xval_cumaltive_at_ypoint(pdf, zarr, 0.15825)
    s68 = (s2 - s1) / 2.0
    return s68

def main_program(in_file):
    print in_file
    file_base = '.'.join(in_file.split('.')[0:-1])
    out_fits = in_file + '.fits'

    #get redshifts from .bpz file
    for j in open(file_base + '.bpz', 'r'):
        t = j.split('=')
        if len(t) > 1:
            if t[0] == '##DZ':
                DZ = float(t[1])
            if t[0] == '##ZMAX':
                ZMAX = float(t[1])
            if t[0] == '##ZMIN':
                ZMIN = float(t[1])
                break

    #slight fix to include ZMAX + DZ
    z_default = np.arange(ZMIN, ZMAX + DZ, DZ)


    pointpred_in = np.loadtxt(file_base + '.bpz')
    mode = copy.copy(pointpred_in[:, 1])
    del pointpred_in

    orig_in = np.loadtxt(file_base + '.cat')
    zspec = copy.copy(orig_in[:, -1])
    del orig_in

    prob_in = np.loadtxt(file_base + '.probs')
    ID = copy.copy(prob_in[:, 0])
    probs = prob_in[:, 1:]

    mean, sigma = np.array([get_meanand_sig(el, z_default) for el in probs])
    median = np.array([get_median(el, z_default) for el in probs])
    mc = np.array([get_mc(el, z_default) for el in probs])
    sig68 = np.array([get_sig68(el, z_default) for el in probs])

    del probs, prob_in

    tbhdu = pyfits.BinTableHDU.from_columns([
                pyfits.Column(name='COADD_OBJECTS_ID', format='K', array=ID),
                pyfits.Column(name='MEAN_Z', format='D', array=mean),
                pyfits.Column(name='MODE_Z', format='D', array=mode),
                pyfits.Column(name='MEDIAN_Z', format='D', array=median),
                pyfits.Column(name='Z_MC', format='D', array=mc),
                pyfits.Column(name='REDSHIFT', format='D', array=zspec),
                pyfits.Column(name='Z_SIGMA68', format='D', array=sig68),
                pyfits.Column(name='Z_SIGMA', format='D', array=sigma)
                                            ])

    tbhdu.writeto(out_fits)

    print ('{:} written to disk' .format(out_fits))
    return True

if __name__ == '__main__':
    args = sys.argv[1:]

    n_jobs = 4
    if (len(args) < 1):
        print ("%>markus_pbz_out_fits.py Paths*ToBPZPredictions.bpz")
        sys.exit(1)

    if isinstance(args, list) is False:
        args = [args]

    res = Parrallelise(n_jobs=n_jobs, method=main_program, loop=args).run()
