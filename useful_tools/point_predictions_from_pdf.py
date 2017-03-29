#! /usr/bin/env python
"""
This script loads in a DES redshift-wg redshift pdf prediction file formatted as per: 
https://opensource.ncsa.illinois.edu/bitbucket/projects/DESDM/repos/photoz-wg/

and returns a list of point predictions which may then be massaged before calling the photo-z wg validation script.

call like:

%>point_predictions_from_pdfs.pdf PathToPDF.hdf5 [additional,columns,as,comma,separated,list]

will extract measuremnents from pdfs and generate fits files with ending
PathToPDF.point_predictions.fits

"""

import sys
import numpy as np
import random as rdm
import copy
import bh_photo_z_validation as pval
from joblib import Parallel, delayed
import os
import pandas as pd
from astropy.table import Table


def _help():
    print ("point_predictions_from_pdfs.py PathToPDF.hdf5 [additional,columns,as,comma,separated,list]")
    print ("load PathToPDF.hdf5 which must be in the standard DES format. It will extract measuremnents from pdfs and generate fits files")
    print ("if additional,columns are requested, they must exist in the /point_predictions/ extension")
    sys.exit()


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


def pdf_loop(pdf_, z_bin_centers, ind):
    """extract point predictions from pdfs
    z_bin_centers must be the bin centers
    ind -- currently a dummy index
    """
    n_gals = len(pdf_)
    #results arrays for this loop
    z_mode = np.zeros(n_gals) + np.nan
    mean = np.zeros(n_gals) + np.nan
    sigma = np.zeros(n_gals) + np.nan
    median = np.zeros(n_gals) + np.nan
    mc = np.zeros(n_gals) + np.nan
    sig68 = np.zeros(n_gals) + np.nan

    for i in np.arange(n_gals):
        #define summary stats from the margenalised posterior.
        marg_post = pdf_[i] / np.sum(pdf_[i])
        if np.all(np.isfinite(marg_post)):
            ind_max_marg = np.where(marg_post == np.amax(marg_post))[0]
            if len(ind_max_marg) > 1:
                ind_max_marg = np.random.choice(ind_max_marg)
                ind_max_marg = ind_max_marg[0]
            mean[i] = pval.get_mean(marg_post, z_bin_centers)
            sigma[i] = pval.get_sig(marg_post, z_bin_centers)
            median[i] = pval.get_median(marg_post, z_bin_centers)
            mc[i] = pval.get_mc(marg_post, z_bin_centers)
            sig68[i] = pval.get_sig68(marg_post, z_bin_centers)
            z_mode[i] = z_bin_centers[ind_max_marg]

    return {'index': ind, 'MEAN_Z': mean, 'Z_SIGMA': sigma,
            'MEDIAN_Z': median, 'Z_SIGMA68': sig68, 'MODE_Z': z_mode,
            'Z_MC': mc
            }


def loop_pdf(ind):
    """useful wrapper for potential parrallisation"""
    ind, pdfc = lst
    pdfs = PDF[pdfc][ind].as_matrix()
    return pdf_loop(pdfs, ind)


if __name__ == '__main__':
    args = sys.argv[1:]

    if len(args) < 1 or 'help' in args:
        _help()
    pdf_file = args[0]
    extra_cols = []
    if len(args) > 1:
        extra_cols = args[1].split(',')

    #output file name
    pointp_file = pdf_file.replace('.h5', '.point_predictions.fits')
    if os.path.isfile(pointp_file):
        print (pointp_file, " already exists. Exiting")
        sys.exit()
    #extract the pdfs
    #PDF is a global, for parralelisation pruposes.
    PDF = pd.read_hdf(pdf_file, 'pdf_predictions')

    #get number of rows in file -- useful later for parrallisation...
    pdfc = PDF.keys()[0]
    Nrows = len(PDF[pdfc])

    #from pdf cols of interst
    pdfc = [i for i in PDF.keys() if 'pdf_' in i]

    #sort and extract bin centers.
    pdfc = list(np.sort(pdfc))
    z_bin_centers = [float(i.replace('pdf_', '')) for i in pdfc]

    #get pdf quantaties, ordered by z_bin_centers
    pdfs = PDF[pdfc].as_matrix()

    #measure point prediction quantities
    res = pdf_loop(pdfs, z_bin_centers, None)

    #ignore this temporary index. Useful for parrellisation
    del res['index']

    #now get any additional columns if requested
    df = pd.read_hdf(pdf_file, 'point_predictions')
    for i in extra_cols:
        res[i] = df[i]

    #generate astropy table
    d = Table(res)

    #free space
    del res

    #write results to ouput file
    d.write(pointp_file)
    print (pointp_file, ' written to disk')
