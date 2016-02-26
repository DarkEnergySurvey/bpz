#! /usr/bin/env python

"""
Fits file 1-d correlation scripts.
This scripts accepts one (or more) fits files as input and

To-do add unit tests.

Authors: Ben Hoyle (benhoyle1212@gmail.com)
"""

#Identifies properties which correlate with each other.
import sys
import cPickle as pickle
import numpy as np
import bh_photo_z_validation as pval
import minepy as MINE

def show_help():
    """Show the user how to use this file"""
    print "./one_d_correlate.py -h or help to see the help file"
    print "./one_d_correlate.py Path2FitsFile"
    sys.exit()


def correlation_tests(arr1_, arr2_):
    """Calculates a battery of 1-d correlation tests against to input arrays. We remove nans/infs before beginning
    current tests returned as dictionary
    {'MI': The Mutual Information, 'KS': The 1-d KS test D value, 'CC': The Pearsons correlation coefficient}
    """
    ind = np.isfinite(arr1_) * np.isfinite(arr2_)

    #if there are not enough objects to correlate
    if np.sum(ind) < 2:
        return {'MI': np.nan, 'KS': np.nan, 'CC': np.nan}

    #our battery of correlations
    mine = MINE.MINE(alpha=0.6, c=15)
    mine.compute_score(arr1_, arr2_)
    MI_ = mine.mic()
    KS_ = pval.ks_test(arr1_[ind], arr2_[ind])
    CC_ = np.corrcoef(arr1_[ind], arr2_[ind])[0, 1]

    return {'MI': MI_, 'KS': KS_, 'CC': CC_}


if __name__ == "__main__":

    args = sys.argv[1:]

    #show help if asked, or if no input files
    if ('-h' in args) or ('help' in args) or (len(args) == 0):
        show_help()

    from astropy.io import fits
    import yaml
    import glob

    """determine how many files we will be working with"""
    files = []
    if len(args) > 1:
        files = files + glob.glob([i for i in args])
    else:
        files = glob.glob(args[0])

    for dataFile in files:

        out_file_name = 'CorrelationResults_1d_' + dataFile
        hdulist = fits.open(dataFile, memmap=True)
        Nrows = hdulist[1].header['NAXIS2']

        AllCols = hdulist[1].columns
        cols = []
        #only corrlate non-COADD_OBJECTS_ID or non-strings
        for i, c in enumerate(AllCols):
            if ('A' not in c.format) and ('COADD_OBJECTS_ID' not in c.name):
                cols.append(c.name)

        test_results = {}
        for i, c1 in enumerate(cols):
            test_results[c1] = {}
            arr1 = hdulist[1].data[c1]
            for j, c2 in enumerate(cols):
                if (j >= i):
                    arr2 = hdulist[1].data[c2]
                    test_results[c1][c2] = correlation_tests(arr1, arr2)

        res = {'filename': dataFile, 'correlation_results': test_results, 'columns': cols, 'number_rows': Nrows}

        #save results as a pickle file
        pickle.dump(res, open(out_file_name + '.p', 'w'))

        # and? as a yaml file
        with open(out_file_name + '.yaml', 'w') as outfile:
            outfile.write(yaml.dump(res, default_flow_style=False))

        print out_file_name + '.yaml /.p created.'
