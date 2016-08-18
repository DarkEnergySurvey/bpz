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

import numpy as np
from sklearn.metrics import mutual_info_score
import copy
from scipy.stats import entropy


def calc_MI(x, y, bins=None):
    if bins is None:
        bins = 25
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi


def mutual_information(arr1, arr2, bins=None):
    return calc_MI(arr1, arr2, bins=bins)


class MutualInformationCriteria:
    """Determine MIC using scikit-learn"""

    def __init__(self, X, Y):
        self.X, self.Y = X, Y

    def mutual_information_2d(self, bins=None):
        """normalised mutual information score alpha
        http://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html
        """
        return calc_MI(self.X, self.Y, bins=bins) / np.sqrt(calc_MI(self.X, self.X, bins=bins) * calc_MI(self.Y, self.Y, bins=bins))


def mutualInformation(arr1_, arr2_, sample_weight=None):
    r = resample(np.arange(len(arr1_)), sample_weight)
    arr1 = arr1_[r]
    arr2 = arr2_[r]
    return MutualInformationCriteria(arr1, arr2).mutual_information_2d()


def show_help():
    """Show the user how to use this file"""
    print "./one_d_correlate.py -h or help to see the help file"
    print "./one_d_correlate.py file=Path2FitsFile [,max_size=maxNumberOfRowsToRandomSample columns=csv,col,to,sample]"
    sys.exit()


def correlation_tests(arr1_, arr2_, max_size=50000):
    """Calculates a battery of 1-d correlation tests against to input arrays. We remove nans/infs before beginning
    current tests returned as dictionary
    {'MI': The Mutual Information, 'KS': The 1-d KS test D value, 'CC': The Pearsons correlation coefficient}
    """
    ind = np.isfinite(arr1_) * np.isfinite(arr2_)

    #if there are not enough objects to correlate
    if np.sum(ind) < 2:
        return {'MI': np.nan, 'KS': np.nan, 'CC': np.nan}

    if np.sum(ind) > max_size:
        num_to_false = np.sum(ind) - max_size
        indN0 = np.arange(len(ind))[ind]
        np.random.shuffle(indN0)
        indN0 = indN0[0: num_to_false]
        ind[indN0] = False

    #our battery of correlations
    MI_ = mutualInformation(arr1_[ind], arr2_[ind])
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

    argDict = {}
    for i in args:
        if '=' in i:
            ky, vl = i.split('=')
            argDict[ky] = vl
        else:
            if 'files' not in argDict:
                argDict['files'] = i

    if 'max_size' not in argDict:
        argDict['max_size'] = 50000
    files = glob.glob(argDict['files'])
    """determine how many files we will be working with"""

    for dataFile in files:

        out_file_name = 'CorrelationResults_1d_' + dataFile + '_nsample_' + str(argDict['max_size'])
        hdulist = fits.open(dataFile, memmap=True)
        Nrows = hdulist[1].header['NAXIS2']

        AllCols = hdulist[1].columns
        cols2 = []
        #only corrlate non-COADD_OBJECTS_ID or non-strings
        for i, c in enumerate(AllCols):
            if ('A' not in c.format) and ('COADD_OBJECTS_ID' not in c.name):
                cols2.append(c.name)

        if 'columns' in argDict:
            cols1 = argDict['columns'].split(',')
        else:
            cols1 = cols2

        test_results = {}
        for i, c1 in enumerate(cols1):
            test_results[c1] = {}
            arr1 = hdulist[1].data[c1]
            for j, c2 in enumerate(cols2):
                if (j >= i):
                    arr2 = hdulist[1].data[c2]
                    test_results[c1][c2] = correlation_tests(arr1, arr2, max_size=argDict['max_size'])

        res = {'filename': dataFile, 'correlation_results': test_results, 'columns': cols1, 'columns2': cols2, 'number_rows': Nrows}

        pickle.dump(res, open(out_file_name + '.p', 'w'))
        print out_file_name + '.p created'
