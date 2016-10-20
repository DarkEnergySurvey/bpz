#! /usr/bin/env python
import numpy as np
import sys
import cPickle as pickle
import glob
from lib import Parrallelise
import lib
import bh_photo_z_validation as pval
from astropy.table import Table
import pandas as pd
from scipy.stats import mstats, gaussian_kde

"""
Author Ben Hoyle
Removes all the hassle of obtaining pdfs from a trained machine (using standardised_machine_learning.py with tree methods to generate DES formatted data predictions)

%>des_ml_output.py MLA.p files.* COADD_OBJECTS_ID,REDSHIFT

ToDo:
currently hard coded the parrelisation, using
n_jobs = 3

"""


#loop over all trees to get a redshift dist
def allPredictions(cl, inp):

    pdf = np.zeros((len(inp), len(cl.estimators_)))
    if hasattr(cl, 'estimator_weights_'):
        weights = cl.estimator_weights_ / np.sum(cl.estimator_weights_)
        weightsInt = np.array(len(weights) * 3 / np.sum(weights) * weights, dtype=int)
        wpdf = np.zeros((len(inp), np.sum(weightsInt)))
        del pdf

    sm = 0
    #in case the file is large, split into smaller pieces to save memory.
    if len(inp) > 1e4:
        ind_arr = np.array_split(np.arange(len(inp)), int(len(inp)/1e4) + 2)
    else:
        ind_arr = [np.arange(len(inp))]

    for i, DT in enumerate(cl.estimators_):
        if isinstance(DT, (list, tuple, np.ndarray)):
            DT = DT[0]
        DTz = np.zeros(len(inp))

        for ind_ in ind_arr:
            DTz[ind_] = DT.predict(inp[ind_])

        if hasattr(cl, 'estimator_weights_'):
            #make a weighted pdf
            for j in range(weightsInt[i]):
                wpdf[:, sm] = DTz
                sm += 1
        else:
            pdf[:, i] = DTz

    if hasattr(cl, 'estimator_weights_'):
        return wpdf[:, 0:sm]
    return pdf


def preds_to_pdf_mc(pred_, zbins_):
    pdfs = np.zeros((len(pred_), len(zbins_)), dtype=float)
    z_mc = np.zeros(len(pred_), dtype=float)
    for i in np.arange(len(pred_)):
        kds = gaussian_kde(pred_[i])
        pdfs[i] = kds.evaluate(zbins_)

        z_ = kds.resample(size=1)
        cnt = 0
        while z_ < 0 and cnt < 4:
            z_ = kds.resample(size=1)
            cnt += 1
        if z_ < 0:
            z_ = -1
        z_mc[i] = z_

    return pdfs, z_mc


def makePredictions(testFile):

    # load the test data
    inputs, ExtraColsDict = lib.dataSet(testFile, features['input'], ExtraCols).loadData()

    #remove any nans in feature values.
    ind = np.ones(len(inputs), dtype=bool)
    for i in np.arange(len(inputs[0, :])):
        ind = ind * (np.isfinite(inputs[:, i]) == True)

    if np.sum(ind) < len(ind):
        inputs = inputs[ind]
        for i in ExtraColsDict:
            ExtraColsDict[i] = ExtraColsDict[i][ind]
    else:
        del ind

    MEAN_Z = np.zeros(len(inputs), dtype=float)

    if len(inputs) > 1e4:
        ind_arr = np.array_split(np.arange(len(inputs)), int(len(inputs)/1e4) + 2)
    else:
        ind_arr = [np.arange(len(inputs))]

    for ind_ in ind_arr:
            MEAN_Z[ind_] = res_best['clf'].predict(inputs[ind_])

    # predictions for New Objects
    ExtraColsDict['MEAN_Z'] = MEAN_Z

    #get distributions of all points
    dfs = allPredictions(res_best['clf'], inputs)
    # free space
    del inputs
    print np.shape(dfs)

    MEDIAN_Z = np.zeros(len(dfs), dtype=float)
    for ind_ in ind_arr:
        MEDIAN_Z[ind_] = np.median(dfs[ind_], axis=1)

    ExtraColsDict['MEDIAN_Z'] = MEDIAN_Z
    ExtraColsDict['Z_SIGMA68'] = pval.sigma_68(dfs, axis=1)
    ExtraColsDict['Z_SIGMA'] = np.std(dfs, axis=1)

    binCenters = np.arange(50) / 20.0 + 0.025

    #return the pdf and a random draw from it
    pdfs, z_mc = preds_to_pdf_mc(dfs, binCenters)

    #save the mcmc draw from the pdf
    ExtraColsDict['Z_MC'] = z_mc

    # get the mode
    dfs = np.round_(dfs, decimals=3)
    ExtraColsDict['MODE_Z'] = np.array([mstats.mode(dfs[i], axis=None)[0][0] for i in np.arange(len(dfs))])

    #free space
    del dfs

    #save the point predictions file
    d_ = Table(ExtraColsDict)
    point_file = testFile + mlaTxt + '.DES.point.predictions.fits'
    pdf_file = testFile.replace('.fit', '') + mlaTxt + '.DES.pdf.hdf5'
    d_.write(point_file, format='fits')
    del d_

    for i in ['Z_MC', 'MEAN_Z', 'MEDIAN_Z', 'Z_SIGMA68', 'Z_SIGMA']:
        ExtraColsDict[i] = np.array(ExtraColsDict[i], dtype=np.float16)

    df = pd.DataFrame()
    for i in ExtraColsDict:
        df[i] = ExtraColsDict[i]

    df.to_hdf(pdf_file, 'point_predictions')

    for i in ['Z_MC', 'MEAN_Z', 'MEDIAN_Z', 'Z_SIGMA68', 'Z_SIGMA', 'MODE_Z']:
        del df[i]

    #make smaller data sizes
    pdfs = np.array(pdfs, dtype=np.float16)

    for i, val in enumerate(binCenters):
        df['pdf_%0.2f' % val] = pdfs[:, i]

    df.to_hdf(pdf_file, 'pdf')
    print " "
    print "info written to "
    print pdf_file
    print point_file
    del pdfs

if __name__ == "__main__":
    # get config file name
    args = sys.argv[1:]

    print 'input arguments', args

    if ('-h' in args) or len(args) < 3 or ('help' in args):
        print "incorrect arguments: ", args
        print "des_ml_output.py MLA.p files.* COADD_OBJECTS_ID,REDSHIFT"
        sys.exit()

    mlaFile = args[0]
    ExtraCols = args[-1].split(',')
    files = args[1:-1]
    if isinstance(files, list) is False:
        files = list(files)
        print 'converted to list'

    mlaTxt = mlaFile.split('.p')[0].split('/')[-1]
    # load model
    res_best = pickle.load(open(mlaFile, "r"))
    #we don't need xval results
    if 'results' in res_best:
        del res_best['results']

    # make the features file for the test data
    features = {}
    features['input'] = res_best['features']['input']

    n_jobs = 3

    if n_jobs > len(files):
        n_jobs = len(files)

    if len(files) == 1:
        r = makePredictions(files[0])
    else:
        res1 = Parrallelise(n_jobs=n_jobs, method=makePredictions, loop=files).run()
