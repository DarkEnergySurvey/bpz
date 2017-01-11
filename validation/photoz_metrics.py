#! /usr/bin/env python
import numpy as np
numpy = np
import sys
import os
import yaml
import bh_photo_z_validation as pval

#helper wrapper functions
import vlfn
from scipy import stats
import glob
import textwrap
import cPickle as pickle
import random
import string
import inspect

"""
Photo-z validation codes

Update Aug 2016, to include WL metrics

Authors: Ben Hoyle

-input:
photoz_metrics.py data/PointPredictions1.fits data/PointPredictions*.fits
or
photoz_metrics.py data/pdfPredictions*.hdf5
or a mix of the two
photoz_metrics.py data/pdfPredictions*.hdf5 data/PointPredictions*.fits
or you can make more fine tuned validations using a configuration yaml file
photoz_metrics.py config.yaml

-help
Also see the ipython notebook, called ValidationScriptExample.ipynb
if you run
./photoz_metrics.py
an example configuration file will been written to the directory.

-outputs:
"""

#determine path to enclosing directory
pathname = os.path.dirname(sys.argv[0])
path = os.path.abspath(pathname)


def get_function(function_string):
    import importlib
    module, function = function_string.rsplit('.', 1)
    module = importlib.import_module(module)
    function = getattr(module, function)
    return function


#write a config file, if called without a .fits of .hdf5 file
def writeExampleConfig():
    if os.path.isfile('exampleValidation.yaml') is False:
        f = open('exampleValidation.yaml', 'w')
        txt = textwrap.dedent("""

#paths to file locations. will assume '.fits' as point predictions '.hdf5' as pdf predictions
#add more files to list to compare multiple files
filePaths: ['Y1A1_GOLD101_Y1A1trainValid_14.12.2015.valid.fits', 'Y1A1_GOLD101_Y1A1trainValid_14.12.2015.valid.hdf5']

#1) OPTIONAL Which metrics and tolerance should we measure either a list of metrics, such as
# and or a precomputed collection of group metrics and tolerances
#set blank, or delete this line to not use these preconfigured metrics/bins/tolerances
standardPredictions: [/testConfig/photoz.yaml, /testConfig/weak_lensing.yaml]

# what will the path/ and or/base file name of the results be called?
resultsFilePrefix:

#2) EITHER 1) AND OR OPTIONAL Tests here:
#And or / additionally choose your own metrics, as list
#remove these if not required
#these are the point prediction tests
point:
    #what is the true redshift that we will compare with?
    truths: REDSHIFT

    #shall we use weights when calculating metrics, if so specify here.
    #:  no weights, or list of weights
    weights: [IN_SCIENCE_SAMPLE, WL_WEIGHT]

    #what metrics do we want to measure. 
    #see validation/vlfn.py for a guide of how to add your own.
    metrics: {MEAN_Z: [vlfn.median, vlfn.median_1pz, 
                         vlfn.sigma_68_1pz, vlfn.sigma_68,
                         vlfn.outFrac_2sigma68_1pz, 
                         vlfn.outFrac_3sigma68_1pz,
                         vlfn.outFrac_2sigma68, 
                         vlfn.outFrac_3sigma68], 
              Z_MC : [vlfn.wl_metric, vlfn.delta_sigma_crit]
              }

    #for WL we also care about measuring some quantities as a function of lens Z
    #we need to pass these in as extra parameters
    extra_params: {'vlfn.delta_sigma_crit': [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]}

    #Finally do we want to also measure the metrics in some "bins". 
    #we define the column_name: 'string of bins / string of function that makes bins'
    bins: {MEAN_Z: [0.0, 0.1, 0.2, 0.43, 0.63, 0.9, 1.3]}
    #'[0, 0.1, 0.2, 0.43, 0.63, 0.9, 1.3]']
    #'[0, 0.1, 0.2, 0.39, 0.45, 0.58, 0.75, 1.3]']

    #Should we calculate errors on each metric? if yes state how
    #you can include as many different error functions as you like.
    #e.g. pval.bootstrap_mean_error_binned == boostrap resampled errors
    error_function:
   
#these are the pdf tests
pdf:
    #we can examine individual redshift pdfs. Remove this part you don't need to compare
    individual:
        truths: REDSHIFT

        #let's perform the test found in Bordoloi et al 2012
        metrics: [bh_photo_z_validation.Bordoloi_pdf_test]
        tolerance: [0.7]

        #show we calculate the metric in some user specified bins?
        bins: [MAG_DETMODEL_I: '[ 17.5, 19, 22, 25]']

        #shall we use weights when calculating metrics, if so specify here.
        weights: WEIGHT

        #how will we calculate an error on this test? Take care when changing this
        error_function: [bh_photo_z_validation.bootstrap_mean_error_pdf_point]

    #or shall we compare against stacked pdfs
    stacks:
        truths: REDSHIFT
        #we convert truths to a distribution by choosing these bins
        truth_bins: [REDSHIFT: 'numpy.arange(5)*0.33']

        #which additional bins shall we use to calculate metrics?
        bins: [MAG_DETMODEL_I: '[ 17.5, 19, 22, 25]']
        metrics: [bh_photo_z_validation.ks_test, bh_photo_z_validation.npoisson, bh_photo_z_validation.log_loss]
        tolerance: [0.7, 20, 50]
        #shall we use weights when calculating metrics, if so specify here. e.g. WEIGHTS_LSS
        weights: WEIGHT
""")
        f.write(txt)
        f.close()


def fileType(filename, _dict):
    if '.fit' in filename:
        print '.fits file found, will compute point prediction metrics'
        _dict['point'].append(filename)

    if '.hdf5' in filename:
        print '.hdf5 file found, will computer pdf metrics'
        _dict['pdf'].append(filename)

    return _dict


def load_yaml(filename):
    try:
        d = yaml.load(open(filename, 'r'))
        return d
    except:
        print "error loading yaml file " + filename
        print "check format here http://yaml-online-parser.appspot.com/"
        print "aborting"
        sys.exit()



#get the galaxy weights
def get_weights(_d, _ky):
    #set all objects equal weight, unless defined
    weights = _d[_ky]
    return weights / np.sum(weights)


def load_file(f, cols):
    okay, d = pval.valid_file(f, cols)
    if okay is False:
        print "Aborting because"
        print "error reading file: " + f
        print "error message: " + d
        sys.exit()
    return d

#get input arguments
args = sys.argv[1:]
print args

if '.p' not in args[0] or '.fits' == args[0][-4:]:
    results_file_name = join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5)) + '.p'
    input_files = args
else:
    results_file_name = args[0]
    input_files = args[1:]

if results_file_name[-2:] != '.p':
    results_file_name += '.p'

#poplate the lists of files for point predictions, and pdf predictions
files = {'point': [], 'pdf': []}

#load the files we will use
for arg in input_files:
    # are these standard .fits and .hdf5 files?
    files = fileType(arg, files)

if len(files['point']) + len(files['pdf']) < 1:
    print "DES photoz validation code"
    print "usage like"
    print "photoz_metrics.py ResultsFileName.p data/PointPredictions*.fits"
    print "or"
    print "photoz_metrics.py ResultsFileName.p data/pdfPredictions*.hdf5"
    writeExampleConfig()
    sys.exit()

#which sets of metrics + tests shall we perform
testProperties = {'point': [], 'pdf': []}

import string
import random

#nothing is specified, using the standard tests
test_path = path + '/testConfig/photoz.yaml'
p = load_yaml(test_path)
for ptype in testProperties:
    if pval.key_not_none(p, ptype):
        testProperties[ptype] = p[ptype]

#First point predictions
ptype = 'point'

#do we have any files of this type to work with?
if len(files[ptype]) > 0:

    #results dictionary
    res = {'test_config': testProperties[ptype]}

    #obtain the tests and required cols
    tst = testProperties[ptype]

    #get all the columns we are gonna test on
    reqcols = tst['metrics'].keys()

    #loop over all files
    for f in files[ptype]:

        #load a file, and complain if it's not formatted correctly.
        d = load_file(f, reqcols)

        res[f] = {}

        #which redshift do we need for this metric test?
        for photoz in tst['metrics']:
            res[f][photoz] = {}

            z_truth = np.array(d[tst['truths']])
            z_pred = np.array(d[photoz])

            #what is the metric test?
            for metric in tst['metrics'][photoz]:

                res[f][photoz][metric] = {}

                #convert metric name to metric function
                metric_function = pval.get_function(metric)

                #do I have to pass any additional arguments to this function?
                extra_params = pval.get_extra_params(tst, metric)

                #what weighting scheme shall we apply?
                for wght in tst['weights']:
                    res[f][photoz][metric][wght] = {}

                    #get the data weights
                    weights = np.array(d[wght], dtype=float)

                    res[f][photoz][metric][wght]['value'] = vlfn.process_function(metric_function, z_truth,
                        z_pred, weights=weights, extra_params=extra_params)


                    #shall we calculate binning statiscs?
                    if pval.key_not_none(tst, 'bins'):
                        binning = tst['bins']

                    res[f][photoz][metric][wght]['bins'] = {}
                    for ky in binning:
                        bin_vals = binning[ky]

                        res[f][photoz][metric][wght]['bins'][ky] = {}
                        res[f][photoz][metric][wght]['bins'][ky]['bin_center'] = []
                        res[f][photoz][metric][wght]['bins'][ky]['value'] = []

                        for bbn in range(len(bin_vals)-1):
                            ind_bn = (d[ky] <= bin_vals[bbn + 1]) * (d[ky] > bin_vals[bbn])
                            if np.sum(ind_bn) > 1 and np.sum(weights[ind_bn]) > 0:
                                res[f][photoz][metric][wght]['bins'][ky]['bin_center'].append(np.mean(d[ky][ind_bn]))
                                res[f][photoz][metric][wght]['bins'][ky]['value'].append(vlfn.process_function(metric_function, z_truth[ind_bn], z_pred[ind_bn], weights=weights[ind_bn], extra_params=extra_params))

    pickle.dump(res, open(results_file_name, 'w'))


""" ==========================
Now compute metrics on pdfs ==
==============================
"""

ptype = 'pdf'

#do we have any files of this type?
if len(files[ptype]) > 0:

    res = {'test_config': testProperties[ptype]}
    #obtain the tests and required cols
    tests = testProperties[ptype]

    #check these test are "valid"
    cont = pval.valid_tests(tests)

    reqcols = tests['metrics'].keys()

    #loop over all files
    for f in files[ptype]:
        d = load_file(f, reqcols)

        res[f] = {}

        zcols = [c for c in d.keys() if 'pdf_' in c]
        #pdfs are quoted as bin centers.
        pdf_z_center = np.array([float(c.split('f_')[-1]) for c in zcols])

        pdf = np.array(d[zcols])

        for tsts in tests:

            "do we have a test name for this code"
            if test_name is None:
                test_name = 'Test_randid' + str(np.random.randint(0, 1000))
                if pval.key_not_none(tst, 'test_name'):
                    test_name = tst['test_name']

            res[f] = {}

            if pval.key_not_none(tsts, 'individual'):

                tst = tsts['individual']

                #set standard bins, or use those in the test file
                truths = np.array(d[tst['truths']])
                weights = get_weights(tst, 'weights', d)

                for metric in tst['metrics']:
                    metric_function = get_function(metric)
                    res[f][metric] = {}
                    res[f][metric]['VALUE'] = np.asscalar(metric_function(pdf, pdf_z_center, truths))

                    #calculate error on statistic
                    if pval.key_not_none(tst, 'error_function'):
                        for errf in tst['error_function']:
                            bserr = get_function(errf)(pdf, pdf_z_center, truths, weights, metric_function)
                            res[f][metric]['MEAN_BS' + errf] = np.asscalar(bserr['mean'])
                            res[f][metric]['SIGMA_BS' + errf] = np.asscalar(bserr['sigma'])

                        if pval.key_not_none(tests, 'bins'):
                            binning = tests['bins']
                            res[f][metric]['binned_result'] = {}
                            for ky, bin_vals in binning:
                                ## remove to file testing location
                                data_to_bin = np.array(d[ky])

                                res[f]['result'][photoz]['binned_result'][ky] = {}
                                res[f]['result'][photoz]['binned_result'][ky]['bin_column'] = ky
                                res[f]['result'][photoz]['binned_result'][ky]['bin_values'] = bin_vals

                                binstats = pval.binned_pdf_point_stats(data_to_bin, bin_vals, pdf, pdf_z_center, truths, weights, metric_function)
                                res[f]['result'][photoz]['binned_result'][ky]['BIN_CENTERS'] = [np.asscalar(binstats[vv]['weighted_bin_center']) for vv in binstats]
                                res[f]['result'][photoz]['binned_result'][ky]['VALUE'] = [np.asscalar(binstats[vv]['weighted_value']) for vv in binstats]

                    """ to do, add errors boot strap to this pdf=point binned stats"""

            """Complete pdf - point comparions, now do pdf - pdf comparisons"""
            if pval.key_not_none(tsts, 'stacks'):
                #perform stacks tests
                tst = tsts['stacks']
                truth_col = tst['truths']
                truths = np.array(d[truth_col]).ravel()
                weights = get_weights(tst, 'weights', d).ravel()

                truth_dist = pval.dist_pdf_weights(truths, pdf_z_center, weights=weights)

                if np.any(truth_dist == 0):
                    print 'KDE have some 0 values. This is dangerous!'
                    print truth_dist, pdf_z_center
                    print "aborting"
                    sys.exit()

                #bin centers are defined as the <z> value in the bin, not center of bin.
                #truth_bins_centers = stats.binned_statistic(truths, truths, bins=truth_bins_edges, statistic=np.mean).statistic
                #turn distribution into a pdfs [? remove this ?]
                #truth_pdf = pval.normalisepdfs(truth_dist, truth_bins_centers)

                #stack the pdfs (and re-normalise)
                stacked_pdf = pval.stackpdfs(pdf, weights=weights)
                stacked_pdf = pval.normalisepdfs(stacked_pdf, pdf_z_center)

                for metric in tst['metrics']:
                    func_ = get_function(metric)
                    res[f][metric] = {}
                    res[f][metric]['VALUE'] = np.asscalar(func_(truth_dist, stacked_pdf))

                    if pval.key_not_none(tst, 'bins'):

                        binning = tst['bins']
                        res[f][metric]['binned_result'] = {}
                        for binDict in binning:
                            bnCol = binDict.keys()[0]
                            bin_vals = eval(binDict[bnCol])

                            res[f][metric]['binned_result'][bnCol] = {}
                            res[f][metric]['binned_result'][bnCol]['bin_column'] = bnCol

                            #this uses the binned_stats function
                            """http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.binned_statistic.html
                            """
                            binned_stats = pval.binned_statistic_dist1_dist2(np.array(d[bnCol]), bin_vals, truths, pdf, pdf_z_center, func_, weights=weights)

                            res[f][metric]['binned_result'][bnCol]['BIN_CENTERS'] = [np.asscalar(binned_stats[vv]['weighted_bin_center']) for vv in binned_stats]

                            res[f][metric]['binned_result'][bnCol]['VALUE'] = [np.asscalar(binned_stats[vv]['weighted_value']) for vv in binned_stats]

    pickle.dump(res, open(results_file_name, 'w'))

