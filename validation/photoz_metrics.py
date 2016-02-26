#! /usr/bin/env python
import numpy as np
numpy = np
import sys
import os
import yaml
import bh_photo_z_validation as pval
from scipy import stats
import glob
import textwrap
import cPickle as pickle


"""
Photo-z validation codes

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
#test name.
test_name: MyExampleTest1

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
    #which photo-z predictions do we want to test
    predictions: [MODE_Z, MEAN_Z, Z_MC]
    
    #what is the true redshift that we will compare with?
    truths: Z_SPEC
    
    #should we calculated weighted metrics where available?
    weights: WEIGHTS

    #what metrics do we want to measure. "numpy.std" is the standard deviation from numpy
    # and "bh_photo_z_validation.sigma_68" is the sigma_68 metric found in the bh_photo_z_validation.py file
    metrics: [numpy.std, numpy.median, bh_photo_z_validation.sigma_68, bh_photo_z_validation.outlier_fraction]
    
    #do we want to assign an accetable tolerance to each of these tests?
    tolerance: [0.4, 0.001, 0.02, 5]
    
    #Finally do we want to also measure the metrics in some "bins".
    #we define the column_name: 'string of bins / string of function that makes bins'
    bins: [MAG_DETMODEL_I: '[10, 15, 20, 25, 30]', MODE_Z: 'numpy.linspace(0, 2, 20)']

    #Should we calculate errors on each metric? if yes state how
    #you can include as many different error functions as you like. Take care when changing this.
    error_function: [bh_photo_z_validation.bootstrap_mean_error]

#these are the pdf tests
pdf: 
    #we can examine individual redshift pdfs. Remove this part you don't need to compare
    individual:
        truths: Z_SPEC

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
        truths: Z_SPEC
        #we convert truths to a distribution by choosing these bins
        truth_bins: [Z_SPEC: 'numpy.arange(5)*0.33']

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
def get_weights(_dict, _ky, _d):
    #set all objects equal weight, unless defined
    if pval.key_not_none(_dict, _ky) is False:
        print "you have not set any weights for this test"
        print "continuing with weights=1"
        weights = np.ones(len(_d))
    else:
        weights = _d[tst['weights']]

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

#poplate the lists of files for point predictions, and pdf predictions
files = {'point': [], 'pdf': []}

#load the files we will use
for arg in args:
    # are these standard .fits and .hdf5 files?
    files = fileType(arg, files)

    #do we also have a yaml configuration file?
    if '.yaml' in arg:

        config = load_yaml(arg)

        if 'filePaths' in config:
            if pval.key_not_none(config, 'filePaths'):
                for i in config['filePaths']:
                    f = glob.glob(i)
                    for ii in f:
                        files = fileType(ii, files)


if len(files['point']) + len(files['pdf']) < 1:
    print "DES photoz validation code"
    print "usage like"
    print "photoz_metrics.py data/PointPredictions1.fits data/PointPredictions*.fits"
    print "or"
    print "photoz_metrics.py data/pdfPredictions*.hdf5"
    print "or a mix of the two"
    print "photoz_metrics.py data/pdfPredictions*.hdf5 data/PointPredictions*.fits"
    print "or you can make more fine tuned validations using a configuration yaml file"
    print "photoz_metrics.py config.yaml "
    print "an example file has been written to this directory."
    writeExampleConfig()
    sys.exit()


#which sets of metrics + tests shall we perform
testProperties = {'point': [], 'pdf': []}

#if we have loaded config file, then use photo-z + WL metrics
if 'config' in locals():
    #loop over point and pdf
    for ptype in testProperties:
        if pval.key_not_none(config, ptype):
            testProperties[ptype].append(config[ptype])
        if pval.key_not_none(config, 'standardPredictions'):
            for i in config['standardPredictions']:
                p = load_yaml(path + i)

                #now prepare tests from precompile
                if pval.key_not_none(p, ptype):
                    testProperties[ptype].append(p[ptype])
else:
    import string
    import random
    config = {'test_name': 'test_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))}

    #if nothing specified, use the standard tests
    for i in glob.glob(path + '/testConfig/*.yaml'):
        p = load_yaml(i)
        for ptype in testProperties:
            if pval.key_not_none(p, ptype):
                testProperties[ptype].append(p[ptype])

#results file prefix
resultsFilePrefix = ''
if pval.key_not_none(config, 'resultsFilePrefix'):
    resultsFilePrefix = config['resultsFilePrefix']

#results dictionary
res = {}

test_name = None
if pval.key_not_none(config, 'test_name'):
    test_name = config['test_name']

#First point predictions
ptype = 'point'

#do we have any files of this type to work with?
if len(files[ptype]) > 0:
    #results dictionary
    res[ptype] = {}

    #obtain the tests and required cols
    tests = testProperties[ptype]

    #check these test are "valid"
    cont = pval.valid_tests(tests)

    reqcols = pval.required_cols(tests, ptype)

    #loop over all files
    for f in files[ptype]:

        #load a file, and complain if it's not formatted correctly.
        d = load_file(f, reqcols)

        res[ptype][f] = {}

        #calculate all unweighted metrics for deltaz and deltaz/(1+z)
        for tst in tests:

            if test_name is None:
                test_name = 'Test_randid' + str(np.random.randint(0, 1000))
                if pval.key_not_none(tst, 'test_name'):
                    test_name = tst['test_name']

            res[ptype][f][test_name] = {}

            #should we calculate an error on these metrics
            error_function = pval.key_not_none(tst, 'error_function')

            err_metric = {}
            if error_function:
                for ef in tst['error_function']:
                    #turn error function.string into a function
                    err_metric[ef.split('.')[-1]] = pval.get_function(ef)

            for photoz in tst['predictions']:
                res[ptype][f][test_name][photoz] = {}
                diff = pval.delta_z(d[tst['truths']], d[photoz])
                diff_1pz = pval.delta_z_1pz(d[tst['truths']], d[photoz])

                points = {'delta_z': diff, 'diff_1pz': diff_1pz}

                for metric in tst['metrics']:

                    #set all objects equal weight, unless defined
                    weights = get_weights(tst, 'weights', d)

                    res[ptype][f][test_name][photoz][metric] = {}

                    #turn string into function
                    metric_function = pval.get_function(metric)

                    #which residuals shall we employ?
                    for diffpp in points.keys():
                        res[ptype][f][test_name][photoz][metric][diffpp] = {}
                        res[ptype][f][test_name][photoz][metric][diffpp]['VALUE'] = np.asscalar(metric_function(points[diffpp]))

                        #calculate errors on these metrics
                        for ef in err_metric:
                            bstamp_mean_err = err_metric[ef](points[diffpp], weights, metric_function)
                            res[ptype][f][test_name][photoz][metric][diffpp]['MEAN_' + ef] = np.asscalar(bstamp_mean_err['mean'])
                            res[ptype][f][test_name][photoz][metric][diffpp]['SIGMA_' + ef] = np.asscalar(bstamp_mean_err['sigma'])

                        #shall we calculate binning statiscs?
                        if pval.key_not_none(tst, 'bins'):
                            binning = tst['bins']

                            res[ptype][f][test_name][photoz][metric][diffpp]['bins'] = {}
                            for binDict in binning:
                                ky = binDict.keys()[0]
                                bin_vals = eval(binDict[ky])

                                res[ptype][f][test_name][photoz][metric][diffpp]['bins'][ky] = {}
                                #this uses the binned_stats function
                                """http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.binned_statistic.html
                                """

                                #calculate the unweighted statistics in each bin
                                bn_stats = stats.binned_statistic(d[ky], points[diffpp], bins=bin_vals, statistic=metric_function)

                                #determine the center of each bin
                                bn_cntr_sts = stats.binned_statistic(d[ky], d[ky], bins=bin_vals, statistic=np.mean)

                                res[ptype][f][test_name][photoz][metric][diffpp]['bins'][ky]['BIN_CENTERS'] = [np.asscalar(vv) for vv in bn_cntr_sts.statistic]
                                res[ptype][f][test_name][photoz][metric][diffpp]['bins'][ky]['VALUE'] = [np.asscalar(vv) for vv in bn_stats.statistic]

                                #calculate the mean and error by bootstrapping
                                bn_bs_stats = pval.bootstrap_mean_error_binned(d[ky], points[diffpp], weights, bin_vals, metric_function)

                                #calculate the bin 'centers' by boot strapping
                                bn_bs_cnters = pval.bootstrap_mean_error_binned(d[ky], d[ky], weights, bin_vals, np.mean)

                                res[ptype][f][test_name][photoz][metric][diffpp]['bins'][ky]['BIN_CENTERS_MEAN_BS'] = [np.asscalar(vv) for vv in bn_bs_cnters['mean']]

                                res[ptype][f][test_name][photoz][metric][diffpp]['bins'][ky]['BIN_CENTERS_SIGMA_BS'] = [np.asscalar(vv) for vv in bn_bs_cnters['sigma']]

                                res[ptype][f][test_name][photoz][metric][diffpp]['bins'][ky]['MEAN_BS'] = [np.asscalar(vv) for vv in bn_bs_stats['mean']]
                                res[ptype][f][test_name][photoz][metric][diffpp]['bins'][ky]['SIGMA_BS'] = [np.asscalar(vv) for vv in bn_bs_stats['sigma']]

    #save this output to a file
    with open('point_' + resultsFilePrefix + '.yaml', 'w') as outfile:
        outfile.write(yaml.dump(res[ptype], default_flow_style=False))

    pickle.dump(res[ptype], open('point_' + resultsFilePrefix + '.p', 'w'))


""" ==========================
Now compute metrics on pdfs ==
==============================
"""

ptype = 'pdf'

#do we have any files of this type?
if len(files[ptype]) > 0:
    #results dictionary
    res[ptype] = {}

    #obtain the tests and required cols
    tests = testProperties[ptype]

    #check these test are "valid"
    cont = pval.valid_tests(tests)

    reqcols = pval.required_cols(tests, ptype)

    #loop over all files
    for f in files[ptype]:
        d = load_file(f, reqcols)

        res[ptype][f] = {}

        zcols = [c for c in d.keys() if 'pdf_' in c]
        pdf_z_edge = np.array([float(c.split('f_')[-1]) for c in zcols])
        pdf_z_center = pdf_z_edge + np.append((pdf_z_edge[1:] - pdf_z_edge[0: -1]) / 2.0, (pdf_z_edge[-1] - pdf_z_edge[-2]) / 2.0)

        pdf = np.array(d[zcols])

        for tsts in tests:

            "do we have a test name for this code"
            if test_name is None:
                test_name = 'Test_randid' + str(np.random.randint(0, 1000))
                if pval.key_not_none(tst, 'test_name'):
                    test_name = tst['test_name']

            res[ptype][f][test_name] = {}

            if pval.key_not_none(tsts, 'individual'):

                tst = tsts['individual']

                #set standard bins, or use those in the test file
                truths = np.array(d[tst['truths']])
                weights = get_weights(tst, 'weights', d)

                for metric in tst['metrics']:
                    metric_function = get_function(metric)
                    res[ptype][f][test_name][metric] = {}
                    res[ptype][f][test_name][metric]['VALUE'] = np.asscalar(metric_function(pdf, pdf_z_center, truths))

                    #calculate error on statistic
                    if pval.key_not_none(tst, 'error_function'):
                        for errf in tst['error_function']:
                            bserr = get_function(errf)(pdf, pdf_z_center, truths, weights, metric_function)
                            res[ptype][f][test_name][metric]['MEAN_BS' + errf] = np.asscalar(bserr['mean'])
                            res[ptype][f][test_name][metric]['SIGMA_BS' + errf] = np.asscalar(bserr['sigma'])

                        if pval.key_not_none(tests, 'bins'):
                            binning = tests['bins']
                            res[ptype][f][test_name][metric]['binned_result'] = {}
                            for binDict in binning:
                                ky = binDict.keys()[0]
                                ## remove to file testing location
                                bin_vals = eval(binning[ky])
                                data_to_bin = np.array(d[ky])

                                res[ptype][f]['result'][photoz]['binned_result'][ky] = {}
                                res[ptype][f]['result'][photoz]['binned_result'][ky]['bin_column'] = ky
                                res[ptype][f]['result'][photoz]['binned_result'][ky]['bin_values'] = bin_vals

                                binstats = pval.binned_pdf_point_stats(data_to_bin, bin_vals, pdf, pdf_z_center, truths, weights, metric_function)
                                res[ptype][f]['result'][photoz]['binned_result'][ky]['BIN_CENTERS'] = [np.asscalar(binstats[vv]['weighted_bin_center']) for vv in binstats]
                                res[ptype][f]['result'][photoz]['binned_result'][ky]['VALUE'] = [np.asscalar(binstats[vv]['weighted_value']) for vv in binstats]

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
                    res[ptype][f][test_name][metric] = {}
                    res[ptype][f][test_name][metric]['VALUE'] = np.asscalar(func_(truth_dist, stacked_pdf))

                    if pval.key_not_none(tst, 'bins'):

                        binning = tst['bins']
                        res[ptype][f][test_name][metric]['binned_result'] = {}
                        for binDict in binning:
                            bnCol = binDict.keys()[0]
                            bin_vals = eval(binDict[bnCol])

                            res[ptype][f][test_name][metric]['binned_result'][bnCol] = {}
                            res[ptype][f][test_name][metric]['binned_result'][bnCol]['bin_column'] = bnCol

                            #this uses the binned_stats function
                            """http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.binned_statistic.html
                            """
                            binned_stats = pval.binned_statistic_dist1_dist2(np.array(d[bnCol]), bin_vals, truths, pdf, pdf_z_center, func_, weights=weights)

                            res[ptype][f][test_name][metric]['binned_result'][bnCol]['BIN_CENTERS'] = [np.asscalar(binned_stats[vv]['weighted_bin_center']) for vv in binned_stats]

                            res[ptype][f][test_name][metric]['binned_result'][bnCol]['VALUE'] = [np.asscalar(binned_stats[vv]['weighted_value']) for vv in binned_stats]

    #save this output to a file
    with open('pdf_' + resultsFilePrefix + '.yaml', 'w') as outfile:
        outfile.write(yaml.dump(res[ptype], default_flow_style=False))

    pickle.dump(res[ptype], open('pdf_' + resultsFilePrefix + '.p', 'w'))

