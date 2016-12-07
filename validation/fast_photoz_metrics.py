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
import random
import string
from photoz_metrics_fn import perform_tests_fast

"""
Photo-z validation codes

Update Aug 2016, to include WL metrics

Authors: Ben Hoyle

-input:
photoz_metrics.py data/PointPredictions1.fits data/PointPredictions*.fits

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
    truths: REDSHIFT
    
    #should we calculated weighted metrics where available?
    weights: WEIGHTS

    #what metrics do we want to measure. "numpy.std" is the standard deviation from numpy
    # and "bh_photo_z_validation.sigma_68" is the sigma_68 metric found in the bh_photo_z_validation.py file
    metrics_diffz: [numpy.std, numpy.median, bh_photo_z_validation.sigma_68, bh_photo_z_validation.outlier_fraction]
    
    #should we measure some metrics on z_truth and z_predict? e.g. bh_photo_z_validation.wl_metric=|<z1>-<z2>|
    metrics_z1_z2: [bh_photo_z_validation.wl_metric]

    #do we want to assign an accetable tolerance to each of these tests?
    tolerance:
    
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
resultsFilePrefix = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
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

            res[ptype][f][test_name] = perform_tests_fast(d, tst)
            
    #save this output to a file
    with open('point_' + resultsFilePrefix + '.yaml', 'w') as outfile:
        outfile.write(yaml.dump(res[ptype], default_flow_style=False))

    pickle.dump(res[ptype], open('point_' + resultsFilePrefix + '.p', 'w'))
