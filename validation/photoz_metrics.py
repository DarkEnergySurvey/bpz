#! /usr/bin/env python
import numpy as np
import sys
import os
import yaml
import sml_mla as ml
import bh_photo_z_validation as pval
import photo_z_metrics as pzm
from scipy import stats
import glob
import textwrap
import inspect

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
        txt = textwrap.dedent("""#paths to file locations. will assume '.fits' as point predictions '.hdf5' as pdf predictions
#add more files to list to compare multiple files
filePaths: ['tests/data/validPointPrediction.fits', 'tests/data/validHDF.hdf5']

#Which metrics and tolerance should we measure either a list of metrics, such as
# and or a precomputed collection of group metrics and tolerances
#set blank, or delete this line to not use these preconfigured metrics/bins/tolerances
standardPredictions: [/testConfig/photoz.yaml, /testConfig/weak_lensing.yaml]

#And or / additionally choose your own metrics, as list
#remove these if not required
point:
    predictions: ['MODE_Z', 'MEAN_Z', 'Z_MC']
    truths: 'Z_SPEC'
    bins: [Z_SPEC: 'np.linspace(0,3,10)', MAG_DETMODEL_I: '[10, 15, 20, 25, 30]']
    metrics: [numpy.std, numpy.median, bh_photo_z_validation.sigma_68, bh_photo_z_validation.outlier_fraction]
    tolerance:


#these are the pdf tests
pdf:
    bins: [Z_SPEC: 'np.linspace(0,3,10)', I_MAG: '[10, 15, 20, 25, 30]']
    metrics: [bh_photo_z_validation.kstest, bh_photo_z_validation.npoisson, bh_photo_z_validation.log_loss]
    tolerance:
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


#check the existence and non-nullness of a key
def key_not_none(_dict, ky):
    if ky in _dict:
        if _dict[ky] is not None:
            return True
    return False


def load_yaml(filename):

    try:
        d = yaml.load(open(filename, 'r'))
        return d

    except:
        print "error loading yaml file " + filename
        print "check format here http://yaml-online-parser.appspot.com/"
        print "aborting"
        sys.exit()


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
            if key_not_none(config, 'filePaths'):
                for i in config['filePaths']:
                    files = fileType(i, files)


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
        if key_not_none(config, ptype):
            testProperties[ptype].append(config[ptype])
        if key_not_none(config, 'standardPredictions'):
            for i in config['standardPredictions']:
                p = load_yaml(path + i)

                #now prepare tests from precompile
                if key_not_none(p, ptype):
                    testProperties[ptype].append(p[ptype])
else:
    #if nothing specified, use the standard tests
    for i in glob.glob(path + '/testConfig/*.yaml'):
        p = load_yaml(i)
        for ptype in testProperties:
            if key_not_none(p, ptype):
                testProperties[ptype].append(p[ptype])

#results dictionary
res = {}

#First point predictions
ptype = 'point'
#do we have any files of this type?
if len(files[ptype]) > 0:
    #results dictionary
    res[ptype] = {}

    tests = testProperties[ptype]
    #loop over all files
    for f in files[ptype]:
        okay, d = pval.valid_file(f)
        if okay is False:
            print "Aborting because"
            print "error reading file: " + f
            print "error message: " + d
            sys.exit()

        res[ptype][f] = {}

        #calculate all unweighted metrics for deltaz and deltaz/(1+z)
        for testNum, tst in enumerate(tests):
            res[ptype][f][testNum] = {}

            for photoz in tst['predictions']:
                res[ptype][f][testNum][photoz] = {}
                diff = pval.delta_z(d[tst['truths']], d[photoz])
                diff_1pz = pval.delta_z_1pz(d[tst['truths']], d[photoz])

                points = {'delta_z': diff, 'diff_1pz': diff_1pz}

                for metric in tst['metrics']:
                    res[ptype][f][testNum][photoz][metric] = {}

                    for diffpp in points.keys():
                        res[ptype][f][testNum][photoz][metric][diffpp] = {}
                        res[ptype][f][testNum][photoz][metric][diffpp]['global'] = get_function(metric)(points[diffpp])

                        if key_not_none(tst, 'bins'):
                            binning = tst['bins']

                            res[ptype][f][testNum][photoz][metric][diffpp]['bins'] = {}
                            for binDict in binning:
                                ky = binDict.keys()[0]
                                if ky not in d.keys():
                                    print "You asked to bin in " + ky
                                    print "but it does not exist in file " + f
                                    print "aborting"
                                    sys.exit()
                                try:
                                    bin_vals = eval(binDict[ky])
                                except:
                                    print "unable to build the bins, please check syntax: " + binDict[ky]
                                    print "Aborting"
                                    sys.exit()

                                res[ptype][f][testNum][photoz][metric][diffpp]['bins'][ky] = {}
                                #this uses the binned_stats function
                                """http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.binned_statistic.html
                                """
                                res[ptype][f][testNum][photoz][metric][diffpp]['bins'][ky] = stats.binned_statistic(d[ky], points[diffpp], bins=bin_vals, statistic=get_function(metric))


    #now print out results in a nice format
    print "filename,testSample,ReshiftPointEstimate,metric,redidualOrRedshiftScaled,value"
    for f in res[ptype]:
        for testNum in res[ptype][f]:
            for photoz in res[ptype][f][testNum]:
                for metric in res[ptype][f][testNum][photoz]:
                    for diffp in res[ptype][f][testNum][photoz][metric]:
                        print f + ',' + str(testNum) + ',' + photoz + ',' + metric + ',' + diffp + ',' + str(res[ptype][f][testNum][photoz][metric][diffp]['global'])

    print "filename,testSample,ReshiftPointEstimate,metric,redidualOrRedshiftScaled,binColumn,"
    for f in res[ptype]:
        for testNum in res[ptype][f]:
            for photoz in res[ptype][f][testNum]:
                for metric in res[ptype][f][testNum][photoz]:
                    for diffp in res[ptype][f][testNum][photoz][metric]:
                        for binCols in res[ptype][f][testNum][photoz][metric][diffp]['bins']:
                            g = ''#print f + ',' + str(testNum) + ',' + photoz + ',' + metric + ',' + diffp + ',' + str(binCols) + ',' + str(res[ptype][f][testNum][photoz][metric][diffp]['bins'][binCols])

""" 
To do. work with pdfs
"""

1/0

#I don't know how to deal with pdf metrics yet
ptype = 'pdf1'
#do we have any files of this type?
if len(files[ptype]) > 0:
    #results dictionary
    res[ptype] = {}

    tests = testProperties[ptype]
    #loop over all files
    for f in files[ptype]:
        okay, d = pval.valid_file(f, args={'tomographic_bins': np.arange(50)})
        if okay is False:
            print "Aborting because"
            print "error reading file: " + f
            print "error message: " + d
            sys.exit()

        res[ptype][f] = {}
        pdf = d['PDF']

        #calculate all unweighted metrics for deltaz and deltaz/(1+z)
        for m, tst in enumerate(tests):
            res[ptype][f][testNum] = {}
            for metric in tst['metrics']:
                res[ptype][f][testNum][metric] = get_function(metric)(pdf)

                if key_not_none(tests, 'bins'):
                    binning = tests['bins']
                    res[ptype][f][testNum][metric]['binned_result'] = {}
                    for binDict in binning:
                        if key_not_none(d, binDict) is False:
                            print "You asked to bin in " + binDict
                            print "but it does not exist in file " + f
                            print "aborting"
                            sys.exit()

                        try:
                            bin_vals = eval(binning[binDict])
                        except:
                            print "unable to build the bins, please check syntax: " + binning[binDict]
                            print "Aborting"
                            sys.exit()

                        res[ptype][f]['result'][photoz]['binned_result'][binDict] = {}
                        res[ptype][f]['result'][photoz]['binned_result'][binDict]['bin_column'] = binDict
                        res[ptype][f]['result'][photoz]['binned_result'][binDict]['bin_values'] = bin_vals

                        #this uses the binned_stats function
                        """http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.binned_statistic.html
                        """
                        res[ptype][f]['result'][photoz]['binned_result'][binDict]['bin_stats'] = stats.binned_statistic(d[binDict], diff, bins=bin_vals, statistic=get_function(metrics))

