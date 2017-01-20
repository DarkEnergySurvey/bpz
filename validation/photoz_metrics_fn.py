#! /usr/bin/env python

"""photo-z tests to perform.
call like:
import photoz_metrics_fn as pzt
import yaml

tst_yaml = yaml.load(open('PathToTestConfig/photoz.yaml', 'r'))
tst_yaml = tst_yaml['point']

d = Table.load('PhotoZPredictionsFile.fits')

res = pzt.perform_tests(d, tst_yaml)

"""
import numpy as np
numpy = np
import bh_photo_z_validation as pval
import vlfn
from scipy import stats
import cPickle as pickle


def get_function(function_string):
    import importlib
    module, function = function_string.rsplit('.', 1)
    module = importlib.import_module(module)
    function = getattr(module, function)
    return function

#get the galaxy weights
def get_weights(_dict, _ky, _d):
    #set all objects equal weight, unless defined
    if pval.key_not_none(_dict, _ky) is False:
        #print ("you have not set any weights for this test")
        #print ("continuing with weights=1")
        weights = np.ones(len(_d))
    else:
        weights = _d[tst['weights']]

    return weights / np.sum(weights)

"""
-- this function is broken because of change to .yaml file.
def perform_tests(d, test_dict):

    #results dictionary

    res = {}

    #should we calculate an error on these metrics
    error_function = pval.key_not_none(test_dict, 'error_function')

    err_metric = {}
    if error_function:
        for ef in test_dict['error_function']:
            #turn error function.string into a function
            err_metric[ef.split('.')[-1]] = pval.get_function(ef)

    for photoz in test_dict['predictions']:
        res[photoz] = {}
        res[photoz]['metrics_z1_z2'] = {}
        res[photoz]['metrics_diffz'] = {}

        z_truth = np.array(d[test_dict['truths']])
        z_pred = np.array(d[photoz])
        
        for metric in test_dict['metrics_z1_z2']:

            res[photoz]['metrics_z1_z2'][metric] = {}

            weights = get_weights(test_dict, 'weights', d)

            #turn string into function
            metric_function = pval.get_function(metric)

            res[photoz]['metrics_z1_z2'][metric] = {}
            res[photoz]['metrics_z1_z2'][metric]['VALUE'] = np.asscalar(metric_function(z_truth, z_pred, weights=weights))


            #shall we calculate binning statiscs?
            if pval.key_not_none(test_dict, 'bins'):
                binning = test_dict['bins']

                res[photoz]['metrics_z1_z2'][metric]['bins'] = {}
                for binDict in binning:
                    ky = binDict.keys()[0]
                    bin_vals = eval(binDict[ky])

                    res[photoz]['metrics_z1_z2'][metric]['bins'][ky] = {}

                    bn_stat = np.zeros(len(bin_vals)-1) -1 
                    bn_cntr_sts = np.zeros(len(bin_vals)-1) -1
                    for bbn in range(len(bin_vals)-1):
                        ind_bn = (d[ky] <= bin_vals[bbn + 1]) * (d[ky] > bin_vals[bbn])
                        if np.sum(ind_bn) > 1:
                            bn_cntr_sts[bbn] = np.mean(d[ky][ind_bn])
                            bn_stat[bbn] = metric_function(z_truth[ind_bn], z_pred[ind_bn], weights=weights[ind_bn])

                    res[photoz]['metrics_z1_z2'][metric]['bins'][ky]['BIN_CENTERS'] = [np.asscalar(vv) for vv in bn_cntr_sts]
                    res[photoz]['metrics_z1_z2'][metric]['bins'][ky]['VALUE'] = [np.asscalar(vv) for vv in bn_stat]

        #calculate stats on diff=z1-z2 and diff_1pz=(z1-z2)/(1+z1)
        diff = pval.delta_z(d[test_dict['truths']], d[photoz])
        diff_1pz = pval.delta_z_1pz(d[test_dict['truths']], d[photoz])

        points = {'delta_z': diff, 'diff_1pz': diff_1pz}

        for metric in test_dict['metrics_diffz']:

            res[photoz]['metrics_diffz'][metric] = {}

            #set all objects equal weight, unless defined
            weights = get_weights(test_dict, 'weights', d)


            #turn string into function
            metric_function = pval.get_function(metric)

            #which residuals shall we employ?
            for diffpp in points.keys():
                res[photoz]['metrics_diffz'][metric][diffpp] = {}
                res[photoz]['metrics_diffz'][metric][diffpp]['VALUE'] = np.asscalar(metric_function(points[diffpp]))

                #calculate errors on these metrics
                for ef in err_metric:
                    bstamp_mean_err = err_metric[ef](points[diffpp], weights, metric_function)
                    res[photoz]['metrics_diffz'][metric][diffpp]['MEAN_' + ef] = np.asscalar(bstamp_mean_err['mean'])
                    res[photoz]['metrics_diffz'][metric][diffpp]['SIGMA_' + ef] = np.asscalar(bstamp_mean_err['sigma'])

                #shall we calculate binning statiscs?
                if pval.key_not_none(test_dict, 'bins'):
                    binning = test_dict['bins']

                    res[photoz]['metrics_diffz'][metric][diffpp]['bins'] = {}
                    for binDict in binning:
                        ky = binDict.keys()[0]
                        bin_vals = eval(binDict[ky])

                        res[photoz]['metrics_diffz'][metric][diffpp]['bins'][ky] = {}
                        #this uses the binned_stats function
                        #http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.stats.binned_statistic.html
                        

                        #calculate the unweighted statistics in each bin
                        bn_stats = stats.binned_statistic(d[ky], points[diffpp], bins=bin_vals, statistic=metric_function)

                        #determine the center of each bin
                        bn_cntr_sts = stats.binned_statistic(d[ky], d[ky], bins=bin_vals, statistic=np.mean)

                        res[photoz]['metrics_diffz'][metric][diffpp]['bins'][ky]['BIN_CENTERS'] = [np.asscalar(vv) for vv in bn_cntr_sts.statistic]
                        res[photoz]['metrics_diffz'][metric][diffpp]['bins'][ky]['VALUE'] = [np.asscalar(vv) for vv in bn_stats.statistic]

                        #calculate the mean and error by bootstrapping
                        bn_bs_stats = pval.bootstrap_mean_error_binned(d[ky], points[diffpp], weights, bin_vals, metric_function)

                        #calculate the bin 'centers' by boot strapping
                        bn_bs_cnters = pval.bootstrap_mean_error_binned(d[ky], d[ky], weights, bin_vals, np.mean)

                        res[photoz]['metrics_diffz'][metric][diffpp]['bins'][ky]['BIN_CENTERS_MEAN_BS'] = [np.asscalar(vv) for vv in bn_bs_cnters['mean']]

                        res[photoz]['metrics_diffz'][metric][diffpp]['bins'][ky]['BIN_CENTERS_SIGMA_BS'] = [np.asscalar(vv) for vv in bn_bs_cnters['sigma']]

                        res[photoz]['metrics_diffz'][metric][diffpp]['bins'][ky]['MEAN_BS'] = [np.asscalar(vv) for vv in bn_bs_stats['mean']]
                        res[photoz]['metrics_diffz'][metric][diffpp]['bins'][ky]['SIGMA_BS'] = [np.asscalar(vv) for vv in bn_bs_stats['sigma']]
    return res
"""

def perform_tests_fast(d, tst):
    """perform_tests_fast performs a fast set of tests without boot strap resampling, or any error resampling """
    #results dictionary
    res = {}

     #get all the columns we are gonna test on
    reqcols = tst['metrics'].keys()

    for i in reqcols:
        if i not in d.keys():
            print "missing column ", i
            sys.exit()

    #which redshift do we need for this metric test?
    for photoz in tst['metrics']:
        res[photoz] = {}

        z_truth = np.array(d[tst['truths']])
        z_pred = np.array(d[photoz])

        #what is the metric test?
        for metric in tst['metrics'][photoz]:

            res[photoz][metric] = {}

            #convert metric name to metric function
            metric_function = pval.get_function(metric)

            #do I have to pass any additional arguments to this function?
            extra_params = pval.get_extra_params(tst, metric)

            #what weighting scheme shall we apply?
            for wght in tst['weights']:
                res[photoz][metric][wght] = {}

                #get the data weights
                weights = np.array(d[wght], dtype=float)

                res[photoz][metric][wght]['value'] = vlfn.process_function(metric_function, z_truth,
                    z_pred, weights=weights, extra_params=extra_params)


                #shall we calculate binning statiscs?
                if pval.key_not_none(tst, 'bins'):
                    binning = tst['bins']

                res[photoz][metric][wght]['bins'] = {}
                for ky in binning:
                    bin_vals = binning[ky]

                    res[photoz][metric][wght]['bins'][ky] = {}
                    res[photoz][metric][wght]['bins'][ky]['bin_center'] = []
                    res[photoz][metric][wght]['bins'][ky]['value'] = []

                    for bbn in range(len(bin_vals)-1):
                        ind_bn = (d[ky] <= bin_vals[bbn + 1]) * (d[ky] > bin_vals[bbn])
                        if np.sum(ind_bn) > 1 and np.sum(weights[ind_bn]) > 0:
                            res[photoz][metric][wght]['bins'][ky]['bin_center'].append(np.mean(d[ky][ind_bn]))
                            res[photoz][metric][wght]['bins'][ky]['value'].append(vlfn.process_function(metric_function, z_truth[ind_bn], z_pred[ind_bn], weights=weights[ind_bn], extra_params=extra_params))

    return res