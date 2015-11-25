import pandas as pd
from astropy.table import Table
import numpy as np
import os
import scipy as sp
from scipy.stats import ks_2samp, binned_statistic
from scipy import interpolate
import sys
from weighted_kde import gaussian_kde
import collections

"""
Authors: Ben Hoyle, Christopher Bonnet

To do:

 define npoisson function

"""

""" ==========================
Cool tools ===================
==============================
"""


def get_function(function_string):
    import importlib
    module, function = function_string.rsplit('.', 1)
    module = importlib.import_module(module)
    function = getattr(module, function)
    return function


#random selects from an (n)d array
def random_choice(arr, axis=None):
    rc = np.shape(arr)
    if len(rc) == 1:
        return arr[np.random.randint(0, rc[0])]
    else:
        indr = np.random.randint(0, rc[1], size=rc[0])
        return arr[np.arange(rc[0]), indr]


#a wrapper around scipy mode
def mode(arr, axis=None):
    return sp.mode(arr, axis=axis)[0]


#extracts the output of a pdf from a decision tree, in this the first column has the highest weight
#this is ordered okay (even for weighted-pdfs) because the first tree has the most weight
def highest_weight(pdf, axis=None):
    return pdf[:, 0]


""" ===========================
Error function checking tools =
===============================
"""


def bootstrap_mean_error(arr, weight, func, Nsamples=None):

    #draw this many samples
    if Nsamples is None:
        Nsamples = 150

    val = np.zeros(Nsamples)
    #what weight do each data have
    p = weight*1.0 / np.sum(weight)

    for i in np.arange(Nsamples):
        #call the function and pass in a bootstrapped sample
        val[i] = func(np.random.choice(arr, size=len(arr), replace=True, p=p))

    #Error is the std of all samples
    return {'mean': np.mean(val), 'sigma': np.std(val)}


def jacknife_error(arr, weight, func, Nsamples=None):
    return False


def bootstrap_mean_error_binned(x, arr, weight, bins, func, Nsamples=None):
    if Nsamples is None:
        Nsamples = 500

    val = np.zeros((Nsamples, len(bins)-1))
    #what weight do each data have
    p = weight*1.0 / np.sum(weight)
    ind = np.arange(len(arr))

    for i in np.arange(Nsamples):
        indrs = np.random.choice(ind, size=len(ind), replace=True, p=p)

        #call the function and pass in a bootstrapped sample
        val[i] = binned_statistic(x[indrs], arr[indrs], bins=bins, statistic=func).statistic

    #Error is the std of all samples
    return {'mean': np.mean(val, axis=0), 'sigma': np.std(val, axis=0)}


def bootstrap_mean_error_pdf_point(_pdf, _bins, _point, _weights, func, Nsamples=None):
    """boot strap error, _pdf is shape (ngals, zbins), _bins = zbins beginings, _weights = galaxy weights
    func = metric function, must accept (pdf[ngals,zbins], zbins, points(ngals)"""

    if Nsamples is None:
        Nsamples = 200

    val = np.zeros(Nsamples)
    #what weight do each data have
    p = _weights * 1.0 / np.sum(_weights)
    ind = np.arange(len(_pdf))

    for i in np.arange(Nsamples):
        indrs = np.random.choice(ind, size=len(ind), replace=True, p=p)

        #call the function and pass in a bootstrapped sample
        val[i] = func(_pdf[indrs], _bins, _point[indrs])

    #Error is the std of all samples
    return {'mean': np.mean(val), 'sigma': np.std(val)}



""" ========================
Data format checking tools =
============================
"""


#check the existence and non-nullness of a key
def key_not_none(_dict, ky):
    if ky in _dict:
        if _dict[ky] is not None:
            return True
    return False


#extract all columns we are asking to work with
def required_cols(_dict, pointOrPdf):
    cols = []
    for dd in _dict:
        if pointOrPdf == 'point':
            if key_not_none(dd, 'point'):
                cols += extract_cols(dd['point'])

        if pointOrPdf == 'pdf':
            if key_not_none(dd, 'pdf'):
                for ttype in ['individual', 'stacks']:
                    if key_not_none(dd['pdf'], ttype):
                        cols += extract_cols(dd['pdf'][ttype])

    return [i for i in np.unique(cols)]


#extract all columns we are asking to work with
def extract_cols(_dict):
    #must have coadd_objects_id
    cols = ['COADD_OBJECTS_ID']
    for i in _dict:
        if key_not_none(_dict, 'predictions'):
            [cols.append(c) for c in _dict['predictions'] if c not in cols]

        if key_not_none(_dict, 'truths'):
            [cols.append(c) for c in [_dict['truths']] if c not in cols]

        if key_not_none(_dict, 'weights'):
            [cols.append(c) for c in [_dict['weights']] if c not in cols]

        if key_not_none(_dict, 'bins'):
            [cols.append(c.keys()[0]) for c in _dict['bins'] if c.keys()[0] not in cols]

        if key_not_none(_dict, 'truth_bins'):
            [cols.append(c.keys()[0]) for c in _dict['truth_bins'] if c.keys()[0] not in cols]

        if key_not_none(_dict, 'metric_bins'):
            [cols.append(c.keys()[0]) for c in _dict['metric_bins'] if c.keys()[0] not in cols]

    return cols


def keytst(_tst):
    for bn in ['bins', 'metric_bins']:
        if key_not_none(_tst, bn):
            for binkyv in _tst[bn]:
                try:
                    bins = eval(binkyv[binkyv.keys()[0]])
                except:
                    print "unable to generate bins in this test"
                    print ' ' + bn + ' ' + binkyv.keys()[0]
                    sys.exit()


def valid_tests(tsts):
    """ check the tests are all valid. """
    for tst in tsts:
        for ttype in tst:
            keytst(ttype)
            for pdftype in ['individual', 'stacks']:
                if key_not_none(ttype, pdftype):
                    keytst(ttype[pdftype])

    return True


def valid_hdf(filename, cols):
    """ Checks that the hdf file is a valid file for our purposes"
    """
    #is args set?
    if cols is None:
        return False, 'you must have at least some defined cols'
    #does the file exist
    if os.path.exists(filename) is False:
        return False, 'file does not exist'

    #can I read the file into a buffer
    try:
        df = pd.read_hdf(filename, 'pdf')
    except:
        return False, 'cannot open file using pandas'

    #is the object a pandas dataframe?
    if type(df) != pd.core.frame.DataFrame:
        return False, 'pandas dataframe not standard'

    #does it have the required columns?
    for i in cols:
        if i not in df:
            return False, 'missing column ' + i

    #does if have the correct number of tomographic bins
    #for i in range(len(args['tomographic_bins'])):
    #    if 'pdf_' + str(i) not in df:
    #        return False, 'missing column ' + 'pdf_' + str(i) + ' of ' + filename

    return True, df


def valid_fits(filename, cols):
    """ checks if the fits file is readable and formatted correctly
    """
    #does file exist
    if os.path.exists(filename) is False:
        return False, 'file does not exist'

    #can we load the file
    try:
        df = Table.read(filename)
    except:
        return False, 'fits file unreadable'

   #is the object a pandas dataframe?
    if type(df) != Table:
        return False, 'astropy table not standard'

    #are all required columns in this file
    for i in cols:
        if i not in df.keys():
            return False, 'missing column ' + i + ' of ' + filename

    return True, df


def valid_file(filename, cols):

    if ('.fit' in filename) or ('.csv' in filename):
        return valid_fits(filename, cols)
    if '.hdf5' in filename:
        return valid_hdf(filename, cols)
    return False, 'currently unable to read file'


""" ================================
Point prediction metrics and tools =
====================================
"""


def delta_z(z_spec, z_phot):
    return z_spec - z_phot


def delta_z_1pz(z_spec, z_phot):
    return delta_z(z_spec, z_phot) / (1 + z_spec)


def sigma_68(arr, axis=None):
    upper, lower = np.percentile(arr, [84.075, 15.825], axis=axis)
    return (upper - lower) / 2.0


def sigma_95(arr, axis=None):
    upper, lower = np.percentile(arr, [97.7, 2.3], axis=axis)
    return (upper - lower) / 2.0


def outlier_rate(arr, outR=None):
    """assumes frac outliers >0.15
    """
    if outR is None:
        outR = 0.15
    return np.sum(np.abs(arr) > outR)*1.0/len(arr)


def outlier_fraction(arr):
    return outlier_rate(arr)*1e2


""" ===================
pdf metrics and tools =
=======================
"""


def normalize_pdf(pdf, z):
    """
    returns normalized pdf
    """
    area = np.trapz(pdf, x=z)
    return pdf / area


def log_loss(act, pred):
    epsilon = 1e-15
    pred = sp.maximum(epsilon, pred)
    pred = sp.minimum(1 - epsilon, pred)
    ll = sum(act * sp.log(pred) + sp.subtract(1, act) * sp.log(sp.subtract(1, pred)))
    ll = ll * -1.0 / len(act)
    return -1 * ll


def kulbachLeiber_bins(arr1, arr2):
    return -1.0 / len(arr1) * np.sum(arr1 * np.log(arr2 + 1e-6))


# what is this test?
def npoisson(dndspecz, dndphotz):
    """ according to https://cdcvs.fnal.gov/redmine/projects/des-photoz/wiki/DC6_Photo-z_Challenge"""
    #unnormalised dn/dz (true number in bins)
    nznorm = (dndphotz - dndspecz) / np.sqrt(dndspecz)
    nzrms = np.sqrt(np.mean(nznorm * nznorm))
    return nzrms


def ks_test(arr1, arr2):
    D, pval = ks_2samp(arr1, arr2)
    return D


def ks_test_prob(arr1, arr2):
    D, pval = ks_2samp(arr1, arr2)
    return pval


def eval_pdf_point(pdf, bins, point):
    val = np.zeros(len(point))
    for i in np.arange(len(pdf)):
        f = interpolate.interp1d(pdf[i], bins)
        val[i] = f(point[i])
    return val


def stackpdfs(pdfs):
    """ numpy shape np.array( (galaxy, pdfbins)) """
    stk_pdf = np.sum(pdfs, axis=0)
    return stk_pdf


def normalisepdfs(pdfs, x):
    """ numpy shape (galaxy, pdfbins) """

    if len(np.shape(pdfs)) > 1:
        smm = np.trapz(pdfs, x, axis=1)
    else:
        smm = np.trapz(pdfs, x)

    #tricks to allow array (ngal, pdfbins) to be divided by array (ngal)
    npdfs = (pdfs.T / smm).T

    return npdfs


def integrate_dist_bin(dfs, x, minval, maxval):
    """Please note! we assume that the minval and maxval are (or are close to) x bin edges. If this is not the case
    This routine will be approximate
    """
    ind = (x >= minval) * (x <= maxval)

    if len(np.shape(dfs)) > 1:
        smm = np.trapz(dfs[:, ind], x[ind], axis=1)
    else:
        smm = np.trapz(dfs[ind], x[ind])

    return smm


def cumaltive_to_point(dfs, x, points):
    """ This assumes bin widths are small enough to ignore.
    Expected shape: numpy dfs shape (galaxy, bins), x begins at 0, ends at just before 2 """
    """Note: all points < x[0] are set to x[0]"""
    """Note: all points > x[-1] are set to x[-1]"""

    if isinstance(points, collections.Iterable):
        points[points > x[-1]] = x[-1]
        points[points < x[0]] = x[0]
    else:
        if points > x[-1]:
            points = x[-1]
        if points < x[0]:
            points = x[0]

    #calcaulte the cumaltive distribution
    if len(np.shape(dfs)) > 1:
        cum = np.cumsum(dfs, axis=1)
        cum = (cum.T / cum[:, -1]).T

        #determine which index is just <= each point value]
        indx = np.arange(len(x))
        ind_max = np.array([np.amax(indx[x <= p]) for p in points])

        #return these cfd values at the points of interst
        return cum[np.arange(len(dfs)), ind_max]
    else:
        cum = np.cumsum(dfs)
        cum = cum / float(cum[-1])
        return cum[x <= points][-1]


def xval_cumaltive_at_ypoint(dfs, x, point):
    """returns the xvals of dfs(ngal, nxbins), x(bins in xdir) point(1), for a point
    which sits on the y-axis of the cdf"""

    """Note: all points < x[0] are set to x[0]"""
    """Note: all points > x[-1] are set to x[-1]"""

    if point > x[-1]:
        point = x[-1]
    if point < x[0]:
        point = x[0]

    if len(np.shape(dfs)) > 1:
        cum = np.cumsum(dfs, axis=1)
        cum = (cum.T / cum[:, -1]).T

        #determine which index is just <= each point value]
        xarr = np.array([np.amax(x[c <= point]) for c in cum])

        #return these x values  correspondong to the y-axis = point interst cdf
        return xarr
    else:
        cum = np.cumsum(dfs)
        cum = cum / float(cum[-1])
        return x[cum <= point][-1]


def gini(sorted_list):
    """ from https://planspacedotorg.wordpress.com/2013/06/21/how-to-calculate-gini-coefficient-from-raw-data-in-python/
    """
    cum = np.cumsum(sorted_list)
    cum = cum * 1.0 / cum[-1]
    cumx = np.cumsum([1] * len(cum))
    cumx = cumx * 1.0 / cumx[-1]
    gini_v = np.mean(np.abs(cumx - cum)) * 2.0
    return gini_v


def Bordoloi_pdf_test(dfs, x, points, n_yaxis_bins=None):
    """determine a test described here:
    http://arxiv.org/pdf/1201.0995.pdf
    Is the pdf flat, in bins of Pval?
    """

    #calculate the cum pdf values up to each point (or spec-z)
    pvals = cumaltive_to_point(dfs, x, points)

    if n_yaxis_bins is None:
        n_yaxis_bins = 100

    #make a histogram of this in steps of 0.01. Must be smaller than 0.03 otherwise
    hp = np.histogram(pvals, bins=np.arange(n_yaxis_bins + 1) / float(n_yaxis_bins))[0]

    #return the giniCoefficient of this distribution. Flat = 0 = great!
    return gini(hp)


def dfs_mode(dfs, x):
    """calcaulte the location of the mode of a heap of dfs"""
    if len(np.shape(dfs)) > 1:
        mx = np.argmax(dfs, axis=1)
        return x[mx]
    else:
        return x[np.argmax(dfs)]


def binned_pdf_point_stats(_data_to_bin, _bin_vals, _pdf, _zbins, _truths, _weights, func, Nsamples=None):
    """ take a binning vector, and bins edges, and calculate the pdf-point statistic. Use Boot strap resampling to use weights keyword. If weights are given, we use bootstrap to estimate the weighted average value of the metric

    _data_to_bin[size ngals], _bin_vals = bins to bin _data_to_bin, _pdf[ngals, nzbins], _zbins the zbins of the pdf
    _truths[ngals] point predictions, _weights[ngals] galaxy weights, func = metric of choice"""

    if Nsamples is None:
        Nsamples = 150
    res = {}
    p = _weights / np.sum(_weights)
    for i in np.arange(len(_bin_vals)-1):
        res[i] = {}

        if i == len(_bin_vals)-1:
            #catch edge case
            ind = (_data_to_bin >= _bin_vals[i]) * (_data_to_bin <= _bin_vals[i+1])
        else:
            ind = (_data_to_bin >= _bin_vals[i]) * (_data_to_bin < _bin_vals[i+1])

        res[i]['weighted_bin_center'] = np.average(_data_to_bin[ind], weights=p[ind])

        #do we need to both reweighting?
        if np.sum(_weights[ind]) != np.sum(ind):

            val = np.zeros(Nsamples)
            indx = np.arange(len(_data_to_bin))[ind]
            for j in np.arange(Nsamples):
                indr = np.random.choice(indx, size=len(indx), replace=True, p=p[ind])
                val[i] = func(_pdf[indr], _zbins, _truths[indr])
            res[i]['weighted_value'] = np.mean(val)
        else:
            res[i]['weighted_value'] = func(_pdf[ind], _zbins, _truths[ind])

    return res


def interpolate_dist(_df1, _bins1, _bins2, kind=None):
    """interpolate df1 at locations _bins2 """
    if kind is None:
        kind = 'linear'
    I = interpolate.interp1d(_bins1, _df1, kind=kind)

    return I(_bins2)


def binned_statistic_dist1_dist2(arr_, bin_vals_, truths_, truth_bins_, pdf_, pdf_z_center_, funv, weights=None):
    """ bins the stacked pdf and truths in bins of arr, then calculates the metric on these distributions"""
    """ metric doesn't work so well with weights at the mo"""

    if weights is None:
        weights = np.ones(len(arr_))

    res = {}
    for i in np.arange(len(bin_vals_)-1):
        ind_ = (arr_ >= bin_vals_[i]) * (arr_ < bin_vals_[i + 1])
        if np.sum(ind_) > 0:
            res[i] = {}
            truth_dist_ = np.histogram(truths_[ind_], bins=truth_bins_)[0]
            truth_bins_centers_ = stats.binned_statistic(truths[ind_], truths[ind_], bins=truth_bins_, statistic=np.mean).statistic
            stacked_pdf_ = pval.stackpdfs(pdf_[ind_])
            stckd_pdfs_at_trth_cntrs_ = pval.interpolate_dist(stacked_pdf_, pdf_z_center_, truth_bins_centers_)

            res[i]['weighted_bin_center'] = np.average(_data_to_bin[ind], weights=p[ind])

            """ Add weights in here  """
            res[i]['weighted_value'] = func(truth_dist, stckd_pdfs_at_trth_cntrs)

    return res

""" ==========================
validation metrics and tools =
==============================
"""


def ld_writedicts(filepath, dictionary):
    f = open(filepath, 'w')
    newdata = dumps(dictionary, 1)
    f.write(newdata)
    f.close()


def ld_readdicts(filepath):
    f = open(filepath, 'r')
    d = load(f)
    f.close()
    return d


def _bias(z_spec, z_phot, weight=None):
    dz1 = (z_spec - z_phot)
    if weight is None:
        bias = np.mean(dz1)
    else:
        bias = np.average(dz1, weights=weight)

    return bias


def _normalize_pdf(pdf, dz):
    """
    returns normalized pdf
    """
    area = np.trapz(pdf, dx=dz)
    return pdf / area


def _sigma(z_spec, z_phot, weight=None):
    dz1 = (z_spec - z_phot)
    if weight is None:
        sigma = np.std(dz1)
    else:
        sigma = _w_std(dz1, weights=weight)
    return sigma


def _percentile(n, percent):
    n = np.sort(n)
    k = (len(n) - 1) * percent
    f = np.floor(k)
    c = np.ceil(k)
    if f == c:
        return n[int(k)]
    d0 = n[int(f)] * (c - k)
    d1 = n[int(c)] * (k - f)
    return d0 + d1



def _w_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((np.abs(values - average)) ** 2, weights=weights)  # Fast and numerically precise
    return np.sqrt(variance)





def mean(df, binning, z_phot, metric='mean', weights=None, tomo_bins=np.array([0, 5.0]), n_resample=50):
    """
    :param df: pandas data-frame
    :param binning: center of redshift bins
    :param metric : 'mean','mode' of the pdf as the point estimate
    :param weights: optional weighting scheme with same len as df
    :param tomo_bins: in which z-bins array exp [0.0, 0.2, 0.6, 1.8]
    :param n_resample : Amount of resamples to estimate mean and variance on the mean.
    :return: pandas data frame with mean estimates
    """

    assert isinstance(df, pd.DataFrame), 'df must be a pandas DataFrame'
    assert isinstance(binning, np.ndarray), 'binning must be a numpy array'
    if weights:
        assert weights in df.columns, str(weights) + ' not in df.columns'
        df[weights] = (df[weights] / df[weights].sum()).values  # normalize weights
    else:
        df[weights] = 1.0 / len(df)  # set uniform weights if none given
    assert isinstance(z_phot, np.ndarray), 'z_phot must be a numpy array'
    assert len(z_phot) == len(df), 'Length of z_phot must be equal to that of df'
    df['phot_sel'] = z_phot  # Make the selection photo-z a part of the DataFrame
    assert 'Z_SPEC' in df.columns, 'The df needs a "z_spec" columns'
    pdf_names = [c for c in df.keys() if 'pdf_' in str(c)]

    if metric == 'mode':
        df['z_phot'] = binning[np.argmax(df[pdf_names].values, axis=1)]
    elif metric == 'mean':
        df['z_phot'] = np.inner(binning, df[pdf_names].values)
    else:
        print 'Metric needs to be either "mode" or "mean", using mean !'
        df['z_phot'] = np.inner(binning, df[pdf_names].values)

    mean_spec_bin_array = []
    err_mean_phot_bin_array = []

    mean_phot_bin_array = []
    err_mean_spec_bin_array = []

    w_mean_spec_bin_array = []
    w_err_mean_spec_bin_array = []

    w_mean_phot_bin_array = []
    w_err_mean_phot_bin_array = []

    df_index = []
    df_count = []
    mean_z_bin = []

    print len(tomo_bins)

    for j in range(len(tomo_bins) - 1):
        sel = (df.phot_sel > tomo_bins[j]) & (df.phot_sel <= tomo_bins[j + 1])
        if sel.sum() > 0:
            df_index.append('z [' + str(tomo_bins[j]) + ', ' + str(tomo_bins[j + 1]) + ']')
            df_count.append(sel.sum())
            mean_z_bin.append((tomo_bins[j] + tomo_bins[j + 1]) / 2.0)
            df_sel = df[sel]

            mean_spec_array = []
            mean_phot_array = []

            w_mean_spec_array = []
            w_mean_phot_array = []

            for i in xrange(n_resample):
                df_sample = df_sel.sample(n=len(df_sel), replace=True, weights=None)
                mean_spec_array.append(df_sample.Z_SPEC.mean())
                mean_phot_array.append(df_sample.z_phot.mean())

                df_sample = df_sel.sample(n=len(df_sel), replace=True, weights=df_sel[weights])
                w_mean_spec_array.append(np.average(df_sample.Z_SPEC))
                w_mean_phot_array.append(np.average(df_sample.z_phot))

            mean_spec = np.mean(mean_spec_array)
            mean_phot = np.mean(mean_phot_array)

            err_mean_spec = np.std(mean_spec_array)
            err_mean_phot = np.std(mean_phot_array)

            w_mean_spec = np.mean(w_mean_spec_array)
            w_mean_phot = np.mean(w_mean_phot_array)

            w_err_mean_spec = np.std(w_mean_spec_array)
            w_err_mean_phot = np.std(w_mean_phot_array)

            mean_spec_bin_array.append(mean_spec)
            mean_phot_bin_array.append(mean_phot)

            err_mean_spec_bin_array.append(err_mean_spec)
            err_mean_phot_bin_array.append(err_mean_phot)

            w_mean_spec_bin_array.append(w_mean_spec)
            w_mean_phot_bin_array.append(w_mean_phot)

            w_err_mean_spec_bin_array.append(w_err_mean_spec)
            w_err_mean_phot_bin_array.append(w_err_mean_phot)

        else:
            print 'Bin ' + str(j) + ' has no objects and it not included in the summary (counting starts at zero)'

    to_pandas = np.vstack((mean_z_bin, df_count,
                           mean_spec_bin_array, err_mean_spec_bin_array,
                           mean_phot_bin_array, err_mean_phot_bin_array,
                           w_mean_spec_bin_array, w_err_mean_spec_bin_array,
                           w_mean_phot_bin_array, w_err_mean_phot_bin_array,
                           )).T

    return_df = pd.DataFrame(to_pandas, columns=['mean_z_bin', 'n_obj',
                                                 'mean_spec', 'err_mean_spec',
                                                 'mean_phot', 'err_mean_phot',
                                                 'w_mean_spec', 'w_err_mean_spec',
                                                 'w_mean_phot', 'w_err_mean_phot',
                                                 ])
    return_df.index = df_index

    return return_df


def weighted_nz_distributions(df, binning, weights=None, tomo_bins=np.array([0, 5.0]), z_phot=None, n_resample=50):
    """
    :param df: pandas data-frame
    :param binning: center of redshift bins
    :param weights: optional weighting scheme with same len as df
    :param tomo_bins: in which z-bins array exp [0.0, 0.2, 0.6, 1.8]
    :param n_resample : Amount of resamples to estimate mean and variance on the mean.
    :return: dictionaries with estimates of weighted n(z) and bootstrap estimates
    """

    assert isinstance(df, pd.DataFrame), 'df must be a pandas DataFrame'
    assert isinstance(binning, np.ndarray), 'binning must be a numpy array'
    if weights:
        assert weights in df.columns, str(weights) + ' not in df.columns'
        df[weights] = (df[weights] / df[weights].sum()).values  # normalize weights
    else:
        df[weights] = 1.0 / float(len(df))  # set uniform weights if none given

    assert isinstance(z_phot, np.ndarray), 'z_phot must be a numpy array'
    assert len(z_phot) == len(df), 'Length of z_phot must be equal to that of df'
    df['phot_sel'] = z_phot  # Make the selection photo-z a part of the DataFrame
    assert 'Z_SPEC' in df.columns, 'The df needs a "Z_SPEC" in df.columns'
    pdf_names = [c for c in df.keys() if 'pdf_' in str(c)]

    phot_iter = {}
    spec_iter = {}

    # In the following section the tomographic bins are treated

    for j in xrange(0, len(tomo_bins) - 1):
        sel = (df.phot_sel > tomo_bins[j]) & (df.phot_sel <= tomo_bins[j + 1])
        if sel.sum() > 0:
            df_sel = df[sel]

            phot_iter[j + 1] = {}
            spec_iter[j + 1] = {}

            for i in xrange(n_resample):
                df_sample = df_sel.sample(n=len(df_sel), replace=True, weights=df_sel[weights])
                kde_w_spec_pdf = gaussian_kde(df_sample['Z_SPEC'].values, bw_method='silverman')
                kde_w_spec_pdf = kde_w_spec_pdf(binning)

                phot_iter[j + 1][i + 1] = _normalize_pdf(df_sample[pdf_names].sum(), binning[1] - binning[0]).values
                spec_iter[j + 1][i + 1] = kde_w_spec_pdf

            phot_iter[j + 1][0] = _normalize_pdf(df_sel[pdf_names].sum(), binning[1] - binning[0]).values
            kde_w_spec_pdf = gaussian_kde(df_sel['Z_SPEC'].values, bw_method='silverman')
            spec_iter[j + 1][0] = kde_w_spec_pdf(binning)

    # In the following section the full n(z) is treated i.e not in tomographic bins

    sel = (df.phot_sel > tomo_bins[0]) & (df.phot_sel <= tomo_bins[len(tomo_bins) - 1])
    df_sel = df[sel]
    phot_iter[0] = {}
    spec_iter[0] = {}

    for i in xrange(n_resample):
        df_sample = df_sel.sample(n=len(df_sel), replace=True, weights=df_sel[weights])
        kde_w_spec_pdf = gaussian_kde(df_sample['Z_SPEC'].values, bw_method='silverman')
        kde_w_spec_pdf = kde_w_spec_pdf(binning)

        phot_iter[0][i + 1] = _normalize_pdf(df_sample[pdf_names].sum(), binning[1] - binning[0]).values
        spec_iter[0][i + 1] = kde_w_spec_pdf

    phot_iter[0][0] = _normalize_pdf(df_sel[pdf_names].sum(), binning[1] - binning[0]).values
    kde_w_spec_pdf = gaussian_kde(df_sel['Z_SPEC'].values, bw_method='silverman')
    spec_iter[0][0] = kde_w_spec_pdf(binning)

    data_for_wl = {'binning': binning, 'phot': phot_iter, 'spec': spec_iter}

    return phot_iter, spec_iter, data_for_wl