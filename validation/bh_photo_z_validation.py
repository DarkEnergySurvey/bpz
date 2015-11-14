import pandas as pd
from astropy.table import Table
import numpy as np
import os
import scipy as sp
from scipy.stats import ks_2samp


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


#a wrapper around numpy mode
def mode(arr, axis=None):
    return np.mode(arr, axis=axis)[0]


def highest_weight(pdf, axis=None):
    return pdf[:, 0]

""" ===========================
Error function checking tools =
===============================
"""


def bootstrap_mean_error(arr, weight, func):

    #draw this many samples
    Nsamples = 500
    val = np.zeros(Nsamples)
    #what weight do each data have
    p = weight*1.0 / np.sum(weight)

    for i in np.arange(Nsamples):
        #call the function and pass in a bootstrapped sample
        val[i] = func(np.random.choice(arr, size=len(arr), replace=True, p=p))

    #Error is the std of all samples
    return {'mean': np.mean(val), 'sigma': np.std(val)}


def jacknife_error(arr, weight, func):
    return False


""" ========================
Data format checking tools =
============================
"""


def valid_hdf(filename, args=None):
    """ Checks that the hdf file is a valid file for our purposes"
    """
    #is args set?
    if args is None:
        return False, 'you must have at least args={tomographic_bins:}'
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
    for i in ['COADD_OBJECTS_ID', 'Z_SPEC', 'WEIGHT']:
        if i not in df:
            return False, 'missing column ' + i

    #does if have the correct number of tomographic bins
    for i in range(len(args['tomographic_bins'])):
        if 'pdf_' + str(i) not in df:
            return False, 'missing column ' + 'pdf_' + str(i) + ' of ' + filename

    return True, df


def valid_fits(filename):
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
    for i in ['MODE_Z', 'MEAN_Z', 'Z_MC', 'COADD_OBJECTS_ID', 'Z_SPEC', 'MAG_DETMODEL_I']:
        if i not in df.keys():
            return False, 'missing column ' + i + ' of ' + filename

    return True, df


def valid_file(filename, args=None):
    if '.fit' in filename:
        return valid_fits(filename)
    if '.hdf5' in filename:
        return valid_hdf(filename, args)
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


def outlier_rate(arr):
    """assumes frac outliers >0.15
    """
    return np.sum(np.abs(arr) > 0.15)*1.0/len(arr)


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
def npoisson(arr1, arr2):
    return 0


def ks_test(arr1, arr2):
    D, pval = ks_2samp(arr1, arr2)
    return D


def ks_test_prob(arr1, arr2):
    D, pval = ks_2samp(arr1, arr2)
    return pval


""" ==========================
validation metrics and tools =
==============================
"""


# should tolerances be on absolute values?
def within_tolerance(val1, val2, tol):
    return np.abs(val1 - val2) < tol






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
    assert 'z_spec' in df.columns, 'The df needs a "z_spec" columns'
    assert 'pdf_0' in df.columns, 'The pdf values must have pdf_i column name'
    pdf_names = ['pdf_' + str(i) for i in range(500) if 'pdf_' + str(i) in df.columns]

    if metric == 'mode':
        df['z_phot'] = binning[np.argmax([df[pdf_names].values], axis=1)][0]
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
                mean_spec_array.append(df_sample.z_spec.mean())
                mean_phot_array.append(df_sample.z_phot.mean())

                df_sample = df_sel.sample(n=len(df_sel), replace=True, weights=df_sel[weights])
                w_mean_spec_array.append(np.average(df_sample.z_spec))
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
    assert 'z_spec' in df.columns, 'The df needs a "z_spec" in df.columns'
    pdf_names = ['pdf_' + str(i) for i in range(500) if 'pdf_' + str(i) in df.columns]

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
                kde_w_spec_pdf = gaussian_kde(df_sample.z_spec.values, bw_method='silverman')
                kde_w_spec_pdf = kde_w_spec_pdf(binning)

                phot_iter[j + 1][i + 1] = _normalize_pdf(df_sample[pdf_names].sum(), binning[1] - binning[0]).values
                spec_iter[j + 1][i + 1] = kde_w_spec_pdf

            phot_iter[j + 1][0] = _normalize_pdf(df_sel[pdf_names].sum(), binning[1] - binning[0]).values
            kde_w_spec_pdf = gaussian_kde(df_sel.z_spec.values, bw_method='silverman')
            spec_iter[j + 1][0] = kde_w_spec_pdf(binning)

    # In the following section the full n(z) is treated i.e not in tomographic bins

    sel = (df.phot_sel > tomo_bins[0]) & (df.phot_sel <= tomo_bins[len(tomo_bins) - 1])
    df_sel = df[sel]
    phot_iter[0] = {}
    spec_iter[0] = {}

    for i in xrange(n_resample):
        df_sample = df_sel.sample(n=len(df_sel), replace=True, weights=df_sel[weights])
        kde_w_spec_pdf = gaussian_kde(df_sample.z_spec.values, bw_method='silverman')
        kde_w_spec_pdf = kde_w_spec_pdf(binning)

        phot_iter[0][i + 1] = _normalize_pdf(df_sample[pdf_names].sum(), binning[1] - binning[0]).values
        spec_iter[0][i + 1] = kde_w_spec_pdf

    phot_iter[0][0] = _normalize_pdf(df_sel[pdf_names].sum(), binning[1] - binning[0]).values
    kde_w_spec_pdf = gaussian_kde(df_sel.z_spec.values, bw_method='silverman')
    spec_iter[0][0] = kde_w_spec_pdf(binning)

    data_for_wl = {'binning': binning, 'phot': phot_iter, 'spec': spec_iter}

    return phot_iter, spec_iter, data_for_wl