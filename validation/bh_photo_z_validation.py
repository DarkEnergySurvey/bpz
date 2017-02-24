import pandas as pd
from astropy.table import Table
import numpy as np

#don't remove double numpy please :D
import numpy as numpy
import os
import scipy as sp
from scipy.stats import ks_2samp, binned_statistic
from scipy import interpolate
from scipy.stats import gaussian_kde
from weighted_kde import gaussian_kde as gss_kde
import sys
#from weighted_kde import gaussian_kde
import collections
from cPickle import dumps, load

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM


"""
Authors: Ben Hoyle, Christopher Bonnet

To do:
 Enable mad() functions to accept axis= keyword for 
 Nd arrays
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
    return sp.stats.mode(arr, axis=axis)[0]

def get_extra_params(tst_, metric_name):
    """check the config test file. If the extra_params keywoard is set
    and it has a key=== metric_name, then return the extra params
    else return None"""

    if tst_['extra_params'] is None:
        return None
    extra_params = key_not_none(tst_, 'extra_params')
    if extra_params is not None:
        if key_not_none(tst_['extra_params'], metric_name):
            extra_params = tst_['extra_params'][metric_name]
        else:
            extra_params = None
    return extra_params


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
    prob = np.array(weight, dtype=float)
    prob *= 1.0 / np.sum(prob, dtype=float)

    for i in np.arange(Nsamples):
        #call the function and pass in a bootstrapped sample
        val[i] = func(np.random.choice(arr, size=len(arr), replace=True, p=prob))

    #Error is the std of all samples
    return {'mean': np.mean(val), 'sigma': np.std(val)}


def jacknife_error(arr, weight, func, Nsamples=None):
    return False


def bootstrap_mean_error_binned(x, arr, weight, bins, func, Nsamples=None):
    if Nsamples is None:
        Nsamples = 500

    val = np.zeros((Nsamples, len(bins)-1))
    #what weight do each data have
    p = weight * 1.0 / np.sum(weight)
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
    prob = _weights * 1.0 / np.sum(_weights, dtype=float)
    prob = prob / np.sum(prob)

    ind = np.arange(len(_pdf))

    for i in np.arange(Nsamples):
        indrs = np.random.choice(ind, size=len(ind), replace=True, p=prob)

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
        cols += extract_cols(dd)

        for ttype in ['individual', 'stacks']:
            if key_not_none(dd, ttype):
                cols += extract_cols(dd[ttype])

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
    """This function attemps to turn a bin string into executed code
    e.g. 'bins': {'Z_SPEC': 'numpy.arange(5)'}
    should turn 'numpy.arange(5)' into numpy.arange(5) array
    It returns False if it cannot compile the string code
    """
    for bn in ['bins', 'truth_bins']:
        if key_not_none(_tst, bn):
            for binkyv in _tst[bn]:
                try:
                    bins = eval(binkyv[binkyv.keys()[0]])
                except:
                    print "unable to generate bins in this test"
                    print ' ' + bn + ' ' + binkyv.keys()[0]
                    print binkyv[binkyv.keys()[0]]
                    sys.exit()

def valid_tests(tsts):
    """ check the tests are all valid. """
    for tst in tsts:
        keytst(tst)
        for pdftype in ['individual', 'stacks']:
            if key_not_none(tst, pdftype):
                keytst(tst[pdftype])

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

    if ('.fit' in filename[-7:]) or ('.csv' in filename[-7:]):
        return valid_fits(filename, cols)
    if '.hdf5' in filename[-5:]:
        return valid_hdf(filename, cols)
    return False, 'currently unable to read file'


""" ================================
Point prediction metrics and tools =
====================================
"""


def wl_metric(z1, z2):
    """Determine the WL metric of choice
    |<z1> - <z2>|
    """
    return np.abs(np.mean(z1) - np.mean(z2))


def delta_sigma_crit(z1, z2, z2weight, z_lens):
    """Determine the int{P_w(photz) * Dds/Ds(z)  dz }
    Where P_w(photz) = sum(z_mc_i * weight_i), and Dds, Ds are angular diam distances
    p_w_phot = sum of weight in each
    z1 = true-zdistribution [N-element array]
    z2 = Z_MC from photo-z code [N-element array]
    z2weight = weight associated to each z2 object [N-element array]
    z_lens = the z_lens redshift [float]
    """
    binCenters = np.linspace(0, 3, 300)
    p_w_phot = np.zeros_like(binCenters)
    p_w_true = np.zeros_like(binCenters)
    for i in range(len(binCenters)-1):
        ind = (z2 > binCenters[i]) * (z2 < binCenters[i+1])
        p_w_phot[i] = np.sum(z2weight[ind])*1.0
        ind = (z1 > binCenters[i]) * (z1 < binCenters[i+1])
        p_w_true[i] = np.sum(z2weight[ind])*1.0

    p_w_phot[p_w_phot < 0] = 0
    
    int_ = np.trapz(p_w_phot, x=binCenters)

    p_w_phot = p_w_phot/int_

    int_ = np.trapz(p_w_true, x=binCenters)

    p_w_true = p_w_true/int_

    #p_w_true = gaussian_kde(z1).evaluate(binCenters)

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    Ds = cosmo.angular_diameter_distance(binCenters)
    Dds = cosmo.angular_diameter_distance_z1z2(z_lens, binCenters)
    DCs_Ds = Dds / Ds

    DCs_Ds[binCenters < z_lens] = 0
    p_w_phot = p_w_phot / np.trapz(p_w_phot, x=binCenters)
    p_w_true = p_w_true / np.trapz(p_w_true, x=binCenters)

    int_phot = np.trapz(p_w_phot * DCs_Ds, x=binCenters)
    int_true = np.trapz(p_w_true * DCs_Ds, x=binCenters)
    return int_true / int_phot


def delta_z(z_spec, z_phot):
    return z_spec - z_phot


def delta_z_1pz(z_spec, z_phot):
    return delta_z(z_spec, z_phot) / (1 + z_spec)


def sigma_68(arr, axis=None):
    """Input: an (multi-dimensional) array
    Optional input: the axis along which to calculate the metric
    Outputs: the 68% spread of data about the median value of the array
    """
    upper, lower = np.percentile(arr, [84.075, 15.825], axis=axis)

    return (upper - lower) / 2.0


def outFrac_2sigma68(arr, axis=None):
    """Input: an (multi-dimensional) array
    Optional input: the axis along which to calculate the metric
    Outputs: the fraction of data with more than 2*68% spread of data
    """
    sigma68 = sigma_68(arr, axis=axis)
    return np.sum(np.abs(arr) > 2 * sigma68) * 1.0 / len(arr)


def outFrac_3sigma68(arr, axis=None):
    """Input: an (multi-dimensional) array
    Optional input: the axis along which to calculate the metric
    Outputs: the 68% spread of data with more than 3*68% spread of data
    """
    sigma68 = sigma_68(arr, axis=axis)
    return np.sum(np.abs(arr) > 3 * sigma68) * 1.0 / len(arr)


def mad(arr, axis=None):
    mad_ = np.median(np.abs(arr - np.median(arr)))
    return mad_


def mean_mad(arr, axis=None):
    mad_ = np.mean(np.abs(arr - np.mean(arr)))
    return mad_


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


def dist_pdf(arr, binCenters):
    if len(np.shape(arr)) > 1:
        pdfs = np.array([gaussian_kde(arr[i]).evaluate(binCenters) for i in np.arange(len(arr))])
        return pdfs
    else:
        pdfs = gaussian_kde(arr).evaluate(binCenters)
    return pdfs


def dist_pdf_weights(arr, binCenters, weights=None):
    """ arr is an array of data e.g. [0,1,2,3,4,4,4,,3,2,1,2] to turn into a [weighted] KDE distribution.
    The values of the distribution at binCenters are returned
    """

    if len(np.shape(arr)) > 1:
        pdfs = np.array([gss_kde(arr[i], weights=weights).evaluate(binCenters) for i in np.arange(len(arr))])
        return pdfs
    else:
        pdfs = gss_kde(arr, weights=weights).evaluate(binCenters)
    return pdfs


def normalize_pdf(pdf, z):
    """
    returns normalized pdf
    """
    area = np.trapz(pdf, x=z)
    return pdf / area


def log_loss(act, pred):
    """https://www.kaggle.com/wiki/LogarithmicLoss"""
    epsilon = 1e-15
    import scipy as sp
    pred = sp.maximum(epsilon, pred)
    pred = sp.minimum(1 - epsilon, pred)
    ll = sum(act * sp.log(pred) + sp.subtract(1, act) * sp.log(sp.subtract(1, pred)))
    ll = ll * -1.0/len(act)
    return ll


def kulbachLeiber_bins(Q, P):
    from scipy.stats import entropy
    """Kullbach- Leibler test for binned [strictly >0] distributions
    See en.wikipedia.org/wiki/Kullback-Leibler_divergence
    For P=measured Q=True distribtions"""

    return entropy(P, Q)


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


def stackpdfs(pdfs, weights=None):
    """ numpy shape np.array( (galaxy, pdfbins)) """

    if weights is not None:
        pdfs_ = np.zeros_like(pdfs)
        weights_ = weights / np.sum(weights, dtype=float)
        for i in np.arange(len(weights)):
            pdfs_[i] = pdfs[i]*weights_[i]

        stk_pdf = np.sum(pdfs_, axis=0)
    else:
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

from statsmodels.distributions.empirical_distribution import ECDF
def cumaltive_to_point(dfs, bincenter, points, N=10000):

    """ Determines the y-axis value (between 0-1) of a CDF determined from a pdf (df) evaluated at the x-axis values given by points
    Expected shape: numpy dfs shape (galaxy, z-bins),  
    bincenter center of dfs bins!

    """


    if len(np.shape(dfs)) > 1:
        point_ = points
        if len(np.shape(points)) == 0:
            point_ = [points] * len(dfs)
        #do some iterative magick!
        xarr = np.array([cumaltive_to_point(dfs[c], bincenter, point_[c]) for c in range(len(dfs))])
        return xarr
    else:

        z_mc = get_mc(dfs, bincenter, N=N)
        ec = ECDF(z_mc)
        c = ec(points)
        return c

def xval_cumaltive_at_ypoint(dfs, bincenter, point, k=3):
    """returns the xvals of dfs(ngal, nxbins), x(bins in xdir) point(1), for a point
    which sits on the y-axis of the cdf. We interpolate the cdf for precision"""

    """Note: all points < x[0]-dx are set to x[0] - dx"""
    """Note: all points > x[-1]+dx are set to x[-1] + dx"""


    if len(np.shape(dfs)) > 1:
        #do some iterative magick!
        xarr = np.array([xval_cumaltive_at_ypoint(c, bincenter, point) for c in dfs])
        return xarr
    else:

        #df value sits in a bin, with a given bin center.
        #append to first and last bin a df value to encompass the bin
        delta_bin = (bincenter[1]-bincenter[0]) / 2.0

        binEdges = np.append(bincenter[0] - delta_bin, bincenter)
        binEdges = np.append(binEdges, bincenter[-1] + delta_bin)

        dfs_bins = np.append(np.append(0, dfs), 0)

        point = np.amin((point, np.amax(binEdges)))
        point = np.amax((point, np.amin(binEdges)))

        cum = np.cumsum(dfs_bins)

        #ensure cum==1 at end
        cum = cum / float(cum[-1])

        #note spline interpolation has crazy results!
        return interpolate.interp1d(cum, binEdges)(point)


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
    The x is bin centers of dfs
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
    """ take a binning vector, and bins values, and calculate the pdf-point statistic. Use Boot strap resampling to use weights keyword. If weights are given, we use bootstrap to estimate the weighted average value of the metric

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

        #do we need to bother reweighting?
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


def binned_statistic_dist1_dist2(arr_, bin_vals_, truths_, pdf_, pdf_z_center_, func_, weights=None):
    """ bins the stacked pdf and truths in bins of arr, then calculates the metric on these distributions"""
    """ metric doesn't work so well with weights at the mo"""

    if weights is None:
        weights = np.ones(len(arr_))

    p = weights / np.sum(weights)

    res = {}
    for i in np.arange(len(bin_vals_)-1):
        ind_ = (arr_ >= bin_vals_[i]) * (arr_ < bin_vals_[i + 1])
        if np.sum(ind_) > 0:
            res[i] = {}

            truth_dist_ = gss_kde(truths_[ind_], weights=p[ind_]).evaluate(pdf_z_center_)

            stacked_pdf_ = stackpdfs(pdf_[ind_], weights=p[ind_])
            stacked_pdf_ = normalisepdfs(stacked_pdf_, pdf_z_center_)

            res[i]['weighted_bin_center'] = np.average(arr_[ind_], weights=p[ind_])

            """ Add weights in here  """
            res[i]['weighted_value'] = func_(truth_dist_, stacked_pdf_)

    return res


""" --- tools for extracting properties from pdfs ----"""
def get_mc(pdf_, binCenter, N=1):
    import random as rdm
    # renorm incase there is probability at higher-z that we've cut, or some error.
    if np.sum(pdf_) > 0:
        pdf = pdf_ / np.sum(pdf_)
        #print pdf
        #rnd = int(rdm.random() * (np.power(2, 31)))
        #np.random.seed(rnd)
        zbn_center = np.random.choice(binCenter, p=pdf, size=N)
        dz = (binCenter[1] - binCenter[0])
        return zbn_center + (np.random.uniform(size=N) - 0.5) * dz
    else:
        return np.nan

def get_mean(pdf, zarr):
    if np.sum(pdf) > 0:
        zm = np.average(zarr, weights=pdf)
        return zm
    else:
        return -1.

def get_sig(pdf, zarr):
    if np.sum(pdf) > 0:
        zm = np.average(zarr, weights=pdf)
        sig = np.sqrt(np.average((zarr-zm)*(zarr-zm), weights=pdf))
        return sig
    else:
        return -1.

def get_mean_and_sig(pdf, zarr):
    if np.sum(pdf) > 0:
        zm = np.average(zarr, weights=pdf)
        sig = np.sqrt(np.average((zarr-zm)*(zarr-zm), weights=pdf))
        return zm, sig
    else:
        return -1.

def get_median(pdf, zarr):
    if np.sum(pdf) > 0:
        return xval_cumaltive_at_ypoint(pdf, zarr, 0.5)
    else:
        return -1.

def get_sig68(pdf, zarr):
    s2 = xval_cumaltive_at_ypoint(pdf, zarr, 0.84075)
    s1 = xval_cumaltive_at_ypoint(pdf, zarr, 0.15825)
    s68 = (s2 - s1) / 2.0
    return s68

#to do, uniti test these bad boys!
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


