#! /usr/bin/env python

# Gneiting, T. and Raftery, A. E. (2004).
import copy
import numpy as np
import random as rdm
import sys
from joblib import Parallel, delayed
import bh_photo_z_validation as pval
from scipy.stats import entropy
import os 

def _help():
    print ("call like: ")
    print ("pdf_recalibration.py PathToPDF.h5 out=OutputRSfile.p true_z=XXX bin_col=YYY")
    print ("or...")
    print ("pdf_recalibration.py InputRSfile.p PathToPDF*.h5")
    print (
    "if true_z and bin_col are given PathToPDF.h5  must contain /point_predictions/ with true_z e.g. true_z=Z_SPEC|REDSHIFT")
    print (
    "else, will load InputRSfile.p and apply the pdf rescaling to PathToPDF*.hdf5 and produce a new .h5 files with PathToPDF*.RS.h5")
    sys.exit()


class Parrallelise:
    """ Creates a generic method to perform
    trivially Parrallel operations
    loop is a tuple of things that one wants to use in the declared function
    call like:
    from bh_parallelise import Parrallelise
    p = Parrallelise(n_jobs=2, method=myFunction, loop=[Tuple,Tuple,..]).run()
    Tuple is passed to myFunction and containes all the info required by the function

    Exampe call:

    #load all the files
    files = glob.glob('/Volumes/Untitled/DES/Y1_GOLD_V1_0_1/*.csv')

    arr = []
    for n, f in enumerate(files):
        arr.append(['ALL_DES_SPECTRA.fits', f])

    res1 = Parrallelise(n_jobs=5, method=stiltsMatch, loop=arr).run()

    """

    def __init__(self, n_jobs=-1, method=None, loop=None):
        self._n_jobs = n_jobs
        self._method = method
        self._loop = loop

    def run(self):
        results = Parallel(n_jobs=self._n_jobs)(delayed(self._method)(l) for l in self._loop)
        return results



def apply_rescalingv1(z, zbin_centers_, OffGaussCauch_, numpy_seed=None):
    """pply_rescaling(zarry, zbin_centers, GaussCauch_, NG=500) 
    Take an intial DF, normalise, sample with gaussian smoothing and add a cauchy term to extend tails
    pd_ = input df
    zbinEdges_ = start bins of pd_ [and final bin len +1 compared to zbin_centers and pdf]
    GaussCauch_ = [GaussSmoothingFactor, Cauchy scale parameter, 
    see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cauchy.html#scipy.stats.cauchy ]
    NG = Number of random draws from pd_ before smoothing
    
    """

    from scipy.stats import cauchy

    # set random seed, otherwise differnt results each time it's called with same OffGaussCauch_
    # during minimisation of KL
    np.random.seed(numpy_seed)
    # get bin centers
    z_ = copy.copy(z)

    # sample from this df with probability pd, with some additional smoothing
    z_ += np.random.normal(size=len(z_)) * OffGaussCauch_[1] + OffGaussCauch_[0]

    # ignore unphysical z<0
    z_ = z_[z_ > 0]
    dz = 0.5 * (zbin_centers_[1] - zbin_centers_[0])
    zbinEdges_ = np.append(zbin_centers_ - dz, zbin_centers_[-1] + dz)
    h = np.histogram(z_, bins=zbinEdges_)
    h = h[0] * 1.0
    h[h < 0] = 0
    h = h / np.trapz(h, x=zbin_centers_)
    # generate a gaussian KDE  and atdd a cauchy pdf
    pd = h * OffGaussCauch_[2] + cauchy.pdf(zbin_centers_, loc=np.median(z_),
                                            scale=OffGaussCauch_[3]) * (1 - OffGaussCauch_[2])
    pd[pd < 0] = 0
    # this is not quite correct because we must intergrate from 0 .. max(zbinEdges_)
    # not zbin_centers_, but we don't really know 0--zbin_centers_[1] width, so bah!
    pd = pd / np.trapz(pd, x=zbin_centers_)
    return pd


def make_histogram(c_, N_hist_bins=101):
    """make_histogram(c_, N_hist_bins=101)
    converts a distrbution c_ into a histogram, binned between [0, 1] with N_hist_bins bins
    return 1d array of histograms
    """
    # ignore NaNs (possible due to a flat pdf)
    ind = np.isfinite(c_) == True

    hist = np.histogram(c_[ind], bins=np.linspace(0, 1, num=N_hist_bins, endpoint=True))

    return hist[0], hist[1]


def parr_rescling(inp):
    
    zarr_, zbin_centers_, OffGaussCauch_, z_spec_ = inp
    c = np.zeros(len(z_spec_))
    for i in range(len(z_spec_)):
        # get trial PDFs
        pd = apply_rescalingv1(zarr_[i], zbin_centers_, OffGaussCauch_, numpy_seed=32)
        # calcaulte y-axis value of pdf
        c[i] = pval.cumaltive_to_point(pd, zbin_centers_, z_spec_[i])
    return c


def learn_rescalingv1(OffGaussCauch_, zarr_, z_spec_, zbin_centers_):
    """learn_rescaling(OffGaussCauch, pdfs_, z_spec_, zbin_centers_)
    Takes a set of dfs and returns a proposed set of pdfs
    Then asks is the proposed pdfs consistent with the distributions of z_spec_ truths 
    GaussCauch =  [GaussSmoothingFactor, Cauchy scale parameter, see 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.cauchy.html]
    pdfs_ = n-d array of dfs pdfs_[0] = 0th pdf
    z_spec_ = truth values 
    zbin_centers_ = z bin centers -- we infer bin edges
    
    returns Kullbach Liebler divergence between distributions of histogramed CDF values of 
    targets, and the average value of the histogram
    0 = perfect match == great rescaling
    >0 = imperfct match , pdfs not really pdfs. Use with care.
    
    """

    # hold the y-axis values of CDF at target
    c = np.zeros(len(z_spec_))

    inds = np.array_split(np.arange(len(z_spec_)), int(len(z_spec_) / 8) + 2)

    parr = []
    for ind in inds:
        parr.append([zarr_[ind], zbin_centers_, OffGaussCauch_, z_spec_[ind]])

    c_res = Parrallelise(n_jobs=6, method=parr_rescling, loop=parr).run()
    for i, ind in enumerate(inds):
        c[ind] = c_res[i]

    # turn this to a histogram
    hist, _ = make_histogram(c, 51)
    # scipy.stats.entropy(pk, qk=None, base=None) = KL divergence between histograms of C and <hist C>
    return entropy(hist, [np.mean(hist)] * len(hist))


def learn_rescaling(fil, bin_col, z_true_col, Nspec_per_bin=8000, Nsamples_pdf=2000):
    """take a pandas file, in the format expected by the photo-z wg, and produce learn a rescaling, for objects in bins given by bin_col. Bins are defined to have Nspec_per_bin in them. The truth column is given by z_true_col

    """
    import pandas as pd
    from scipy.optimize import minimize

    df = pd.read_hdf(fil, 'point_predictions')

    # fixed number of objects with spectra per bin
    mg = np.sort(df[bin_col])
    bin_limits = np.unique(np.append(mg[np.array(np.arange(len(mg) * 1.0 / Nspec_per_bin) * Nspec_per_bin, dtype=int)], np.amax(mg)))
    result_z = {'bin_limits': bin_limits, 'bin_col': bin_col}
    z_spec = np.array(df[z_true_col])

    del mg

    # get bin data, and redshift information
    bns = pd.read_hdf(fil, 'info')
    zbin_center = np.array(bns['z_bin_centers'])

    dm = (zbin_center[1] - zbin_center[0]) / 2.0

    # extract the pdfs
    pdf = pd.read_hdf(fil, 'pdf_predictions')
    pdfc = [i for i in pdf.keys() if 'pdf_' in i]
    pdfs = pdf[list(np.sort(pdfc))].as_matrix()

    # results for unscaled pdfs
    c = np.zeros(len(z_spec))

    #generate draws from pdf
    z_array = np.zeros((len(z_spec), Nsamples_pdf))
    pdf_non_0 = np.zeros(len(z_spec), dtype=bool)
    for i in range(len(z_spec)):

        pd_ = pdfs[i] / np.sum(pdfs[i])
        if np.all(np.isfinite(pd_)):
            zarr = pval.get_mc(pd_, zbin_center, N=Nsamples_pdf)
            c[i] = pval.cumaltive_to_point(pd_, zbin_center, z_spec[i])
            z_array[i] = zarr
            pdf_non_0[i] = True

    res_ = minimize(learn_rescalingv1,
                    [ -2.35634796e-04, 7.32995542e-02, 1.00284366e+00, 6.24823595e-04],
                    args=(z_array[pdf_non_0, :], z_spec[pdf_non_0], zbin_center),
                    method='Nelder-Mead',
                    tol=1e-6, options={'maxiter': 500}
                    )
    print res_.x

    for i in range(len(bin_limits) - 1):
        # identify which data sit in the binned column
        ind = np.arange(len(df))[np.array(df[bin_col] < bin_limits[i + 1]) * np.array(df[bin_col] > bin_limits[i]) * pdf_non_0]

        np.random.shuffle(ind)
        tst = ind[0:len(ind) / 2]
        trn = ind[len(ind) / 2:]

        # turn this to a histogram
        hist, _ = make_histogram(c[trn], 51)
        # scipy.stats.entropy(pk, qk=None, base=None) = KL divergence between histograms of C and <hist C>
        trn_kl_b4 = entropy(hist, [np.mean(hist)] * len(hist))

        # turn this to a histogram
        hist, _ = make_histogram(c[tst], 51)
        # scipy.stats.entropy(pk, qk=None, base=None) = KL divergence between histograms of C and <hist C>
        tst_kl_b4 = entropy(hist, [np.mean(hist)] * len(hist))

        res = minimize(learn_rescalingv1, res_.x,
                       args=(z_array[trn], z_spec[trn], zbin_center),
                       method='Nelder-Mead', tol=1e-6, options={'maxiter': 500})

        trn_kl = learn_rescalingv1(res.x, z_array[trn], z_spec[trn], zbin_center)
        tst_kl = learn_rescalingv1(res.x, z_array[tst], z_spec[tst], zbin_center)
        result_z[i] = {'MAX_BIN': bin_limits[i + 1], 'MIN_BIN': bin_limits[i]}

        result_z[i]['res'] = res
        result_z[i]['KL_TRAIN'] = [trn_kl_b4, trn_kl]
        result_z[i]['KL_TEST'] = [tst_kl_b4, tst_kl]
        print ' '
        print i , {'MAX_BIN': bin_limits[i + 1], 'MIN_BIN': bin_limits[i]}
        print result_z[i]['KL_TRAIN']
        print result_z[i]['KL_TEST']
    return result_z


def apply_rescaling(fil, result_z, Nsamples_pdf=8000):
    """apply rescaling to a pdf.h5 as saved in the standard DES format"""

    df = pd.read_hdf(fil, 'point_predictions')
    bns = pd.read_hdf(fil, 'info')
    zbin_center = np.array(bns['z_bin_centers'])
    dm = (zbin_center[1] - zbin_center[0]) / 2.0

    # extract the pdfs
    pdf = pd.read_hdf(fil, 'pdf_predictions')
    pdfc = [i for i in pdf.keys() if 'pdf_' in i]
    pdfs = pdf[list(np.sort(pdfc))].as_matrix()

    new_pdf = np.zeros_like(pdfs)
    bin_col_vals = np.array(df[result_z['bin_col']])

    n_gals = len(pdfs)
    #results arrays for this loop
    mode = np.zeros(n_gals) + np.nan
    mean = np.zeros(n_gals) + np.nan
    sigma = np.zeros(n_gals) + np.nan
    median = np.zeros(n_gals) + np.nan
    mc = np.zeros(n_gals) + np.nan
    sig68 = np.zeros(n_gals) + np.nan

    for i in range(n_gals):
        ind = np.where(result_z['bin_limits'] >= bin_col_vals[i])[0][0]
        ind = np.amax([ind - 1, 0])

        pd_ = pdfs[i] / np.sum(pdfs[i])
        if np.all(np.isfinite(pd_)):
            zarr = pval.get_mc(pd_, zbin_center, N=Nsamples_pdf)

            pd1 = apply_rescalingv1(zarr, zbin_center, result_z[ind]['res'].x)
            new_pdf[i] = pd1

            #measure point predictions from pdf
            marg_post = new_pdf[i] / np.sum(new_pdf[i])
            if np.all(np.isfinite(marg_post)):
                ind_max_marg = np.where(marg_post == np.amax(marg_post))[0]
                if len(ind_max_marg) > 1:
                    ind_max_marg = np.random.choice(ind_max_marg)
                ind_max_marg = ind_max_marg[0]
                mean[i] = pval.get_mean(marg_post, zbin_center)
                sigma[i] = pval.get_sig(marg_post, zbin_center)
                median[i] = pval.get_median(marg_post, zbin_center)
                mc[i] = pval.get_mc(marg_post, zbin_center)
                sig68[i] = pval.get_sig68(marg_post, zbin_center)
                mode[i] = zbin_center[ind_max_marg]

    return {'pdf': new_pdf, 'MEAN_Z': mean, 'Z_SIGMA': sigma,
            'MEDIAN_Z': median, 'Z_SIGMA68': sig68, 'MODE_Z': mode,
            'Z_MC': mc
            }

if __name__ == '__main__':
    args = sys.argv[1:]

    if len(args) < 2:
        _help()

    inArgs = {}
    files = []
    for i in args:
        if '=' in i:
            k, v = i.split('=')
            inArgs[k] = v
        else:
            files.append(i)

    import cPickle as pickle
    import pandas as pd

    if len(inArgs.keys()) > 0:
        #learn the rescaling and save it to a file.
        if ('true_z' not in inArgs) or ('bin_col' not in inArgs):
            _help()
        # learn rescaling
        result_z = learn_rescaling(files[0], inArgs['bin_col'], inArgs['true_z'])
        pickle.dump(result_z, open(inArgs['out'], 'w'))

    else:
        #apply the rescaling to a set of file[s]
        rescal_file = files[0]
        files = files[1:]

        result_z = pickle.load(open(rescal_file, 'r'))

        # check for success! give a massive warning if not converged!!!
        print result_z.keys()
        for i in range(len(result_z['bin_limits']) - 1):
            print (i, result_z[i]['res'].success, result_z['bin_limits'][i], result_z[i]['KL_TEST'])

        for fil in files:
            rs_file = fil.split('.h')[0] + 'ReScaled.h5'

            if os.path.isfile(rs_file):
                os.remove(rs_file)

            resc_pdfs = apply_rescaling(fil, result_z, Nsamples_pdf=4000)
            #{'pdf': new_pdf, 'MEAN_Z': mean, 'Z_SIGMA': sigma,
            #'MEDIAN_Z': median, 'Z_SIGMA68': sig68, 'MODE_Z': mode,
            #'Z_MC': mc
            #}

            df = pd.DataFrame()
            df.to_hdf(rs_file, 'pdf_predictions', append=True)
            df.to_hdf(rs_file, 'point_predictions', append=True)
            df.to_hdf(rs_file, 'info', append=True)

            #copy the info
            df2 = pd.read_hdf(fil, 'info')
            df2.to_hdf(rs_file, key='info', format='table', append=True, complevel=5, complib='blosc')
            
            #now populate the point predictions
            cols = {}
            point_keys = [j for j in resc_pdfs.keys() if j != 'pdf']
            for i in point_keys:
                cols[i] = resc_pdfs[i]

            #get old columns to copy accross
            old_file = pd.read_hdf(fil, 'point_predictions')

            for i in old_file.keys():
                if i not in cols:
                    if i not in ['KL_POST_PRIOR', 'TEMPLATE_TYPE', 'MINCHI2']:
                        cols[i] = old_file[i]

            df2 = pd.DataFrame(cols)
            df2.to_hdf(rs_file, key='point_predictions', append=True)

            #now save the rescaled pdfs
            old_file = pd.read_hdf(fil, 'pdf_predictions')

            #start populating the resampled pdf file
            post_dict = {'MEAN_Z': resc_pdfs['MEAN_Z']}

            df2 = pd.read_hdf(rs_file, 'info')
            z_bin_centers = np.array(df2['z_bin_centers'])

            for ii in np.arange(len(z_bin_centers)):
                post_dict['pdf_{:0.4}'.format(z_bin_centers[ii])] = resc_pdfs['pdf'][:, ii]
            print post_dict.keys()

            pdfc = [i for i in old_file.keys() if 'pdf_' not in i]
            pdfc = [i for i in pdfc if i not in ['KL_POST_PRIOR', 'TEMPLATE_TYPE', 'MINCHI2']]

            #any additional columns
            for i in pdfc:
                if i not in resc_pdfs:
                    post_dict[i] = old_file[i]

            df2 = pd.DataFrame(post_dict)
            df2.to_hdf(rs_file, key='pdf_predictions', append=True,
                       complevel=5, complib='blosc')

            print ("rescaled file written to ", rs_file)
            #free memory
            del df2
