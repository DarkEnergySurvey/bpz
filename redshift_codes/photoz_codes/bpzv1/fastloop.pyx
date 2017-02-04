#cython: boundscheck=False, wraparound=False, nonecheck=False
from libc.math cimport exp 
import numpy as np
import bh_photo_z_validation as pval
import random as rdm

def bpz_loop(int n_gals, double [:] z_bins, double[:, :] f_obs,  double[:] ef_obs, double[:, :, :] f_mod , double[:] prior_mag, GALPRIOR, template_type_dict):

    cdef int i, j, nf

    nf = len(f_mod[0, :, 0])

    cdef double[:] z_max_post = np.zeros(n_gals)
    cdef double[:] mean = np.zeros(n_gals)
    cdef double[:] sigma = np.zeros(n_gals)
    cdef double[:] median = np.zeros(n_gals)
    cdef double[:] mc = np.zeros(n_gals)
    cdef double[:] sig68 = np.zeros(n_gals)
    #cdef double[:] chi2 = np.zeros(n_gals, nf)

    eps = 1e-300
    eeps = np.log(eps)

    for i in range(n_gals):

        f = f_obs[i]
        ef = ef_obs[i]
        foo = np.sum(np.power(f / ef ,2))
        fot = np.sum(f.reshape(1, 1, nf) * f_mod / (ef.reshape(1, 1, nf) ** 2.0), axis=2)

        ftt = np.sum(np.power(f_mod, 2) / np.power(ef.reshape(1, 1, nf), 2.0), axis=2)
        chi2 = foo - np.power(fot, 2) / (ftt + eps)

        ind_mchi2 = np.where(chi2 == np.amin(chi2))
        min_chi2 = chi2[ind_mchi2][0]

        z_min_chi2 = z_bins[ind_mchi2[0]]

        likelihood = exp(-0.5 * np.clip(chi2 - min_chi2, 0., -2 * eeps))

        prior = np.zeros_like(likelihood)
        for j in range(len(f_mod[0, :, 0])):
            prior[:, j] = GALPRIOR.evaluate(prior_mag[i], 
                           frac_prior_from_types=template_type_dict[j])

        #posterior is prior * Likelihood
        posterior = prior * likelihood

        #margenalise over Templates:
        marg_post = np.sum(posterior, axis=1)
        marg_post /= np.sum(marg_post)
        
        ind_max_marg = np.where(marg_post == np.amax(marg_post))[0][0]

        mean[i] = get_mean(marg_post, z_bins)
        sigma[i] = get_sig(marg_post, z_bins)
        median[i] = get_median(marg_post, z_bins)
        mc[i] = get_mc(marg_post, z_bins)
        sig68[i] = get_sig68(marg_post, z_bins)

        z_max_post[i] = z_bins[ind_max_marg]

    cols = {'MEAN_Z': mean, 'Z_SIGMA': sigma, 'MEDIAN_Z': median,
           'Z_MC': mc, 'Z_SIGMA68':sig68, 'z_max_post': z_max_post}
    return cols



def get_mc(pdf, zarr):
    # renorm incase there is probability at higher-z that we've cut, or some error.
    if np.sum(pdf) > 0:
        targ_prob = rdm.random()
        return pval.xval_cumaltive_at_ypoint(pdf, zarr, targ_prob)
    else:
        return -1.

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
        return pval.xval_cumaltive_at_ypoint(pdf, zarr, 0.5)
    else:
        return -1.

def get_sig68(pdf, zarr):
    s2 = pval.xval_cumaltive_at_ypoint(pdf, zarr, 0.84075)
    s1 = pval.xval_cumaltive_at_ypoint(pdf, zarr, 0.15825)
    s68 = (s2 - s1) / 2.0
    return s68

