#! /usr/bin/env python

"""BPZ in 2017

author: benhoyle1212@gmail.com

about: This codes implements some of the key features in the original Bayesean Probability Redshift BPZ (Benitez 200X), with a more transparent user interface, and the native .fits support.

    The original BPZ is still very useful. e.g determining photometric offsets is not implemented

To do:
Deal with missing / unseen data properly

"""

import sys
from joblib import Parallel, delayed
import numpy as np
import inspect
import bh_photo_z_validation as pval
import random as rdm


def _help():
    print "new version BPZ"
    print "call like"
    print "bpzv1.py pathToConfigFile.yaml pathTofits*.fits"
    writeExampleConfig()
    sys.exit()


def writeExampleConfig():
    import os.path
    import os
    import textwrap
    path = os.getcwd() + '/'
    if os.path.isfile(path + 'exampleBPZConfig.yaml') is False:
        f = open(path + 'exampleBPZConfig.yaml', 'w')
        txt = textwrap.dedent("""
#redshift bins min, max, width
redshift_bins: [0.01, 3.5, 0.01]

#either set here, or we will determine this from the code.
#Default: determine from code
BPZ_BASE_DIR: 
AB_DIR: 
SED_DIR: 
FILTER_DIR: 

#spectra list. This *must* match the sed_type below.
sed_list: [El_B2004a.sed, Sbc_B2004a.sed, Scd_B2004a.sed, El_B2004a.sed, Im_B2004a.sed, SB3_B2004a.sed, El_B2004a.sed, SB2_B2004a.sed]

#Either E/S0 Spiral or Irr (elliptical, spiral, Irregular). The SEDs will be rearranged such that their is exact ascending order E/S0->Spiral->Irr
sed_type: [E/S0, Spiral, Spiral, E/S0, Irr, Irr, E/S0, Irr]

#go crazy and reorder all spectra types? Note: this is properly unphysical,
#due to interpolation reordering above!
rearrange_spectra: False

# prior name, any set you like. See sed_proir_file.py for details.
prior_name: sed_prior_file.cosmos_Laigle

#expect i-band mag. e.g. MAG_AUTO_I
PRIOR_MAGNITUDE: MAG_I

#work with MAGS [True] or FLUX [False]. If left blank [default] the code infers this 
#from the presence of MAG or mag in the XXX of filters: ky: {MAG_OR_FLUX: 'XXX'}
INPUT_MAGS: 

#minimum magnitude error
MIN_MAGERR: 0.001

# Objects not observed
mag_unobs: -99 

#Objects not detected 
mag_undet: 99

#this construct which magnitudes / or FLUXES map to which filters
filters: {
    DECam_2014_g.res: {MAG_OR_FLUX: MAG_G, ERR: MAGERR_G, AB_V: AB, zp_error: 0.02, zp_offset: 0.0},
    DECam_2014_r.res: {MAG_OR_FLUX: MAG_R, ERR: MAGERR_R, AB_V: AB, zp_error: 0.02, zp_offset: 0.0},
    DECam_2014_i.res: {MAG_OR_FLUX: MAG_I, ERR: MAGERR_I, AB_V: AB, zp_error: 0.02, zp_offset: 0.0},
    DECam_2014_z.res: {MAG_OR_FLUX: MAG_Z, ERR: MAGERR_Z, AB_V: AB, zp_error: 0.02, zp_offset: 0.0}
    #DECam_2014_Y.res: {MAG: MAG_MOF_Y, ERR: MAGERR_MOF_Y, AB_V: AB, zp_error: 0.02, zp_offset: 0.0}
    }


#this construct which magnitudes / or FLUXES map to which filters
filters: {
    DECam_2014_g.res: {MAG_OR_FLUX: FLUX_G, ERR: FLUX_ERR_G, AB_V: AB, zp_error: 0.02, zp_offset: 0.0},
    DECam_2014_r.res: {MAG_OR_FLUX: FLUX_R, ERR: FLUX_ERR_R, AB_V: AB, zp_error: 0.02, zp_offset: 0.0},
    DECam_2014_i.res: {MAG_OR_FLUX: FLUX_I, ERR: FLUX_ERR_I, AB_V: AB, zp_error: 0.02, zp_offset: 0.0},
    DECam_2014_z.res: {MAG_OR_FLUX: FLUX_Z, ERR: FLUX_ERR_Z, AB_V: AB, zp_error: 0.02, zp_offset: 0.0}
    #DECam_2014_Y.res: {MAG: MAG_MOF_Y, ERR: MAGERR_MOF_Y, AB_V: AB, zp_error: 0.02, zp_offset: 0.0}
    }


#which magnitude will we use for flux normalisation?
normalisation_filter: DECam_2014_i.res

#this is the id column. Don't mess around! use it.
ID: COADD_OBJECTS_ID

#if these columns [in the stated case] don't exist a warning will be made, but the code will run.
ADDITIONAL_OUTPUT_COLUMNS: [REDSHIFT, R11, R22, MAG_I, MAGERR_I]

#do you wanna output a suffix for a filename
output_file_suffix:
#do we also want pdfs to be produced?
output_pdfs: True

#N_INTERPOLATE_TEMPLATES: Blank means No
INTERP: 8

#should we parralise the loops?
n_jobs: 5

#print some information to screen
verbose: True 
""")
        f.write(txt)
        print ("An example file exampleBPZConfig.yaml has been written to disk")

def key_not_none(d, ky):
    if ky in d:
        return d[ky] is not None
    return False

def load_file(file_name):
    z, f = [], []
    for i in open(file_name, 'r'):
        tmp = i.replace('  ', ' ').replace('\n', '').split(' ')
        z.append(float(tmp[0]))
        f.append(float(tmp[1]))
    return np.array(z), np.array(f)

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


def e_mag2frac(errmag):
    """Convert mag error to fractionary flux error"""
    return 10. ** (.4 * errmag) - 1.


"""
if __name__ == '__main__':
    args = sys.argv[1:]

    if len(args) < 2 or 'help' in args or '-h' in args:
        _help()
"""


def parr_loop(lst):
    """parralisation loop for joblib
    lst is a [] containing ind, f_obs_, ef_obs_, prior_mag_, f_mod, gal_mag_type_prior, z_bins, config 
    f_mod = tempalate fluxes [redshift, template, band]

    returns a dictionary with
    {'ind': ind_, 'mean': mean, 'sigma': sigma, 'median': median, 'sig68': sig68, 'z_max_post': z_max_post,
            'KL_post_prior': KL_post_prior, 'pdfs_': pdfs_}
    for later combination
    """

    from scipy.stats import entropy

    ind_, f_obs_, ef_obs_, prior_mag_, f_mod, gal_mag_type_prior, z_bins, config = lst

    n_gals = len(ind_)
    #some small value to truncate probs.
    eps = 1e-300
    eeps = np.log(eps)

    #results arrays for this loop
    z_max_post = np.zeros(n_gals) + np.nan
    mean = np.zeros(n_gals) + np.nan
    sigma = np.zeros(n_gals) + np.nan
    median = np.zeros(n_gals) + np.nan
    mc = np.zeros(n_gals) + np.nan
    sig68 = np.zeros(n_gals) + np.nan
    KL_post_prior = np.zeros(n_gals) + np.nan
    min_chi2_arr = np.zeros(n_gals) + np.nan
    maxL_template_ind = np.zeros(n_gals, dtype=int) - 1000
    pdfs_ = np.zeros((n_gals, len(z_bins))) + np.nan

    for i in np.arange(n_gals):

        foo = np.sum(np.power(f_obs_[i] / ef_obs_[i], 2))

        nf = len(f_mod[0, 0, :])

        f = f_obs_[i].reshape(1, 1, nf)
        ef = ef_obs_[i].reshape(1, 1, nf)

        #this is a slow part of the code!
        fot = np.sum(np.divide(f * f_mod,  np.power(ef, 2)), axis=2)

        #this is the slowest part of the code!
        ftt = np.sum(np.power(np.divide(f_mod, ef), 2), axis=2)

        chi2 = foo - np.power(fot, 2) / (ftt + eps)

        ind_mchi2 = np.where(chi2 == np.amin(chi2))

        min_chi2 = chi2[ind_mchi2][0]

        min_chi2_arr[i] = min_chi2

        z_min_chi2 = z_bins[ind_mchi2[0]][0]

        likelihood = np.exp(-0.5 * np.clip(chi2 - min_chi2, 0., -2 * eeps))

        prior = np.zeros_like(likelihood)
        pr_mg = gal_mag_type_prior.keys()

        ind_mag_p = np.argmin(np.abs(prior_mag_[i] - np.array(pr_mg)))

        for j in np.arange(len(f_mod[0, :, 0])):
            prior[:, j] = gal_mag_type_prior[pr_mg[ind_mag_p]][j]

        #posterior is prior * Likelihood
        posterior = prior * likelihood

        #get the maximum posterior, and determine which template this is
        ind_max = np.where(posterior == np.amax(posterior))[1]

        #if many "best maf posteriors" then choose one at random.
        if len(ind_max) > 1:
            ind_max = np.random.choice(ind_max)
        maxL_template_ind[i] = ind_max

        #margenalise over Templates in Prior and posterior:
        marg_post = np.sum(posterior, axis=1)
        marg_post /= np.sum(marg_post)
        marg_prior = np.sum(prior, axis=1)
        marg_prior /= np.sum(marg_prior)

        KL_post_prior[i] = entropy(marg_post, marg_prior)

        ind_max_marg = np.where(marg_post == np.amax(marg_post))[0][0]

        #define summary stats from the margenalised posterior.
        mean[i] = pval.get_mean(marg_post, z_bins)
        sigma[i] = pval.get_sig(marg_post, z_bins)
        median[i] = pval.get_median(marg_post, z_bins)
        mc[i] = pval.get_mc(marg_post, z_bins)
        sig68[i] = pval.get_sig68(marg_post, z_bins)
        z_max_post[i] = z_bins[ind_max_marg]

        if key_not_none(config, 'output_pdfs'):
            pdfs_[i] = marg_post

    verbose = key_not_none(config, 'verbose')
    if verbose:
        print ('loop complete', config['n_jobs'])

    return {'ind': ind_, 'mean': mean, 'sigma': sigma, 'median': median, 'sig68': sig68, 'z_max_post': z_max_post,
            'KL_post_prior': KL_post_prior, 'pdfs_': pdfs_, 'mc': mc,
            'min_chi2': min_chi2_arr, 'maxL_template_ind': maxL_template_ind}


def main(args):

    import yaml
    from galaxy_type_prior import GALAXYTYPE_PRIOR
    import os
    import copy
    from astropy.io import fits as pyfits
    import copy
    import time
    import pandas as pd

    #some timing tests
    t = time.time()
    config = yaml.load(open(args[0]))

    verbose = key_not_none(config, 'verbose')

    files = args[1:]
    if isinstance(files, list) is False:
        files = [files]

    #Set up bpz paths, or infer them from the location of this script.
    if key_not_none(config, 'BPZ_BASE_DIR') is False:
        config['BPZ_BASE_DIR'] = '/' + '/'.join([i for i in os.path.realpath(__file__).split('/')[0:-1]]) + '/'

    if key_not_none(config, 'AB_DIR') is False:
        config['AB_DIR'] = config['BPZ_BASE_DIR'] + '../../templates/ABcorr/'

    if key_not_none(config, 'FILTER_DIR') is False:
        config['FILTER_DIR'] = config['BPZ_BASE_DIR'] + 'FILTER/'

    if key_not_none(config, 'SED_DIR') is False:
        config['SED_DIR'] = config['BPZ_BASE_DIR'] + 'SED/'

    if key_not_none(config, 'ID') is False:
        print ("you must provide an ID column in the config file!")
        sys.exit()

    output_file_suffix = ''
    if key_not_none(config, 'output_file_suffix'):
        output_file_suffix = config['output_file_suffix']

    #load redshift bins
    z_bins = np.arange(config['redshift_bins'][0], config['redshift_bins'][1] + config['redshift_bins'][2], config['redshift_bins'][2])

    #load filters.
    filters = config['filters'].keys()

    MAG_OR_FLUX = [config['filters'][i]['MAG_OR_FLUX'] for i in filters]
    MAG_OR_FLUXERR = [config['filters'][i]['ERR'] for i in filters]

    if config['INPUT_MAGS'] is None:
        config['INPUT_MAGS'] = ('mag' in MAG_OR_FLUX[0]) or ('MAG' in MAG_OR_FLUX[0])

    #jazz with spectra_list. Make all combinations of everything
    if config['rearrange_spectra']:
        spectra = []
        sed = []
        for i, s1 in enumerate(config['sed_list']):
            for j, s2 in enumerate(config['sed_list']):
                spectra.append(s1)
                spectra.append(s2)
                sed.append(config['sed_type'][i])
                sed.append(config['sed_type'][j])

        config['sed_list'] = spectra
        config['sed_type'] = sed

    #each sed (in order) is a number between 1 -> Num seds.
    #interpolatated sed are fractional quantites between, e.g. 1.4 = 0.6 of sed1 and 0.4 of sed2 
    sed_float_list = ((np.arange(len(config['sed_list'])) + 1.0) / len(config['sed_list'])) * len(config['sed_list'])

    #identify the filter we will normalise to.
    ind_norm = np.where(np.array(filters) == config['normalisation_filter'])[0][0]

    #get flux offsets, and errors.
    zp_offsets = np.array([config['filters'][i]['zp_offset'] for i in filters])
    zp_offsets = np.power(10, -.4 * zp_offsets)

    zp_error = np.array([config['filters'][i]['zp_error'] for i in filters])
    zp_frac = e_mag2frac(zp_error)


    #prepare the prior
    GALPRIOR = GALAXYTYPE_PRIOR(z=z_bins, tipo_prior=config['prior_name'],
                    mag_bins=np.arange(18, 24.1, 0.1),
                    template_type_list=config['sed_type'])

    """Populate the model fluxes here """
    #prepare to populate the model fluxes
    f_mod = np.zeros((len(z_bins), len(config['sed_list']), len(filters)), dtype=float)

    #which template types does each SED correspond to
    template_type_dict = {}
    for i, sed in enumerate(config['sed_list']):
        sed_ = sed.replace('.sed', '')

        #this dummy dictionary allows for new "galaxy SED types"
        #to be defined in the future.
        dict_ = {}
        for gal_typ in np.unique(config['sed_type']):
            dict_[gal_typ] = 0.0
        dict_[config['sed_type'][i]] = 1.0
        template_type_dict[i] = copy.copy(dict_)

        for j, fltr in enumerate(filters):
            fltr_ = fltr.replace('.res', '')
            file_name = config['AB_DIR'] + '.'.join([sed_, fltr_, 'AB'])
            z_, f_ = load_file(file_name)

            #interpolate to the exact redshifts we care about
            f_mod[:, i, j] = np.interp(z_bins, z_, f_)

            #clip, following original bpz
            f_mod[:, i, j] = np.clip(f_mod[:, i, j], 0, 1e300)

    #interp between SEDs
    if key_not_none(config, 'INTERP'):

        #how many interpolated points?
        num_interps = (len(config['sed_list'])-1) * config['INTERP'] + len(config['sed_list'])

        #generate some dummy indicies that are intergers spaced between 0 -- num_interps
        index = np.linspace(0, num_interps, len(config['sed_list']), endpoint=True, dtype=int)

        #generate values between dummy indicies. These will be interpolated
        ind_int = np.arange(num_interps + 1)

        sed_float_list = np.interp(ind_int, index, sed_float_list)

        #Frist interpolate the galaxy SED types. E.g.
        #{'E/S0': 1 'Spiral': 0 'Irr': 0} -> {'E/S0': 0.5 'Spiral': 0.5 'Irr': 0}
        #->{'E/S0': 0 'Spiral': 1 'Irr': 0}
        template_type_dict_interp_ = {}
        for i in ind_int:
            template_type_dict_interp_[i] = copy.copy(template_type_dict[0])

        for gal_typ in np.unique(config['sed_type']):
            vals = np.array([template_type_dict[i][gal_typ] for i in range(len(config['sed_type']))])
            intp_vals = np.interp(ind_int, index, vals)
            for i in ind_int:
                template_type_dict_interp_[i][gal_typ] = intp_vals[i]

        #save as original template_type_dict
        template_type_dict = copy.copy(template_type_dict_interp_)
        del template_type_dict_interp_

        #Iterpolate the fluxes, between the templates for each filter
        f_mod_iterp = np.zeros((len(z_bins), len(ind_int), len(filters)), dtype=float)

        for i in range(len(z_bins)):
            for j, fltr in enumerate(filters):
                f_mod_iterp[i, :, j] = np.interp(ind_int, index, f_mod[i, :, j])

        #save as original f_mod
        f_mod = copy.copy(f_mod_iterp)
        del f_mod_iterp

    if key_not_none(config, 'save_template_fluxes'):
        import cPickle as pickle
        pickle.dump(
                    {'fluxes': f_mod, 'template_type': template_type_dict,
                    'filter_order': filters,
                    'filters_dict': config['filters'], 'z_bins': z_bins
                    }, open('template_fluxes.p', 'w')
                    )
        print ("template fluxes written to template_fluxes.p")
    #fast access to prior dictionary
    gal_mag_type_prior = GALPRIOR.prepare_prior_dictionary_types(template_type_dict)
    mags_bins = np.array(gal_mag_type_prior.keys(), dtype=float)

    #now load each file in turn.
    for fil in files:

        #get corresponding magnitudes
        orig_table = pyfits.open(fil)[1].data
        orig_cols = orig_table.columns
        n_gals = len(orig_table[orig_cols.names[0]])
        prior_mag = np.array(orig_table[config['PRIOR_MAGNITUDE']]*100, dtype=int)
        ID = np.array(orig_table[config['ID']])

        ADDITIONAL_OUTPUT_COLUMNS = []

        if key_not_none(config, 'ADDITIONAL_OUTPUT_COLUMNS'):

            #warn if all other requested columns are not in table
            for i, cl in enumerate(config['ADDITIONAL_OUTPUT_COLUMNS']):
                if cl not in orig_cols.names:
                    print ('Warning {:} not found in {:}. Continuing'.format(cl, fil))
                else:
                    ADDITIONAL_OUTPUT_COLUMNS.append(cl)

        f_obs = np.array([np.array(orig_table[i]) for i in MAG_OR_FLUX]).T
        ef_obs = np.array([np.array(orig_table[i]) for i in MAG_OR_FLUXERR]).T

        ind_process = np.ones(len(f_obs), dtype=bool)
        if config['INPUT_MAGS']:
            for i in range(len(MAG_OR_FLUX)):
                ind_process *= (f_obs[:, i] != config['mag_unobs'])

            #we are dealing with MAGS
            f_obs = np.power(10, -0.4 * f_obs)
            ef_obs = (np.power(10, (0.4 * ef_obs)) - 1) * f_obs

        else:
            #we are dealing with FLUXEs
            for i in np.arange(len(MAG_OR_FLUX)):
                ind_process *= (ef_obs[:, i] > 0)

        #add photoemtric offset error
        ef_obs = np.sqrt(ef_obs * ef_obs + np.power(zp_frac * f_obs, 2))

        #apply photo-z offsets
        f_obs = f_obs * zp_offsets
        ef_obs = ef_obs * zp_offsets

        #get normalised flux column
        norm_col = config['filters'][config['normalisation_filter']]['MAG_OR_FLUX']
        ind_norm_flux = np.where([i == norm_col for i in MAG_OR_FLUX])[0][0]

        if ind_norm != ind_norm_flux:
            print ("problem the template and real fluxes are out of order!: != ", ind_norm, ind_norm_flux)
            sys.exit()

        if key_not_none(config, 'output_pdfs'):
            pdf_file = fil.split('.fits')[0] + '.pdf.h5'
            df = pd.DataFrame()
            df.to_hdf(pdf_file, 'pdf_predictions', append=True)
            df.to_hdf(pdf_file, 'point_predictions', append=True)
            df.to_hdf(pdf_file, 'info', append=True)
            df2 = pd.DataFrame({'z_bin_centers': z_bins})
            df2.to_hdf(pdf_file, key='info', append=True)
            store = pd.HDFStore(pdf_file)

        nf = len(filters)

        #prepare for trivial parralisation using job_lib see  Parrallelise
        #above for an example. Split into 50k chunks
        ind = np.arange(n_gals)[ind_process]

        parr_lsts = []
        if key_not_none(config, 'n_jobs'):
            parr_lsts = []
            ind_ = np.array_split(ind, int(len(ind) / 50000) + 2)
            for ind1 in ind_:
                parr_lsts.append([ind1, f_obs[ind1], ef_obs[ind1], prior_mag[ind1], f_mod, gal_mag_type_prior, z_bins, config])

            res1 = Parrallelise(n_jobs=config['n_jobs'], method=parr_loop, loop=parr_lsts).run()
        else:
            #we do not want to parralise. let's send all the data required to the same function in one go.
            parr_lsts = [ind, f_obs, ef_obs, prior_mag, f_mod, gal_mag_type_prior, z_bins, config]

            #this must be kept as a list, so that we can loop over it, as it it was parrellised
            res1 = [parr_loop(parr_lsts)]

        #free space
        del parr_lsts

        #results files
        z_max_post = np.zeros(n_gals) + np.nan
        mean = np.zeros(n_gals) + np.nan
        sigma = np.zeros(n_gals) + np.nan
        median = np.zeros(n_gals) + np.nan
        mc = np.zeros(n_gals) + np.nan
        sig68 = np.zeros(n_gals) + np.nan
        KL_post_prior = np.zeros(n_gals) + np.nan
        min_chi2 = np.zeros(n_gals) + np.nan
        template_type = np.zeros(n_gals, dtype=float) + np.nan
        if key_not_none(config, 'output_pdfs'):
            pdfs_ = np.zeros((n_gals, len(z_bins))) + np.nan

        #let's combine all the results from the parrallel (or not) jobs
        for res in res1:
            ind_ = res['ind']
            z_max_post[ind_] = res['z_max_post']
            mean[ind_] = res['mean']
            sigma[ind_] = res['sigma']
            median[ind_] = res['median']
            mc[ind_] = res['mc']
            sig68[ind_] = res['sig68']
            KL_post_prior[ind_] = res['KL_post_prior']
            min_chi2[ind_] = res['min_chi2']
            template_type[ind_] = sed_float_list[res['maxL_template_ind']]
            if key_not_none(config, 'output_pdfs'):
                pdfs_[ind_] = res['pdfs_']

        #free up space
        del res1
        del res

        #saving point predictions as .fits files
        cols = {'MEAN_Z': mean, 'Z_SIGMA': sigma, 'MEDIAN_Z': median,
                'Z_MC': mc, 'Z_SIGMA68': sig68, 'KL_POST_PRIOR': KL_post_prior,
                'TEMPLATE_TYPE': template_type, 'MINCHI2': min_chi2}

        new_cols = pyfits.ColDefs([pyfits.Column(name=col_name, array=cols[col_name], format='D') for col_name in cols.keys()])

        id_cols = pyfits.ColDefs([pyfits.Column(name='ID', array=ID, format='K')])

        if len(ADDITIONAL_OUTPUT_COLUMNS) > 0:
            add_cols = pyfits.ColDefs([pyfits.Column(name=col_name, array=orig_table[col_name], format='D') for col_name in ADDITIONAL_OUTPUT_COLUMNS])
            hdu = pyfits.BinTableHDU.from_columns(id_cols + add_cols + new_cols)
        else:
            hdu = pyfits.BinTableHDU.from_columns(id_cols + new_cols)

        fname = fil.replace('.fits', '.BPZ' + output_file_suffix + '.fits')
        hdu.writeto(fname)

        #free memory
        del new_cols, add_cols

        #save pdf files
        if key_not_none(config, 'output_pdfs'):

            if len(ADDITIONAL_OUTPUT_COLUMNS) > 0:
                for col_name in ADDITIONAL_OUTPUT_COLUMNS:
                    cols[col_name] = np.array(orig_table[col_name])

            #split into manageable write chunks to save RAM [otherwise blows up with >1M rows!]
            inds = np.array_split(np.arange(n_gals), int(n_gals/200000) + 2)
            for ind in inds:
                cols_ = {'ID': ID[ind]}
                for j in cols.keys():
                    cols_[j] = cols[j][ind]

                df2 = pd.DataFrame(cols_)
                df2.to_hdf(pdf_file, key='point_predictions', format='table', append=True, complevel=5, complib='blosc')

                #free memory
                del cols_
                del df2
                if verbose:
                    print 'entering pdf'
                post_dict = {'KL_POST_PRIOR': KL_post_prior[ind], 'MEAN_Z': mean[ind], 'ID': ID[ind]}

                for ii in np.arange(len(z_bins)):
                    post_dict['pdf_{:0.4}'.format(z_bins[ii] + (z_bins[0]+z_bins[1])/2.0)] = pdfs_[ind, ii]
                if verbose:
                    print 'generating DataFrame'

                df2 = pd.DataFrame(post_dict)

                if verbose:
                    print 'writing pdf'

                df2.to_hdf(pdf_file, key='pdf_predictions', format='table', append=True, complevel=5, complib='blosc')

                #free memory
                del df2
                del post_dict
                if verbose:
                    print 'leaving pdf'
            del inds
            del cols
            
if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 2 or 'help' in args or '-h' in args:
        _help()

    #args = ['bpzConfig.yaml', 'WL_CLASS.METACAL.rescaled.slr.cosmos.v2._96_200_sampled.fits']
    main(args)
#    """