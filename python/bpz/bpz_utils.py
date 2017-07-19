"""
A series of utility functions and classes needed to run BPZ. These used
to live on the main executable code, but I've moved them here for
simplicity and make them more easy to maintain.

FM, June 2017

"""

import os
import sys
from joblib import Parallel, delayed
import numpy as np
import time
import copy
import pandas as pd
import cPickle as pickle
import bpz.bh_photo_z_validation as pval
from bpz.galaxy_type_prior import GALAXYTYPE_PRIOR
import fitsio
import collections
import logging
import math

max_gal_chunk_size = 50000

# Setting paths
try:
    SED_DIR = os.environ['SED_DIR']
except:
    print "Cannot find SED_DIR on environment, will guess the path"
    SED_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)).split('python')[0],'etc/SED_DIR')

try:
    BPZ_PATH = os.environ['BPZ_DIR']
except:
    print "Cannot find BPZ_DIR on environment, will guess the path"
    BPZ_PATH = os.path.dirname(os.path.realpath(__file__)).split('python')[0]

try:
    AB_DIR = os.environ['AB_DIR']
except:
    print "Cannot find AB_DIR on environment, will guess the path"
    AB_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)).split('python')[0],'etc/AB_DIR')

def create_logger(level=logging.NOTSET,name='default'):
    logging.basicConfig(level=level,
                        format='[%(asctime)s] [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(name)
    return logger

# Create a logger for all functions
LOGGER = create_logger(level=logging.NOTSET,name='BPZ')

def cmdline():

    import argparse
    import yaml

    # 1. We make a proto-parse use to read in the default yaml
    # configuration file, Turn off help, so we print all options in response to -h
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("-c", "--config",help="BPZ yaml config file")
    args, remaining_argv = conf_parser.parse_known_args()
    conf_defaults = yaml.load(open(args.config))

    # 2. This is the main parser
    parser = argparse.ArgumentParser(description="Run BPZv1 over an input catalog",
                                     # Inherit options from config_parser
                                     parents=[conf_parser])
    parser.add_argument("--incat", action="store",default=None,required=True,
                        help="Name of input fits catalog")
    parser.add_argument("--outbpz", action="store",default=None,required=True,
                        help="Name of output bpz fits catalog")
    # The positional arguments
    parser.add_argument("--redshift_bins", action="store",nargs='+',default=None,type=float,
                        help="Redshift bins: min, max and width")
    parser.add_argument("--rearrange_spectra", action="store_true", default=False,
                        help="Go crazy and reorder all spectra types")
    parser.add_argument("--prior_name", action="store",default=None,
                        help="prior name, any set you like (i.e.: bpz.sed_prior_file.des_y1_prior). See sed_proir_file.py for details")
    parser.add_argument("--n_jobs", action="store", default=1, type=int,
                        help="Number of jobs/cpu per run")
    parser.add_argument("--gal_chunk_size", action="store", default=0, type=int,
                        help="Number of galaxies per loop (0=auto)")
    parser.add_argument("--mag_unobs", action="store", default=-99, type=float,
                        help="Objects not observed (default=-99)")
    parser.add_argument("--mag_undet", action="store", default=99, type=float,
                        help="Objects not detected (default=99)")
    parser.add_argument("--normalisation_filter", action="store", default='DECam_2014_i.res',
                        help="Filter name we use for flux normalisation")
    parser.add_argument("--INPUT_MAGS",action="store_true", default=False,
                        help="Work with MAGS [True] or FLUX [False]")
    parser.add_argument("--AB_DIR", action="store",default=AB_DIR,
                        help="Location of AB files")
    parser.add_argument("--PRIOR_MAGNITUDE", action="store",default=None,
                        help="Column Identifier for PRIOR_MAGNITUDE")
    parser.add_argument("--SED_DIR", action="store",default=SED_DIR,
                        help="Path to SED Directory")
    parser.add_argument("--ID", action="store",default='NUMBER',
                        help="ID column to use from input catalog")
    parser.add_argument("--INTERP", action="store",type=int, default=8,
                        help="Interpolate templates")
    parser.add_argument("--MIN_MAGERR", action="store",type=float, default=0.001,
                        help="The minimum magnitude error")
    parser.add_argument("--ADDITIONAL_OUTPUT_COLUMNS", action="store",nargs='+',default=None,
                        help="List of additional output columns")
    parser.add_argument("--output_pdfs", action="store",default=None,
                        help="Name of output pdf/h5 file to be produced (Default=None)")
    # We need to remobe this option
    parser.add_argument("--output_file_suffix", action="store",default='',
                        help="Want output a suffix for a filename")
    parser.add_argument("--output_sed_lookup_file", action="store",default=None,
                        help="Want output the templates as a dictionary")
    
    # Set the defaults of argparse using the values in the yaml config file
    parser.set_defaults(**conf_defaults)
    args = parser.parse_args(args=remaining_argv)
    
    # Update keys/options in case these were not properly undefined
    if not args.output_file_suffix:
        args.output_file_suffix = ''

    # Make sure default paths are not None
    if not args.AB_DIR:
        args.AB_DIR = AB_DIR

    if not args.SED_DIR:
        args.SED_DIR = SED_DIR

    #print "Will use:"
    #for k, v in vars(args).iteritems():
    #    print "%s: %s" % (k, v)

    return args, parser

def load_file(file_name):
    z, f = [], []
    for i in open(file_name, 'r'):
        tmp = i.replace('  ', ' ').replace('\n', '').split(' ')
        z.append(float(tmp[0]))
        f.append(float(tmp[1]))
    return np.array(z), np.array(f)

def e_mag2frac(errmag):
    """Convert mag error to fractionary flux error"""
    return 10. ** (.4 * errmag) - 1.

# Format time
def elapsed_time(t1,verb=False):
    import time
    t2    = time.time()
    stime = "%dm %2.2fs" % ( int( (t2-t1)/60.), (t2-t1) - 60*int((t2-t1)/60.))
    if verb:
        print >>sys.stderr,"Elapsed time: %s" % stime
    return stime

    
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
    files = glob.glob('*.csv')

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

def photoz_loop(lst):

    from scipy.stats import entropy
    
    """parralisation loop for joblib
    lst is a [] containing ind, f_obs_, ef_obs_, prior_mag_, f_mod, gal_mag_type_prior, z_bins, config 
    f_mod = tempalate fluxes [redshift, template, band]

    returns a dictionary with
    {'ind': ind_, 'mean': mean, 'sigma': sigma, 'median': median, 'sig68': sig68, 'z_max_post': z_max_post,'z_minchi2': z_minchi2,
            'KL_post_prior': KL_post_prior, 'pdfs_': pdfs_}
    for later combination
    """

    t0 = time.time()
    ind_, f_obs_, ef_obs_, prior_mag_, f_mod, gal_mag_type_prior, z_bins, config, loop_number = lst
    n_gals = len(ind_)
    #some small value to truncate probs.
    eps = 1e-300
    eeps = np.log(eps)

    #results arrays for this loop
    mean = np.zeros(n_gals) + np.nan
    z_minchi2 = np.zeros(n_gals) + np.nan
    sigma = np.zeros(n_gals) + np.nan
    median = np.zeros(n_gals) + np.nan
    max_z_marg_likelihood = np.zeros(n_gals) + np.nan
    mode = np.zeros(n_gals) + np.nan
    mc = np.zeros(n_gals) + np.nan
    sig68 = np.zeros(n_gals) + np.nan
    KL_post_prior = np.zeros(n_gals) + np.nan
    min_chi2_arr = np.zeros(n_gals) + np.nan
    maxL_template_ind = np.zeros(n_gals, dtype=int) - 1000
    pdfs_ = None
    if config['output_pdfs']:
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
        z_minchi2[i] = z_min_chi2

        likelihood = np.exp(-0.5 * np.clip(chi2 - min_chi2, 0., -2 * eeps))

        prior = np.zeros_like(likelihood)
        pr_mg = gal_mag_type_prior.keys()

        ind_mag_p = np.argmin(np.abs(prior_mag_[i] - np.array(pr_mg)))

        for j in np.arange(len(f_mod[0, :, 0])):
            prior[:, j] = gal_mag_type_prior[pr_mg[ind_mag_p]][j]

        #posterior is prior * Likelihood
        posterior = prior * likelihood
        #posterior = likelihood

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

        marg_likelihood = np.sum(likelihood, axis=1)
        marg_likelihood /= np.sum(marg_likelihood)
        max_z_marg_likelihood[i] = z_bins[np.where(marg_likelihood == np.amax(marg_likelihood))[0][0]]

        KL_post_prior[i] = entropy(marg_post, marg_prior)

        ind_max_marg = np.where(marg_post == np.amax(marg_post))[0][0]

        #define summary stats from the margenalised posterior.
        mean[i] = pval.get_mean(marg_post, z_bins)
        sigma[i] = pval.get_sig(marg_post, z_bins)
        median[i] = pval.get_median(marg_post, z_bins)
        mc[i] = pval.get_mc(marg_post, z_bins)
        sig68[i] = pval.get_sig68(marg_post, z_bins)
        mode[i] = z_bins[ind_max_marg]

        if config['output_pdfs']:
            pdfs_[i] = marg_post

    if config['verbose']:
        LOGGER.info("Total loop %s time: %s" % (loop_number, elapsed_time(t0)))

    return {'ind': ind_, 'mean': mean, 'sigma': sigma, 'median': median, 'sig68': sig68, 'mode': mode, 'z_minchi2': z_minchi2,
            'KL_post_prior': KL_post_prior, 'pdfs_': pdfs_, 'mc': mc,
            'min_chi2': min_chi2_arr, 'max_z_marg_likelihood': max_z_marg_likelihood, 'maxL_template_ind': maxL_template_ind}



def bpz_main():

    t0 = time.time()
    # Get the configuration and get it into a config dictionary 
    args,parser = cmdline()
    config = vars(args)
    
    # Load redshift bins
    z_bins = np.arange(config['redshift_bins'][0], config['redshift_bins'][1] + config['redshift_bins'][2], config['redshift_bins'][2])

    # Load filters.
    filters = config['filters'].keys()

    # Adding band 
    MAG_OR_FLUX = ["%s_%s" % (config['filters'][i]['MAG_OR_FLUX'],config['filters'][i]['band']) for i in filters]
    MAG_OR_FLUXERR = ["%s_%s" % (config['filters'][i]['ERR'],config['filters'][i]['band']) for i in filters]

    if config['INPUT_MAGS'] is None:
        config['INPUT_MAGS'] = ('mag' in MAG_OR_FLUX[0]) or ('MAG' in MAG_OR_FLUX[0])

    # Jazz with spectra_list. Make all combinations of everything
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

    # Each sed (in order) is a number between 1 -> Num seds.
    sed_float_list = np.arange(len(config['sed_list']), dtype=float)
    ind_sed_int = np.arange(len(config['sed_list']), dtype=int)

    # Identify the filter we will normalise to.
    ind_norm = np.where(np.array(filters) == config['normalisation_filter'])[0][0]

    # Get flux offsets, and errors.
    zp_offsets = np.array([config['filters'][i]['zp_offset'] for i in filters])
    zp_offsets = np.power(10, -.4 * zp_offsets)

    zp_error = np.array([config['filters'][i]['zp_error'] for i in filters])
    zp_frac = e_mag2frac(zp_error)
    
    if config['output_sed_lookup_file']:        
        if os.path.isfile(config['output_sed_lookup_file']):
            msg = config['output_sed_lookup_file'], ' already exists! remove this file before continuing'
            sys.exit("ERROR: %s" % msg)
        if config['SED_DIR'] is False:
            sys.exit("ERROR: define SED_DIR in input file")
        full_sed = {}
        for i, sed in enumerate(config['sed_list']):
            full_sed[i] = np.genfromtxt(config['SED_DIR'] + '/' + sed)
            print i,sed, full_sed[i]
    # Prepare the prior
    GALPRIOR = GALAXYTYPE_PRIOR(z=z_bins, tipo_prior=config['prior_name'],
                    mag_bins=np.arange(18, 24.1, 0.1),
                    template_type_list=config['sed_type'])

    """Populate the model fluxes here """
    # Prepare to populate the model fluxes
    f_mod = np.zeros((len(z_bins), len(config['sed_list']), len(filters)), dtype=float)

    # Which template types does each SED correspond to
    template_type_dict = {}
    for i, sed in enumerate(config['sed_list']):
        sed_ = sed.replace('.sed', '')

        # This dummy dictionary allows for new "galaxy SED types" to be defined in the future.
        dict_ = {}
        for gal_typ in np.unique(config['sed_type']):
            dict_[gal_typ] = 0.0
        dict_[config['sed_type'][i]] = 1.0
        template_type_dict[i] = copy.copy(dict_)

        for j, fltr in enumerate(filters):
            fltr_ = fltr.replace('.res', '')
            file_name = os.path.join(config['AB_DIR'], '.'.join([sed_, fltr_, 'AB']))
            z_, f_ = load_file(file_name)

            # Interpolate to the exact redshifts we care about
            f_mod[:, i, j] = np.interp(z_bins, z_, f_)

            # Clip, following original bpz
            f_mod[:, i, j] = np.clip(f_mod[:, i, j], 0, 1e300)

    # Interp between SEDs
    if config['INTERP']:

        # Wow many interpolated points?
        num_interps = (len(config['sed_list'])-1) * config['INTERP'] + len(config['sed_list'])

        # Generate some dummy indicies that are intergers spaced between 0 -- num_interps
        index = np.linspace(1, num_interps, len(config['sed_list']), endpoint=True, dtype=int) - 1

        # Generate values between dummy indicies. These will be interpolated
        ind_sed_int_orig = copy.copy(ind_sed_int)
        ind_sed_int = np.arange(num_interps)

        sed_float_list = np.interp(ind_sed_int, index, sed_float_list)

        # Should we interpolate between the SEDs also?
        if config['output_sed_lookup_file']:
            
            full_sed_ = copy.copy(full_sed)
            # For each interpolated index
            for i in ind_sed_int:
                # If at limits, no interpolation
                if i == 0:
                    full_sed_[i] = full_sed[0]
                elif i == np.amax(ind_sed_int):
                    full_sed_[i] = full_sed[ind_sed_int_orig[-1]]
                else:
                    # interpolate between SEDs
                    # get index of SED1 we will interpolate from
                    # we will interpolate to indexI +1  -> SED2
                    indexI = int(np.floor(sed_float_list[i]))

                    # get fraction of weight for this SED1
                    fractI = sed_float_list[i] - indexI

                    # identify overlapping Lambda ranges! these *maybe* different!
                    L1 = full_sed[ind_sed_int_orig[indexI]][:, 0]
                    L2 = full_sed[ind_sed_int_orig[indexI + 1]][:, 0]

                    lrange0in1 = (L1 >= np.amin(L2)) * (L1 <= np.amax(L2))
                    lrange1in0 = (L2 >= np.amin(L1)) * (L2 <= np.amax(L1))

                    # generate a common set of lamdas
                    delta_l = L2[lrange1in0][1] - L2[lrange1in0][0]

                    common_lrange = np.arange(np.amin(L2[lrange1in0]), np.amax(L2[lrange1in0]) + delta_l, delta_l)

                    # interpolate to this new grid
                    interp1 = np.interp(common_lrange, L1[lrange0in1], full_sed[ind_sed_int_orig[indexI]][lrange0in1, 1])
                    interp2 = np.interp(common_lrange, L2[lrange1in0], full_sed[ind_sed_int_orig[indexI + 1]][lrange1in0, 1])
                    # finally make a linear combination of SED1 + SED2
                    full_sed_[i] = [common_lrange, interp1 * fractI + (1.0 - fractI) * interp2]
            full_sed = copy.copy(full_sed_)
            del full_sed_

        # Frist interpolate the galaxy SED types. E.g.
        # {'E/S0': 1 'Spiral': 0 'Irr': 0} -> {'E/S0': 0.5 'Spiral': 0.5 'Irr': 0}
        # ->{'E/S0': 0 'Spiral': 1 'Irr': 0}
        template_type_dict_interp_ = {}
        for i in ind_sed_int:
            template_type_dict_interp_[i] = copy.copy(template_type_dict[0])

        for gal_typ in np.unique(config['sed_type']):
            vals = np.array([template_type_dict[i][gal_typ] for i in range(len(config['sed_type']))])
            intp_vals = np.interp(ind_sed_int, index, vals)
            for i in ind_sed_int:
                template_type_dict_interp_[i][gal_typ] = intp_vals[i]

        # Save as original template_type_dict
        template_type_dict = copy.copy(template_type_dict_interp_)
        del template_type_dict_interp_

        # Interpolate the fluxes, between the templates for each filter
        f_mod_iterp = np.zeros((len(z_bins), len(ind_sed_int), len(filters)), dtype=float)

        for i in range(len(z_bins)):
            for j, fltr in enumerate(filters):
                f_mod_iterp[i, :, j] = np.interp(ind_sed_int, index, f_mod[i, :, j])

        # save as original f_mod
        f_mod = copy.copy(f_mod_iterp)
        del f_mod_iterp

    if config['output_sed_lookup_file']:
        pickle.dump(
                    {'flux_per_z_template_band': f_mod, 'template_type': template_type_dict,
                    'filter_order': filters,
                    'filters_dict': config['filters'], 'z_bins': z_bins,
                    'SED': full_sed,
                    }, open(config['output_file_suffix'] + config['output_sed_lookup_file'], 'w')
                    )

        print ("template fluxes written to: ", config['output_file_suffix'] + config['output_sed_lookup_file'])
    #fast access to prior dictionary
    gal_mag_type_prior = GALPRIOR.prepare_prior_dictionary_types(template_type_dict)
    mags_bins = np.array(gal_mag_type_prior.keys(), dtype=float)


    # Load in the input catalog file get corresponding magnitudes
    if args.verbose:
        LOGGER.info("Reading input catalog: %s" % config['incat'])
    tab = fitsio.FITS(args.incat)
    orig_table = tab[1].read() # Change to 'OBJECTS' hdu!!!!
    orig_cols_names = orig_table.dtype.names
    n_gals = len(orig_table[orig_cols_names[0]])
    ID = np.array(orig_table[config['ID']])
    prior_mag = np.array(np.round(orig_table[config['PRIOR_MAGNITUDE']], 1) * 100).astype(np.int)

    ADDITIONAL_OUTPUT_COLUMNS = []
    if config['ADDITIONAL_OUTPUT_COLUMNS']:
        # warn if all other requested columns are not in table
        for i, cl in enumerate(config['ADDITIONAL_OUTPUT_COLUMNS']):
            if cl not in orig_cols_names:
                print ('Warning {:} not found in {:}. Continuing'.format(cl, config['incat']))
            else:
                ADDITIONAL_OUTPUT_COLUMNS.append(cl)

    f_obs = np.array([np.array(orig_table[i]) for i in MAG_OR_FLUX]).T
    ef_obs = np.array([np.array(orig_table[i]) for i in MAG_OR_FLUXERR]).T

    ind_process = np.ones(len(f_obs), dtype=bool)
    if config['INPUT_MAGS']:
        for i in range(len(MAG_OR_FLUX)):
            ind_process *= (f_obs[:, i] != config['mag_unobs'])

        # we are dealing with MAGS
        f_obs = np.power(10, -0.4 * f_obs)
        ef_obs = (np.power(10, (0.4 * ef_obs)) - 1) * f_obs

    else:
        # we are dealing with FLUXEs
        for i in np.arange(len(MAG_OR_FLUX)):
            ind_process *= (ef_obs[:, i] > 0)

    # add photoemtric offset error
    ef_obs = np.sqrt(ef_obs * ef_obs + np.power(zp_frac * f_obs, 2))

    # apply photo-z offsets
    f_obs = f_obs * zp_offsets
    ef_obs = ef_obs * zp_offsets

    # get normalised flux column
    norm_filter = config['normalisation_filter']
    norm_band = config['filters'][norm_filter]['band']
    norm_col = "%s_%s" % (config['filters'][norm_filter]['MAG_OR_FLUX'],norm_band)
    ind_norm_flux = np.where([i == norm_col for i in MAG_OR_FLUX])[0][0]

    if ind_norm != ind_norm_flux:
        sys.exit("ERROR:problem the template and real fluxes are out of order!: != ", ind_norm, ind_norm_flux)

    if config['output_pdfs']:
        pdf_file = config['output_pdfs']
        # Clobber previous file if exists
        if os.path.exists(pdf_file):
            if args.verbose: LOGGER.info("Removing pdfs file: %s" % pdf_file)
            os.remove(pdf_file)
        if args.verbose: LOGGER.info("Will write pdfs to file: %s" % pdf_file)
        df = pd.DataFrame()
        df.to_hdf(pdf_file, 'pdf_predictions', append=True)
        df.to_hdf(pdf_file, 'point_predictions', append=True)
        df.to_hdf(pdf_file, 'info', append=True)
        df2 = pd.DataFrame({'z_bin_centers': z_bins})
        df2.to_hdf(pdf_file, key='info', append=True)
        store = pd.HDFStore(pdf_file)

    # Nuber of filters
    nf = len(filters)

    # Prepare for trivial parralisation using job_lib see  Parrallelise above for an example. 
    ind = np.arange(n_gals)[ind_process]

    # Define chunk size
    if config['gal_chunk_size']:
        gal_chunk_size = config['gal_chunk_size']
        LOGGER.info("Will use input chunk_size=%s" % gal_chunk_size)
    else:
        gal_chunk_size = int(len(ind)/config['n_jobs'])
        LOGGER.info("Will use auto chunk_size=%s" % gal_chunk_size)

    # Resize gal_chunk_size if necessary to for to fit n_jobs
    if gal_chunk_size > max_gal_chunk_size:
        scale = math.ceil(gal_chunk_size/max_gal_chunk_size/config['n_jobs'])
        gal_chunk_size = int(gal_chunk_size/(config['n_jobs']*scale))
        LOGGER.info("Re-sizing chunk_size=%s" % gal_chunk_size)
    
    parr_lsts = []
    if config['n_jobs']:
        parr_lsts = []
        nloops = int(len(ind) / gal_chunk_size)
        if nloops < 2:  nloops = 2
        ind_ = np.array_split(ind, nloops)
        k = 1
        for ind1 in ind_:
            LOGGER.info("Preparing loop %s" % k)
            parr_lsts.append([ind1, f_obs[ind1], ef_obs[ind1], prior_mag[ind1], f_mod, gal_mag_type_prior, z_bins, config,k])
            k =  k + 1

        res1 = Parrallelise(n_jobs=config['n_jobs'], method=photoz_loop, loop=parr_lsts).run()
    else:
        # we do not want to parralise. let's send all the data required to the same function in one go.
        parr_lsts = [ind, f_obs, ef_obs, prior_mag, f_mod, gal_mag_type_prior, z_bins, config]

        # this must be kept as a list, so that we can loop over it, as it it was parrellised
        res1 = [photoz_loop(parr_lsts)]

    # free space
    del parr_lsts
    
    # results arrays for point predictions
    mode = np.zeros(n_gals) + np.nan
    z_minchi2 = np.zeros(n_gals) + np.nan
    z_max_marg_like = np.zeros(n_gals) + np.nan
    mean = np.zeros(n_gals) + np.nan
    sigma = np.zeros(n_gals) + np.nan
    median = np.zeros(n_gals) + np.nan
    mc = np.zeros(n_gals) + np.nan
    sig68 = np.zeros(n_gals) + np.nan
    KL_post_prior = np.zeros(n_gals) + np.nan
    min_chi2 = np.zeros(n_gals) + np.nan
    template_type = np.zeros(n_gals, dtype=float) + np.nan
    template_int = np.zeros(n_gals, dtype=int) - 999
    
    if config['output_pdfs']:
        pdfs_ = np.zeros((n_gals, len(z_bins))) + np.nan

    # let's combine all the results from the parrallel (or not) jobs
    for res in res1:
        ind_ = res['ind']
        mode[ind_] = res['mode']
        mean[ind_] = res['mean']
        sigma[ind_] = res['sigma']
        median[ind_] = res['median']
        mc[ind_] = res['mc']
        z_minchi2[ind_] = res['z_minchi2']
        z_max_marg_like[ind_] = res['max_z_marg_likelihood']
        sig68[ind_] = res['sig68']
        KL_post_prior[ind_] = res['KL_post_prior']
        min_chi2[ind_] = res['min_chi2']
        template_type[ind_] = sed_float_list[res['maxL_template_ind']]
        template_int[ind_] = res['maxL_template_ind']
        if config['output_pdfs']:
            pdfs_[ind_] = res['pdfs_']

    # free up space
    del res1
    del res

    # ----------------------------------------
    # Write out_bpz using fitsio instead
    # Built the array/dtype structure
    nrows = len(mean)
    cols = collections.OrderedDict()
    cols = {'MEAN_Z': mean,
            'Z_SIGMA': sigma,
            'MEDIAN_Z': median,
            'Z_MC': mc,
            'Z_SIGMA68': sig68,
            'KL_POST_PRIOR': KL_post_prior,
            'TEMPLATE_TYPE': template_type,
            'MINCHI2': min_chi2,
            'MODE_Z': mode,
            'Z_MINCHI2': z_minchi2,
            'Z_MAXMARG_LIKE': z_max_marg_like}
    dtypes = [(col_name,'f8') for col_name in cols.keys()]
    # Add 'str' type for these ones separately
    dtypes.insert(0,('ID','i8'))
    dtypes.insert(8,('TEMPLATE_ID','i8'))

    # Got Additional columns -- add them as f8
    if len(ADDITIONAL_OUTPUT_COLUMNS) > 0:
        for col_name in ADDITIONAL_OUTPUT_COLUMNS:
            # Figure out the dtype of the extra column and append
            col_dtype = ("%s" % orig_table[col_name].dtype)[1:]
            dtypes.append((col_name,col_dtype))
            cols[col_name] = orig_table[col_name]

    data_out = np.zeros(nrows, dtype=dtypes)
    data_out['ID'] = ID
    data_out['TEMPLATE_ID'] = template_int
    for col_name in cols.keys():
        data_out[col_name] = cols[col_name]

    fitsio.write(config['outbpz'], data_out, extname='OBJECTS', clobber=True)
    if args.verbose:
        LOGGER.info("Wrote output to file: %s" % config['outbpz'])
    # -- End of write  -----

    # free memory
    del data_out

    # Save pdf files
    if config['output_pdfs']:

        if len(ADDITIONAL_OUTPUT_COLUMNS) > 0:
            for col_name in ADDITIONAL_OUTPUT_COLUMNS:
                cols[col_name] = np.array(orig_table[col_name])

        # split into manageable write chunks to save RAM [otherwise blows up with >1M rows!]
        inds = np.array_split(np.arange(n_gals), int(n_gals/200000) + 2)
        for ind in inds:
            cols_ = {config['ID']: ID[ind]}
            for j in cols.keys():
                cols_[j] = cols[j][ind]

            df2 = pd.DataFrame(cols_)
            df2.to_hdf(pdf_file, key='point_predictions', format='table', append=True, complevel=5, complib='blosc')

            # free memory
            del cols_
            del df2
            if args.verbose:
                LOGGER.info('Entering pdf')
            post_dict = {'KL_POST_PRIOR': KL_post_prior[ind], 'MEAN_Z': mean[ind], config['ID']: ID[ind],
            'TEMPLATE_ID': template_int[ind]}

            for ii in np.arange(len(z_bins)):
                post_dict['pdf_{:0.4}'.format(z_bins[ii])] = pdfs_[ind, ii]
            if args.verbose: LOGGER.info('Generating DataFrame')
            df2 = pd.DataFrame(post_dict)
            if args.verbose: LOGGER.info('Writing pdf')

            df2.to_hdf(pdf_file, key='pdf_predictions', format='table', append=True, complevel=5, complib='blosc')

            # free memory
            del df2
            del post_dict
            if args.verbose:
                LOGGER.info('Leaving pdf')
        del inds
    # free space
    del cols

    # Close the h5 file store
    if config['output_pdfs']:
        store.close()

    if args.verbose:
        LOGGER.info("Total BPZ time %s" % elapsed_time(t0))

    # Done
    return 
