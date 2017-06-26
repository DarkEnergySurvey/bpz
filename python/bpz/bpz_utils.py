"""
A series of utility function and classes needed to run BPZ. These used
to live on the main executable code, but I've moved them here for
simplicity.

FM

"""

import os
import sys
from joblib import Parallel, delayed
import numpy as np

# Setting paths
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


def cmdline():

    import argparse
    import yaml

    # 1. We make a proto-parse use to read in the default yaml
    # configuration file, Turn off help, so we print all options in response to -h
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("-c", "--config",help="Specify config file")
    args, remaining_argv = conf_parser.parse_known_args()
    conf_defaults = yaml.load(open(args.config))

    # 2. This is the main parser
    parser = argparse.ArgumentParser(description="Run BPZv1 over an input catalog",
                                     # Inherit options from config_parser
                                     parents=[conf_parser])
    parser.add_argument("--files", action="store",nargs='+',default=None,required=True,
                        help="Name of input fits catalog(s)")
    parser.add_argument("--outbpz", action="store",default=None,required=True,
                        help="Name of output bpz fits catalog")
    # The positional arguments
    parser.add_argument("--prior_name", action="store",default=None,
                        help="prior name, any set you like (i.e.: bpz.sed_prior_file.des_y1_prior). See sed_proir_file.py for details")
    parser.add_argument("--n_jobs", action="store", default=1, type=int,
                        help="Number of jobs/cpu per run")
    parser.add_argument("--gal_chunk_size:", action="store", default=0, type=int,
                        help="Number of galaxies per loop (0=auto)")
    parser.add_argument("--AB_DIR", action="store",default=AB_DIR,
                        help="Location of AB files")
    parser.add_argument("--ID", action="store",default='NUMBER',
                        help="ID column to use from input catalog")
    
    # Set the defaults of argparse using the values in the yaml config file
    parser.set_defaults(**conf_defaults)
    args = parser.parse_args(args=remaining_argv)

    
    # Update keys in case these were undefined
    #if not args.BPZ_BASE_DIR:
    #    try:
    #        args.BPZ_BASE_DIR = os.environ['BPZ_DIR']
    #    except:
    #        print "Cannot find BPZ_DIR on environment, will guess the path"
    #        args.BPZ_BASE_DIR = os.path.dirname(os.path.realpath(__file__)).split('python')[0]
    #    print "# Will set BPZ_BASE_DIR to %s" % args.BPZ_BASE_DIR

    print vars(args)
    print "----"
    for k, v in vars(args).iteritems():
        print k, v
    return args




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

def e_mag2frac(errmag):
    """Convert mag error to fractionary flux error"""
    return 10. ** (.4 * errmag) - 1.

def help():
    print "new version BPZ"
    print "call like"
    print "bpzv1.py pathToConfigFile.yaml pathTofits*.fits"
    writeExampleConfig()
    sys.exit()

# Format time
def elapsed_time(t1,verb=False):
    import time
    t2    = time.time()
    stime = "%dm %2.2fs" % ( int( (t2-t1)/60.), (t2-t1) - 60*int((t2-t1)/60.))
    if verb:
        print >>sys.stderr,"Elapsed time: %s" % stime
    return stime

def writeExampleConfig():
    import textwrap

    example_yaml = os.path.join(os.getcwd(),'exampleBPZConfig.yaml')
    if os.path.isfile(example_yaml) is False:

        try:
            BPZ_PATH = os.environ['BPZ_DIR']
        except:
            print "Cannot find BPZ_DIR on environment, will guess the path"
            BPZ_PATH = os.path.dirname(os.path.realpath(__file__)).split('python')[0]
        f = open(example_yaml, 'w')
        txt = textwrap.dedent(
"""
#redshift bins min, max, width
redshift_bins: [0.01, 3.5, 0.01]

#either set here, or we will determine this from the code.
#Default: determine from code
BPZ_BASE_DIR:
#for SV this is %s/etc/AB_BPZ_ORIG
#for Y1 this is %s/etc/AB_BPZ_HIZ
AB_DIR: %s/etc/AB_BPZ_HIZ

#--------for Y1 this is----------
#spectra list. This *must* match the sed_type below. They must exist as expected in AB_DIR/*.AB
sed_list: [El_B2004a.sed, Sbc_B2004a.sed, Scd_B2004a.sed,Im_B2004a.sed, SB3_B2004a.sed, SB2_B2004a.sed]

#Either E/S0 Spiral or Irr (elliptical/Spherical, spiral, Irregular). The SEDs will be interpolated in the order of the list. They *should* be interpolated as E/S0->Spiral->Irr
sed_type: [E/S0, Spiral, Spiral, Irr, Irr, Irr]

#--------for SV this is----------
#spectra list. This *must* match the sed_type below. They must exist as expected in AB_DIR/*.AB
#sed_list: [El_B2004a.sed, Sbc_B2004a.sed, Scd_B2004a.sed, Im_B2004a.sed, SB3_B2004a.sed, SB2_B2004a.sed, ssp_25Myr_z008.sed,ssp_5Myr_z008.sed]

#Either E/S0 Spiral or Irr (elliptical, spiral, Irregular). The SEDs will be interpolated in the order of the list. They *should* be interpolated as E/S0->Spiral->Irr
#sed_type: [E/S0, Spiral, Spiral, Irr, Irr, Irr, Irr, Irr]


#go crazy and reorder all spectra types? Note: this is properly unphysical,
#due to interpolation reordering above!
rearrange_spectra: False

# prior name, any set you like. See sed_proir_file.py for details.
#for SV use sed_prior_file.des_sva1_prior
prior_name: sed_prior_file.des_y1_prior

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
output_pdfs:

#N_INTERPOLATE_TEMPLATES: Blank means No
INTERP: 8

#Should we output the templates as a dictionary:
#if yes, provide a pickle file path.
#if this file aleady exists, the code will stop.
output_sed_lookup_file:

SED_DIR: %s/etc/SED

#should we parralise the loops?
n_jobs: 5

#print some information to screen
verbose: True
""" % tuple([BPZ_PATH]*4))
        f.write(txt)
        print "An example file exampleBPZConfig.yaml has been written to disk"
    else:
        print "Example file already exists"
    return
    
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

    #import bpz.bpz_utils as bpz_utils
    import bpz.bh_photo_z_validation as pval
    from scipy.stats import entropy
    import time
    
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

    verbose = False
    if key_not_none(config, 'verbose'):
        verbose = config['verbose']

    if verbose:
        #print 'loop %s complete' % loop_number
        print "# Total loop %s time: %s" % (loop_number, elapsed_time(t0))

    return {'ind': ind_, 'mean': mean, 'sigma': sigma, 'median': median, 'sig68': sig68, 'mode': mode, 'z_minchi2': z_minchi2,
            'KL_post_prior': KL_post_prior, 'pdfs_': pdfs_, 'mc': mc,
            'min_chi2': min_chi2_arr, 'max_z_marg_likelihood': max_z_marg_likelihood, 'maxL_template_ind': maxL_template_ind}


