#! /usr/bin/env python

"""BPZ in 2017

author: benhoyle1212@gmail.com

about: This codes implements some of the key features in the original Bayesean Probability Redshift BPZ (Benitez 200X), with a more transparent user interface, and the native .fits support. 

    The original BPZ is still very useful. e.g determining photometric offsets is not implemented

To do:
Deal with missing / unseen data properly

"""

import sys

def help():
    print "new version BPZ"
    print "call like"
    print "bpzv1.py pathToConfigFile.yaml pathTofits*.fits"
    sys.exit()


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

@profile
def main(args):
    import yaml
    import numpy as np
    from galaxy_type_prior import GALAXYTYPE_PRIOR
    import os
    import copy
    from astropy.io import fits as pyfits
    import random as rdm
    import copy
    from joblib import Parallel, delayed
    import bh_photo_z_validation as pval
    from fastloop import bpz_loop

    config = yaml.load(open(args[0]))

    files = args[1:]
    if isinstance(files, list) is False:
        files = [files]

    if key_not_none(config, 'BPZ_BASE_DIR') is False:
        config['BPZ_BASE_DIR'] = '/' + '/'.join([i for i in os.path.realpath(__file__).split('/')[0:-1]]) + '/'

    if key_not_none(config, 'AB_DIR') is False:
        config['AB_DIR'] = config['BPZ_BASE_DIR'] + '../../templates/ABcorr/'

    if key_not_none(config, 'FILTER_DIR') is False:
        config['FILTER_DIR'] = config['BPZ_BASE_DIR'] + 'FILTER/'

    if key_not_none(config, 'SED_DIR') is False:
        config['SED_DIR'] = config['BPZ_BASE_DIR'] + 'SED/'

    output_file_suffix = ' '
    if key_not_none(config, 'SED_DIR'):
        output_file_suffix = config['output_file_suffix']

    #load redshift bins
    z_bins = np.arange(config['redshift_bins'][0], config['redshift_bins'][1] + config['redshift_bins'][2], config['redshift_bins'][2])

    #load filters. Direct copy past from bpz.py
    filters = config['filters'].keys()

    MAG_OR_FLUX = [config['filters'][i]['MAG_OR_FLUX'] for i in filters]
    MAG_OR_FLUXERR = [config['filters'][i]['ERR'] for i in filters]

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
        #print 'num_interps', num_interps
        #1 /0
        #generate some dummy indicies that are intergers spaced between 0 -- num_interps
        index = np.linspace(0, num_interps, len(config['sed_list']), endpoint=True, dtype=int)
        
        #generate values between dummy indicies. These will be interpolated 
        ind_int = np.arange(num_interps + 1)

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


    # some small value to truncate probs.
    eps = 1e-300
    eeps = np.log(eps)

    #now load each file in turn.
    for fil in files:

        #get corresponding magnitudes
        orig_table = pyfits.open(fil)[1].data
        orig_cols = orig_table.columns
        n_gals = len(orig_table[orig_cols.names[0]])
        prior_mag = orig_table[config['PRIOR_MAGNITUDE']]

        #warn if all other requested columns are not in table
        for i in config['ADDITIONAL_OUTPUT_COLUMNS']:
            if i not in orig_cols.names:
                print ('Warning {:} not found in {:}. Continuing'.format(i, fil))

        f_obs = np.array([np.array(orig_table[i]) for i in MAG_OR_FLUX]).T
        ef_obs = np.array([np.array(orig_table[i]) for i in MAG_OR_FLUXERR]).T

        ind_process = np.ones(len(f_obs), dtype=bool)
        if config['MAGS']:
            for i in range(len(MAG_OR_FLUX)):
                ind_process *= (f_obs[:, i] != config['mag_unobs'])

            #we are dealing with MAGS
            f_obs = np.power(10, -0.4 * f_obs)
            ef_obs = (np.power(10, (0.4 * ef_obs)) - 1) * f_obs

        else:
            #we are dealing with FLUXEs
            for i in np.arange(len(MAG_OR_FLUX)):
                ind_process *= (f_obs[:, i] > 0)
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

        nf = len(filters)

        #results files
        #cols = {'MEAN_Z': mean, 'Z_SIGMA': sigma, 'MEDIAN_Z': median,
        #    'Z_MC': mc, 'Z_SIGMA68':sig68}
        i = 2000
        cols = bpz_loop(i, z_bins, f_obs, ef_obs, f_mod, prior_mag, GALPRIOR, template_type_dict)
        #  int n_gals, double [:] z_bins, double[:, :] f_obs,  double[:] ef_obs, double[:, :, :] f_mod , double[:] prior_mag, GALPRIOR, template_type_dict)

        delta_z = res['z_max_post'] - orig_table['REDSHIFT'][0: i]
        print 'mean, std, len', np.median(delta_z), np.mean(delta_z), np.std(delta_z), len(delta_z) 
        1 / 0

        if False:
            print config['ADDITIONAL_OUTPUT_COLUMNS']

            new_cols = pyfits.ColDefs([pyfits.Column(name=col_name, array=cols[col_name], format='D') for col_name in cols.keys()])
            add_cols = pyfits.ColDefs([pyfits.Column(name=col_name, array=orig_table[col_name], format='D') for col_name in config['ADDITIONAL_OUTPUT_COLUMNS']])

            hdu = pyfits.BinTableHDU.from_columns(add_cols + new_cols)

            fname = fil.replace('.fits', 'BPZ.'+output_file_suffix+'.fits')
            hdu.writeto(fname)

args = sys.argv[1:]

if len(args) < 2 or 'help' in args or '-h' in args:
    help()

if __name__ == '__main__':
    main(args)
    
