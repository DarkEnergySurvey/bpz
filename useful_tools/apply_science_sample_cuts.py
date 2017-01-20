#! /usr/bin/env python
"""This helper function should be run, once the redshift predictions have been made. It will massages data according to each specified science working group
Author: Ben Hoyle

Call like:

%>apply_science_sample_cuts.py YL|LSS|WL PathToPredictions.*.fits

This will impose the post processing selection

Note, these columns are required to be in the files"
COADD+MOF: MEAN_Z Z_SIGMA_68
depending on science sample + photometry
for COADD: MAG_AUTO_G MAG_AUTO_R MAG_AUTO_I MAG_AUTO_Z WEIGHT M1 M2 MAGERR_AUTO_I
for MOF: MAG_MOF_G MAG_MOF_R MAG_MOF_I MAG_MOF_Z mcal_e1_1p mcal_e1_1m mcal_e2_2p mcal_e2_2m WEIGHT
"""

import sys
import numpy as np

# we knew it from the start


def help():
    print "%>apply_science_sample_cuts.py YL|LSS|WL PathToPredictions.*.fits"
    print "This will impose the post processing selection after the photo-z's are inplace"
    print "Note, these columns are required to be in the files"
    print "COADD+MOF: MEAN_Z Z_SIGMA_68"
    print "for COADD: MAG_AUTO_G MAG_AUTO_R MAG_AUTO_I MAG_AUTO_Z WEIGHT M1 M2 MAGERR_AUTO_I"
    print "for MOF: MAG_MOF_G MAG_MOF_R MAG_MOF_I MAG_MOF_Z mcal_e1_1p mcal_e1_1m mcal_e2_2p mcal_e2_2m WEIGHT"
    sys.exit()


def remove_crazy_colors(d, PHOT_, indicies=False):
    """remove crazy colors in photometry """

    if PHOT_ == 'MOF':
        cols = [['MAG_MOF_G', 'MAG_MOF_R'], ['MAG_MOF_R', 'MAG_MOF_I'], ['MAG_MOF_I', 'MAG_MOF_Z']]

    if PHOT_ == 'COADD':
        cols = [['MAG_AUTO_G', 'MAG_AUTO_R'], ['MAG_AUTO_R', 'MAG_AUTO_I'],
                ['MAG_AUTO_I', 'MAG_AUTO_Z']]

    ind = np.ones(len(d), dtype=bool)
    for m1, m2 in cols:
        ind *= (np.array(d[m1] - d[m2]) > -1) * (np.array(d[m1] - d[m2]) < 4)

    if indicies:
        return ind
    return d[ind]


def lensing_weight_cosmic_shear(obj):
    # return lensing weight (shape noise only) for an object from the matched shear catalogs

    if('m1' in obj.columns):
        # it's im3shape
        return obj['weight'] * (1. + (obj['m1'] + obj['m2']) / 2.)

    if('mcal_e1_1m' in obj.columns):
        # it's metacal
        return 0.5 * obj['weight'] * [(obj['mcal_e1_1p'] - obj['mcal_e1_1m']) / 0.02 + (obj['mcal_e2_2p'] - obj['mcal_e2_2m']) / 0.02]
    print("this is not a proper shape object. i am failing you")
    sys.exit()


def lensing_weight_stacked_shear(obj, zlens):
    """Calculate the lengin weight for different zlens"""
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    
    Ds = cosmo.angular_diameter_distance(obj['MEAN_Z'])
    Dds = cosmo.angular_diameter_distance_z1z2(zlens, obj['MEAN_Z'])
    
    return lensing_weight_cosmic_shear(obj) * Dds / Ds


def apply_cuts(d, sample):

    """d = astropy data table
    sample = science sample WL | LSS | GOLD
    """

    cols = {}
    if sample == 'LSS':
        #LSS sample defintion
        cols['IN_LSS_SAMPLE'] = np.array((d['MEAN_Z'] < 1.0) * (d['MEAN_Z'] > 0.6))

    elif sample == 'WL':
        #WL weights
        for z_l in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
            cols['WEIGHT_ZLENS_{:}'.format(z_l)] = lensing_weight_stacked_shear(obj, z_l)

    elif sample == 'Y1':
        #Y1 sample definition
        cols['IN_Y1_SAMPLE'] = np.array((d['MEAN_Z'] > 0) * (d['SIGMA_68'] < 1.0))

        if ('MAG_AUTO_I' in obj.columns):
            cols['IN_Y1_SAMPLE'] *= np.array((d['MAG_AUTO_I'] < 23.5) * (d['MAGERR_AUTO_I'] < 0.3))
            cols['IN_Y1_SAMPLE'] *= remove_crazy_colors(d, 'COADD', indicies=True)
            cols['IN_Y1_SAMPLE'] *= np.array((d['MAG_AUTO_I'] > 16))
        else:
            cols['IN_Y1_SAMPLE'] *= remove_crazy_colors(d, 'MOF', indicies=True)
            cols['IN_Y1_SAMPLE'] *= np.array((d['MAG_MOF_I'] < 23.5) * (d['MAG_MOF_I'] > 16))

    return cols

if __name__ == "__main__":

    from astropy.io import fits as pyfits
    import os

    args = sys.argv[1:]

    if (len(args) < 2):
        help()

    sample = args[0]
    files = args[1:]

    #loop over all files
    for i in files:
        print (i)
        orig_table = pyfits.open(i)[1].data
        orig_cols = orig_table.columns
        cols = apply_cuts(orig_table, sample)
        new_cols = pyfits.ColDefs([pyfits.Column(name=col_name, array=cols[col_name], format='D') for col_name in cols.keys()])

        hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)

        fname = i.replace('.fits', sample + '.fits')
        hdu.writeto(fname)
        print ("Science sample {:} identified in {:}".format(sample, i))
        print ("Saving to {:} with new column names {:}".format(i, cols.keys()))
        import os
        os.rename(fname, i)