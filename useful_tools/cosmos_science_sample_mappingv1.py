#! /usr/bin/env python
import sys

def show_help(args):
    print ("""y1_cosmos_science_mappingv1.py CosDataFile.fits ScienceSamplesFile.fits mapping_file.yaml
        ScienceSamplesFile.fits must have columns LSS_CLASS WL_CLASS Y1_CLASS
        mapping_file.yaml should look like

FLUXES = ["FLUX_AUTO_G", "FLUX_AUTO_R", "FLUX_AUTO_I", "FLUX_AUTO_Z"]

flux_to_mag = {"FLUX_MOF_G" :"MAG_MOF_G", "FLUX_MOF_R" : "MAG_MOF_R", "FLUX_MOF_I": "MAG_MOF_I", 
              "FLUX_MOF_Z" : "MAG_MOF_Z"}

FLUXERRS = ["FLUXERR_MOF_G", "FLUXERR_MOF_R", "FLUXERR_MOF_I", "FLUXERR_MOF_Z"]

OTHER_COLS = ["CM_T"]
OTHER_ERRS = ["CM_T_ERR"]
ADDITIONAL_COLS = ["e1", "e2", 'mean_psf_e1_sky', 'mean_psf_e2_sky']

cosmo_redshifts: Cos.DES.Validation.resampled_redshifts_rescaled_pdfs.p

science_sample_file: Y1A1_MOF_SAMPLE_RAND_4pcnt.WL.LSS.Y1.CLASS.cut.fits

science_sample: LSS_CLASS

    where we will map in N-d using nd_abundance_matching.py on [vars, errors] and then determine new science scample like errors for errors.

        """)
    print(args)
    print sys.exit()

import numpy as np
args = sys.argv[1:]

if len(args) < 3:
    show_help(args)

from astropy.table import Table, vstack
import numpy as np
from math import log10
import cPickle as pickle
import yaml

# some useful functions for later
def get_zp(cat, cols, magcols):
# determine zeropoints based on what's in the catalog 
    zp=[]
    for t in zip(cols, magcols):
        if (t[1]==""): # it's not a flux/mag column at all
            zp.append(0.)
            continue
        for o in cat:
            if(o[t[0]]>0): # a proper flux measurement, use to determine ZP
                zp.append(o[t[1]]+2.5*log10(o[t[0]]))
                break # no need to go on down the catalog

    return zp


def error_mask(sgal, reference, matched_columns_error):
# select galaxies in reference for which all errors are smaller or equal than those of galaxy sgal
    #mask=(reference[matched_columns_error[0]]<=sgal[matched_columns_error[0]])
 
    #for p in matched_columns_error[1:]:
    #    mask=np.logical_and(mask,(reference[p]<=sgal[p]))
    arr1 = np.zeros(len(matched_columns_error))
    arr2 = np.zeros((len(reference), len(matched_columns_error)))
    for i, col in enumerate(matched_columns_error):
        arr1[i] = np.array(sgal[col])
        arr2[:, i] = np.array(reference[col])

    mask = error_maskv1(arr1, arr2)
    return mask

def error_maskv1(arr1, arr2, mtch=np.less_equal):
    #replicate arr1 to have shape arr2
    ind_le = mtch(arr2, np.tile(arr1, (len(arr2),1)))
    #find rows in which all conditions are true
    indarr = np.sum(ind_le, axis=1) == len(arr1)
    return indarr


def test_zp_equal(zp1, zp2, target):
# test whether these vectors of zeropoints are actually equal to the target value
    for i in zip(zp1, zp2):
        if(i[0] == 0 and i[1] == 0): # not a magnitude column
            continue
        if(abs(i[0] - target) > 0.001):
            return 0
        if(abs(i[1] - target) > 0.001):
            return 0
    return 1

def flux2mag(flux):
    #37.5 # default for negative fluxes
    mag = np.zeros(len(flux)) + 37.5
    ind = flux > 0.001
    mag[ind] = -2.5 * log10(flux[ind]) + 30
    return  mag

def fluxerr_2_magerr(flux, fluxerr):
    
    # default if negative or insignificant flux detected
    magerr = np.zeros(len(flux)) + 99
    ind = (flux > 0) * (flux > fluxerr)
    magerr[ind] = 1.0857*fluxerr[ind]/flux[ind] 
    return magerr



def test_zp_equal(zp1, zp2, target):
# test whether these vectors of zeropoints are actually equal to the target value
    for i in zip(zp1,zp2):
        if(i[0]==0 and i[1]==0): # not a magnitude column
            continue
        if(abs(i[0]-target)>0.001):
            return 0
        if(abs(i[1]-target)>0.001):
            return 0
    return 1

conf = yaml.load(open(args[2], 'r'))


eli_slr_table= path + "y3a1_fgcm_2_5_nside32_zpshift2_dered.fits"


matched_columns = conf['FLUXES'] + conf['OTHER_COLS']
# associated error columns
matched_columns_error = conf['FLUXERRS'] + conf['OTHER_ERRS']


reference = Table.read(conf['reference_catalog'])
selection = Table.read(conf['selection_catalog'])

#select 'good' object by flags here
selection = selection[selection[conf['SCIENCE_SAMPLE']] == 1]

#ToDo: only select columns of interest from these tables, to save RAM 
print(len(reference),"reference objects,",len(selection),"selected objects")

matched_magnitudes = []
for i in FLUXES:
    matched_magnitudes.append(flux_to_mag[i])
    
zp_reference=get_zp(reference, FLUXES, matched_magnitudes)
zp_selection=get_zp(selection, FLUXES, matched_magnitudes)

print("zeropoints:",zp_reference,zp_selection)

if(test_zp_equal(zp_reference,zp_selection,30)==0):
     