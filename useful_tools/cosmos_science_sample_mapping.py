#! /usr/bin/env python
import sys

def show_help(args):
    print ("""y1_cosmos_science_mapping.py CosDataFile.fits ScienceSamplesFile.fits mapping_file.yaml
        ScienceSamplesFile.fits must have columns LSS_CLASS WL_CLASS Y1_CLASS
        mapping_file.yaml should look like

vars: [FLUX_MOF_G ,FLUX_MOF_R ,FLUX_MOF_I ,FLUX_MOF_Z ,PSF_FLUX_I ,PSF_FLUX_Z ,CM_T ,CM_FLUX_G ,CM_FLUX_R ,CM_FLUX_I ,CM_FLUX_Z ,CM_G_2 ,CM_G_COV_1 ,CM_G_COV_2 ,CM_G_COV_3 ,CM_G_COV_4 ,CM_FRACDEV]

errors: [FLUXERR_MOF_G ,FLUXERR_MOF_R ,FLUXERR_MOF_I ,FLUXERR_MOF_Z ,PSF_FLUX_ERR_G ,PSF_FLUX_ERR_R ,PSF_FLUX_ERR_I ,PSF_FLUX_ERR_Z ,CM_T_ERR ,CM_FLUX_COV_G_G ,CM_FLUX_COV_G_R ,CM_FLUX_COV_G_I ,CM_FLUX_COV_G_Z ,CM_FLUX_COV_R_G ,CM_FLUX_COV_R_R ,CM_FLUX_COV_R_I ,CM_FLUX_COV_R_Z ,CM_FLUX_COV_I_G ,CM_FLUX_COV_I_R ,CM_FLUX_COV_I_I ,CM_FLUX_COV_I_Z ,CM_FLUX_COV_Z_G ,CM_FLUX_COV_Z_R ,CM_FLUX_COV_Z_I ,CM_FLUX_COV_Z_Z ,CM_G_1 ,CM_FRACDEV_ERR]
    
    where we will map in N-d using nd_abundance_matching.py on [vars, errors] and then determine new science scample like errors for errors.

        """)
    print(args)
    print sys.exit()

import numpy as np
args = sys.argv[1:]

if len(args) < 3:
    show_help(args)

import yaml
import nd_abundance_matching as ndab
from astropy.table import Table

conf = yaml.load(open(args[2], 'r'))

d1 = Table.read(args[0])

kys = conf['vars'] + conf['errors']
data1 = np.zeros((len(d1), len(kys)), dtype=float)
for i, ky in enumerate(kys):
    data1[:, i] = np.array(d1[ky])

for sample in ['WL_CLASS', 'LSS_CLASS', 'Y1_CLASS']:

    d2 = Table.read(args[1])
    d2 = d2[d2[sample] == 1]
    d2 = d2[conf['vars'] + conf['errors']]

    data2 = np.zeros((len(d2), len(kys)), dtype=float)
    for i, ky in enumerate(kys):
        data2[:, i] = np.array(d2[ky])

    match_d1_d2 = ndab.match_data_clouds(data1, data2, num_percentiles=200)

    #add erorrs:
    err = sample.split('_')[0] + '_'
    for i, ky in enumerate(conf['errors']):
        d1[err + ky] = np.array([data2[ky][j] for j in match_d1_d2], dtype=float)

d1.write(args[0].replace('.fits') + '_resampleErrs.fits')    