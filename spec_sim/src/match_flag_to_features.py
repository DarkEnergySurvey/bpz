from __future__ import print_function
import numpy as np
import astropy.io.fits as pyfits
import sys

"""
Marz file is in csv
Need to match line number with the table of simulated spectra that have been redshifted, to get ID number. This is in the simulated subset, *field*.fits.
Finally, match to the output from the feature detection.
"""

def match_ids_to_features(ids,data):
    lines = []
    breaks = []
    for id in ids:
        indx = np.where(data['ID'] == id)[0]
        if len(indx) > 1:
            print('warning: duplicate ID')
        indx = indx[0]
        lines.append(data['line_SN'][indx])
        breaks.append(data['break_SN'][indx])
    return lines, breaks

form = 0 # 0 for vipers, zcos; 1 for vvds  

try:
    marz_file = sys.argv[1]
    tmp = marz_file.split('_')
    s = '_'
    file_stem = s.join(tmp[0:2+form])
    field_file = file_stem+'_'+tmp[2+form]+'.fits'
    feature_file = file_stem+'_features.fits'
    out_file = field_file.split('.')[0]+'_feat.fits'
    field_num = tmp[2+form].split('field')[1]
    print(field_num)
except:
    print('use: python match_flag_to_features.py <marz file> <simulated subset> <feature file>')
    sys.exit()


# read marz file
AutoZ, FinZ, QOP = np.loadtxt(marz_file, delimiter=',', usecols=(8, 12, 13), unpack=True)

# read field file
subset_data = pyfits.open(field_file)[4].data

# read feature file
feature_data = pyfits.open(feature_file)[1].data

# get id list
ids = subset_data['ID']

# match ids to get features
line_sn, break_sn = match_ids_to_features(ids, feature_data)

# write out
line_array = np.array(line_sn)
break_array = np.array(break_sn)

### need to add all the marz info

prihdr = pyfits.Header()
prihdu = pyfits.PrimaryHDU(header=prihdr)

# might change to have each line and break use their name from the line list file
c1 = pyfits.Column(name='simID', format='K', array=np.array(ids))
c14 = pyfits.Column(name='field', format='K', array=np.repeat(field_num, len(ids)))
c2 = pyfits.Column(name='Z', format='E', array=subset_data['Z'])
c3 = pyfits.Column(name='OMAG_G', format='E', array=subset_data['OMAG_G'])
c4 = pyfits.Column(name='OMAG_R', format='E', array=subset_data['OMAG_R'])
c5 = pyfits.Column(name='OMAG_I', format='E', array=subset_data['OMAG_I'])
c6 = pyfits.Column(name='OMAG_Z', format='E', array=subset_data['OMAG_Z'])
c7 = pyfits.Column(name='OMAG_Y', format='E', array=subset_data['OMAG_Y'])
c8 = pyfits.Column(name='line_SN', format=str(line_array.shape[1])+'E', array=line_array) 
#c9 = pyfits.Column(name='break_SN', format=str(break_array.shape[1])+'E',array=break_array)
c9 = pyfits.Column(name='break_SN', format='E',array=break_array)
c10 = pyfits.Column(name='AutoZ', format='E', array=AutoZ)
c11 = pyfits.Column(name='FinZ', format='E', array=FinZ)
c12 = pyfits.Column(name='QOP', format='E', array=QOP)
c13 = pyfits.Column(name='Observer', format='3A', array=np.repeat(tmp[3+form], len(ids)))
c15 = pyfits.Column(name='line_num', format='K', array=np.arange(len(ids))+1)
CC = [c1, c14, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c15]
hdu1 = pyfits.BinTableHDU.from_columns(CC, nrows=len(ids))

hdulist = pyfits.HDUList([prihdu, hdu1])
hdulist.writeto(out_file, clobber=True)
