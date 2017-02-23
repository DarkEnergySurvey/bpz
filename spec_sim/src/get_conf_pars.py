from __future__ import print_function
import numpy as np
import astropy.io.fits as pyfits
from det_feature import *
import sys

"""
Some lines are blended - so can we take those into account???
"""

try:
    spec_file = sys.argv[1]
    feat_file = sys.argv[2]
except:
    print('use: python get_conf_pars.py <spec file> <out file>')
    sys.exit()



# read in large spec file from Chihway
spec_data = pyfits.open(spec_file)
n_spec = len(spec_data[1].data['ID'])

# read in list of possible features
feat_type = {'E' : lambda x,y: lines.append([float(x),y]),
             'A' : lambda x,y: lines.append([float(x),y]),
             'EA' : lambda x,y: lines.append([float(x),y]),
             'AE' : lambda x,y: lines.append([float(x),y]),
             'B' : lambda x,y: breaks.append(np.array(x.split(':'),dtype=float)),
             'X' : lambda x,y: None}
lines = []
breaks = []
f = open('features.dat','r')
for line in f:
    if line[0] != "#":
        feat_type[line.split()[2]](line.split()[0],line.split()[2])
print('lines:',lines)
print('breaks:',breaks)
#print(spec_data[2].data['WAVE'])
line_array = []
break_array = []
#for spec in range(50):
for spec in range(n_spec):
    #print(spec_data[1].data['ID'][spec])
    # for each spectrum we want to pass the spectrum, noise and list of features to the det_feature routine, which will return the S/N in each feature.
    wave = spec_data[2].data['WAVE']/(1.+spec_data[1].data['Z'][spec])
    # SN & spec need some regularisation
    tmp_spec = spec_data[1].data['SPEC_OBS'][spec]
    tmp_spec[np.isnan(tmp_spec)] = 0.
    tmp_SN = spec_data[1].data['SPEC_SN'][spec]
    tmp_SN[np.isnan(tmp_SN)] = 1.e-15
    tmp_SN[np.where((tmp_SN==0.))[0]] = 1.e-15
    # if spec_true is exactly zero (due to spec window for instance), then computed noise (from S/N) will be zero and we get errors. Obviously S/N in these cases should be zero (not the 1.e-15 above). So we should artificially set the noise to a very high value.
    tmp_noise = spec_data[1].data['SPEC_TRUE'][spec]/tmp_SN
    tmp_noise[spec_data[1].data['SPEC_TRUE'][spec]==0.] = 1.e15
    SN_lines, SN_breaks = assemble_SN(wave, tmp_spec, tmp_noise, lines, breaks)
    line_array.append(SN_lines)
    break_array.append(SN_breaks)
    #print(spec_data[1].data['Z'][spec], spec_data[1].data['TMAG_I'][spec], SN_lines, SN_breaks)
line_array = np.array(line_array)
break_array = np.array(break_array)
print(line_array.shape,break_array.shape)

prihdr = pyfits.Header()
prihdu = pyfits.PrimaryHDU(header=prihdr)

# might change to have each line and break use their name from the line list file
c1 = pyfits.Column(name='ID', format='K', array=spec_data[1].data['ID'])
c2 = pyfits.Column(name='Z', format='E', array=spec_data[1].data['Z'])
c3 = pyfits.Column(name='OMAG_G', format='E', array=spec_data[1].data['OMAG_G'])
c4 = pyfits.Column(name='OMAG_R', format='E', array=spec_data[1].data['OMAG_R'])
c5 = pyfits.Column(name='OMAG_I', format='E', array=spec_data[1].data['OMAG_I'])
c6 = pyfits.Column(name='OMAG_Z', format='E', array=spec_data[1].data['OMAG_Z'])
c7 = pyfits.Column(name='OMAG_Y', format='E', array=spec_data[1].data['OMAG_Y'])
c8 = pyfits.Column(name='line_SN', format=str(line_array.shape[1])+'E', array=line_array) 
c9 = pyfits.Column(name='break_SN', format=str(break_array.shape[1])+'E',array=break_array)

CC = [c1, c2, c3, c4, c5, c6, c7, c8, c9]
# this is deprecated and needs changing
#hdu1 = pyfits.new_table(CC, nrows=len(spec_data[1].data['ID']))
hdu1 = pyfits.BinTableHDU.from_columns(CC, nrows=len(spec_data[1].data['ID']))

hdulist = pyfits.HDUList([prihdu, hdu1])
hdulist.writeto(feat_file, clobber=True)
