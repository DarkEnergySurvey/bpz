from __future__ import print_function
import numpy as np
import calculate_weights
import astropy.io.fits as pyfits
import sys

"""
Spec file and target file will often be the same, but not always. Coded this way to allow flexibility. Hopefully it doesn't cause problems.
Typically:

python get_weights.py vvds_deep1_sim_flag.fits,1 vvds_deep1_sim_flag.fits,1
python get_weights.py vvds_deep1_sim_flag.fits,1 Buzzard_v1.1_truth_vvds_deep1.fit,1

Won't work quite as well as I hoped: magnitude and ID arrays are slightly different between buzzard file and spec _feat.fits files. Could use try, or auto detect based on file name????

"""
try:
    spec_file = sys.argv[1].split(',')[0]
    spec_ext = int(sys.argv[1].split(',')[1])
    target_file = sys.argv[2].split(',')[0]
    target_ext = int(sys.argv[2].split(',')[1])
    out_file = spec_file.split('.')[0]+'_weight.fits'
except:
    print('use: python get_weights.py <spec file>,<extension> <target file>,<extension>')
    sys.exit()


# read in spec sample
spec = pyfits.open(spec_file)[spec_ext].data

# read in WL (sub-)sample
target = pyfits.open(target_file)[target_ext].data

# construct dimensions
gr_spec = spec['OMAG_G']-spec['OMAG_R']
ri_spec = spec['OMAG_R']-spec['OMAG_I']
iz_spec = spec['OMAG_I']-spec['OMAG_Z']

# can add mask in here - extend this to a function if needed, to get controlled incompleteness.
mask = np.where((spec['QOP']>=3) & (spec['QOP']<5))[0]

spec_data = np.vstack((spec['OMAG_I'][mask],gr_spec[mask],ri_spec[mask],iz_spec[mask]))
print(spec_data.shape)

gr_targ = target['OMAG_G']-target['OMAG_R']
ri_targ = target['OMAG_R']-target['OMAG_I']
iz_targ = target['OMAG_I']-target['OMAG_Z']

targ_data = np.vstack((target['OMAG_I'],gr_targ,ri_targ,iz_targ))

# scale dimensions? - add later


# get weights for training data
print('calculating weights...')
train_weights = calculate_weights.best_lima_knn_weights(spec_data.T, targ_data.T)


# output spec file with weights
prihdr = pyfits.Header()
prihdu = pyfits.PrimaryHDU(header=prihdr)
c1 = pyfits.Column(name='ID', format='K', array=spec['simID'][mask])
c2 = pyfits.Column(name='Z', format='E', array=spec['Z'][mask])
c3 = pyfits.Column(name='GR', format='E', array=gr_spec[mask])
c4 = pyfits.Column(name='RI', format='E', array=ri_spec[mask])
c5 = pyfits.Column(name='IZ', format='E', array=iz_spec[mask])
c6 = pyfits.Column(name='OMAG_I', format='E', array=spec['OMAG_I'][mask])
c7 = pyfits.Column(name='weight', format='E', array=train_weights)
CC = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
hdu1 = pyfits.BinTableHDU.from_columns(CC)

hdulist = pyfits.HDUList([prihdu,hdu1])
hdulist.writeto(out_file, clobber=True)
