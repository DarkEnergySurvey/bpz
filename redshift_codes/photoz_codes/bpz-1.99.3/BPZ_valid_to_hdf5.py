from collections import namedtuple
from operator import itemgetter
from pprint import pformat
import numpy as np
import math as m
import astropy.io.fits as pyfits
import copy
import random as rdm
import pandas as pd

"""
Revamping....
- everything is by line indexing, so this has to hold!
- add optional compression
"""

# some constants
dz = 0.01
z_max = 3.5

# files
prob_file = "valid_103.probs"
bpz_file = "valid_103.bpz"
weight_file = "valid_103_sample_allcols.fits"
hdf5_file = "PHOTOZ_BPZ_Y1_v0.11_VALID_103_z35.hdf5"
fits_file = "PHOTOZ_BPZ_Y1_v0.11_VALID_103_z35.fits"

def get_mc(pdf,zarr):
    # renorm incase there is probability at higher-z that we've cut, or some error.
    if np.sum(pdf) > 0:
        pdf = pdf/np.sum(pdf)
        cum = np.cumsum(pdf)
        targ_prob = rdm.random()
        indx = np.where(cum >= targ_prob)[0]
        zmc_ob = zarr[indx[0]-1] + (zarr[indx[0]] - zarr[indx[0]-1]) * ( (targ_prob - cum[indx[0]-1]) / (cum[indx[0]] - cum[indx[0]-1]) )
        return zmc_ob
    else:
        return -1.

def get_mean(pdf,zarr):
    if np.sum(pdf) > 0:
        zm = np.average(zarr, weights=pdf)
        return zm
    else:
        return -1.

def get_median(pdf,zarr):
    if np.sum(pdf) > 0:
        pdf = pdf/np.sum(pdf)
        cum = np.cumsum(pdf)
        indx = np.where(cum >= 0.5)[0]
        zmed_ob = zarr[indx[0]-1] + (zarr[indx[0]] - zarr[indx[0]-1]) * ( (0.5 - cum[indx[0]-1]) / (cum[indx[0]] - cum[indx[0]-1]) )
        return zmed_ob
    else:
        return -1.

# read prob file
ID = []
pz = []
with open(prob_file,'r') as f:
    for line in f:
        if line[0] != "#":
            ID.append(line.split()[0])
            pz.append(line.split()[1:int(z_max/dz)+1])
print("read prob file...")

# read bpz file
zmode = []
zspec = []
with open(bpz_file,'r') as f:
    for line in f:
        if line[0] != "#":
            zmode.append(line.split()[1])
            zspec.append(line.split()[9])
print("read bpz file...")

# read weights (assumes topcat fits output)
weight_data = pyfits.open(weight_file)[1].data
WLweights = weight_data['WL_VALID_WEIGHTS']
LSSweights = weight_data['LSS_VALID_WEIGHTS']

#ID = np.array(ID, dtype=long)
# get ID from fits file instead
ID = weight_data['COADD_OBJECTS_ID']
pz = np.array(pz, dtype=float)
zmode = np.array(zmode, dtype=float)
#zspec = np.array(zspec, dtype=float)
# zspec is weird too... (fixed now - was reading in wrong file.)
zspec = weight_data['Z']

zarr = np.arange(dz,z_max+dz,dz)

zmc = [get_mc(pz[i],zarr) for i in range(len(ID))]
zmean = [get_mean(pz[i],zarr) for i in range(len(ID))]
zmed = [get_median(pz[i],zarr) for i in range(len(ID))]
zmc = np.array(zmc)
zmean = np.array(zmean)
zmed = np.array(zmed)

#print(ID.shape,zmean.shape,zmode.shape,zmed.shape,zmc.shape,zspec.shape,WLweights.shape,LSSweights.shape,pz.shape)
data = np.vstack((ID,zmean,zmode,zmed,zmc,zspec,WLweights,LSSweights))
data = np.vstack((data,pz.T)).T
print(data.shape)
df = pd.DataFrame.from_records(data)

zbins = (np.arange(int(z_max/dz))+1)/(z_max/dz) * z_max
pdf_names = ['PDF_' + str(j) for j in zbins]
df.columns = ['COADD_OBJECTS_ID','MEAN_Z','MODE_Z','MEDIAN_Z','Z_MC','Z_SPEC','WL_WEIGHTS','LSS_WEIGHTS'] + pdf_names
print(df.head())

print("writing...")
df.to_hdf(hdf5_file, 'PDF', format='table', complib='zlib', complevel=5)


#fits output is deprecated....
col0 = pyfits.Column(name='COADD_OBJECTS_ID', format='K', array=ID)
col1 = pyfits.Column(name='MEAN_Z', format='E', array=np.array(zmean))
col2 = pyfits.Column(name='MODE_Z', format='E', array=np.array(zmode))
col3 = pyfits.Column(name='MEDIAN_Z', format='E', array=np.array(zmed))
col4 = pyfits.Column(name='Z_MC', format='E', array=np.array(zmc))
col5 = pyfits.Column(name='Z_SPEC', format='E', array=np.array(zspec))
col6 = pyfits.Column(name='WL_WEIGHTS', format='E', array=np.array(WLweights))
col7 = pyfits.Column(name='LSS_WEIGHTS', format='E', array=np.array(LSSweights))
pdfstr = str(int(z_max/dz))+'E'
col8 = pyfits.Column(name='PDF', format=pdfstr, array=pz)
newcols = pyfits.ColDefs([col0,col1,col2,col3,col4,col5,col6,col7,col8])

tbhdu = pyfits.BinTableHDU.from_columns(newcols)
tbhdu.writeto(fits_file, clobber=True)

