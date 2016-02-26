import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval


""" ==================
make data for tests ==
======================
"""


def create_data():
    #make things reproducible
    np.random.seed(0)
    import pandas as pd
    import os
    N = 10000
    zbins = np.arange(300)/300.0 * 2 
    zbins_centers = zbins + zbins[1]/2.0
    #write hdf5 files
    df = pd.DataFrame()
    df['COADD_OBJECTS_ID'] = np.arange(N)

    df['Z_SPEC'] = np.abs(np.random.normal(size=N))
    df['Z_SPEC'] = df['Z_SPEC']/np.amax(df['Z_SPEC']) * 1.89
    
    pdfs = np.zeros((N, 300))
    for i in range(N):
        h = np.histogram(np.abs(np.random.normal(size=1e5) * np.random.uniform()*0.2 + (np.random.uniform()-0.5)*0.1 + df['Z_SPEC'][i]), bins=np.append(zbins, 2))[0]
        pdfs[i] = h

    npdfs = pval.normalisepdfs(pdfs, zbins_centers)

    for i, pdf in enumerate(['pdf_' + str(j) for j in zbins_centers]):
        df[pdf] = npdfs[:, i]

    for i in ['MEAN_Z', 'Z_MC', 'MEDIAN_Z', 'MODE_Z', 'weights_valid']:
        df[i] = df['Z_SPEC'] + np.random.uniform(size=N) * 0.1

    df['WEIGHT'] = np.random.dirichlet(np.arange(N) + N)
    df['MAG_DETMODEL_I'] = np.random.uniform(size=N) * 15 + 15
    df.to_hdf('data/validHDF.hdf5', 'pdf')

    #write an invalid pdf
    #deliberate typo here
    df1 = pd.DataFrame()
    df1['COADDED_OBJECTS_ID'] = np.arange(N)
    for i, pdf in enumerate(['pdf_' + str(j) for j in np.arange(50)/50.0 * 2]):
        df1[pdf] = np.random.dirichlet(np.arange(N) + i)
    df1['Z_SPEC'] = np.random.dirichlet(np.arange(N) + N)
    df1['Z_SPEC'] = df1['Z_SPEC']/np.amax(df1['Z_SPEC']) * 2.0
    df1['WEIGHT'] = np.random.dirichlet(np.arange(N) + N)
    df1['MAG_DETMODEL'] = np.random.uniform(size=N) * 15 + 15
    df1.to_hdf('data/invalidHDF.hdf5', 'pdf')

    np.random.seed(0)
    #create the test fits files
    from astropy.table import Table
    d = {}
    d['Z_SPEC'] = np.random.dirichlet(np.arange(N) + N)
    d['Z_SPEC'] = d['Z_SPEC'] / np.amax(d['Z_SPEC']) * 2.0
    d['COADD_OBJECTS_ID'] = np.arange(N)
    d['MAG_DETMODEL_I'] = np.random.uniform(size=N) * 15 + 15
    d['WEIGHTS'] = np.random.uniform(size=N)
    for i in ['MODE_Z', 'MEAN_Z', 'Z_MC', 'MEDIAN_Z' , 'weights_valid']:
        d[i] = np.random.uniform(size=N) * 2
    fit = Table(d)
    if os.path.exists('data/validPointPrediction.fits'):
        os.remove('data/validPointPrediction.fits')
    fit.write('data/validPointPrediction.fits', format='fits')

    d1 = {}
    d1['Z_SPEC'] = np.random.dirichlet(np.arange(N) + N)
    d1['Z_SPEC'] = d1['Z_SPEC']/np.amax(d1['Z_SPEC']) * 2.0
    d1['COADDED_OBJECTS_ID'] = np.arange(N)
    d1['MAG_DETMODEL_I'] = np.random.uniform(size=N) * 15 + 15
    d1['WEIGHTS'] = np.random.uniform(size=N)
    for i in ['MODE_Z', 'Z_MC']:
        d1[i] = np.random.uniform(size=N) * 2
    fit1 = Table(d1)
    if os.path.exists('data/invalidPointPrediction.fits'):
        os.remove('data/invalidPointPrediction.fits')

    fit1.write('data/invalidPointPrediction.fits', format='fits')


create_data()

"""
import pandas as pd
from astropy.table import Table
filename = 'data/ValidHDF'

filename = 'data/invalidHDF'
df = pd.read_hdf(filename, 'pdf')

filename = 'data/invalidPointPrediction.fits'
df = Table.read(filename)
for i in ['MODE_Z', 'MEAN_Z', 'Z_MC', 'COADD_OBJECTS_ID', 'Z_SPEC']:
    print i in df.keys()

"""