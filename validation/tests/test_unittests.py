import numpy as np
import sys

sys.path.append('../')
import bh_photo_z_validation as pval


"""Units tests for reading hdf5 file format"""


def test_hdf1():
    filename = 'data/__nonValidHDF__'
    err, mess = pval.valid_hdf(filename, args=1)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'file does not exist')


def test_hdf2():
    filename = 'data/invalidHDF'
    err, mess = pval.valid_hdf(filename, args=1)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'missing column COADD_OBJECTS_ID')


#now check for wrong number of tomographic bins
def test_hdf3():
    filename = 'data/validHDF'
    err, mess = pval.valid_hdf(filename, args={'tomographic_bins': np.arange(100)})
    print mess
    np.testing.assert_equal(mess, 'missing column ' + 'pdf_50')
    np.testing.assert_equal(err, False)


def test_hdf4():
    import pandas as pd
    filename = 'data/validHDF'
    err, mess = pval.valid_hdf(filename, args={'tomographic_bins': np.arange(50)})
    np.testing.assert_equal(err, True)

    #To do, fix the data type comparison
    #np.testing.assert_equal(type(mess), pd.core.frame.DataFrame)

"""fits file format testing"""


def test_fits1():
    filename = 'data/NotHerevalidPointPrediction.fits'
    err, mess = pval.valid_fits(filename)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'file does not exist')


def test_fits2():
    filename = 'data/invalidPointPrediction.fits'
    err, mess = pval.valid_fits(filename)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'missing column MEAN_Z')


def test_fits3():
    filename = 'data/validPointPrediction.fits'
    err, mess = pval.valid_fits(filename)
    np.testing.assert_equal(err, True)


"""metric unit tests"""


def test_sigma68_1():
    arr1 = np.arange(1000)
    np.random.shuffle(arr1)
    print 2 * pval.sigma_68(arr1)
    np.testing.assert_equal(np.abs(2 * pval.sigma_68(arr1) - 681.8) < 0.1, True)


def test_sigma68_2():
    arr2 = np.arange(1e4)*1.0
    np.random.shuffle(arr2)
    print 2 * pval.sigma_68(arr2)
    np.testing.assert_equal(np.abs(2 * pval.sigma_68(arr2) - 6824) < 1, True)


def test_sigma68_3():
    arr1 = np.random.normal(size=1e7)
    np.testing.assert_equal(np.abs(pval.sigma_68(arr1) - np.std(arr1)) < 0.01, True)


def test_sigma95_1():
    arr1 = np.random.normal(size=1e7)
    np.testing.assert_equal(np.abs(pval.sigma_95(arr1) - 2 * np.std(arr1)) < 0.01, True)


def test_sigma95_2():
    arr1 = np.arange(1000)
    np.random.shuffle(arr1)
    print 2 * pval.sigma_95(arr1)
    np.testing.assert_equal(np.abs(2 * pval.sigma_95(arr1) - 953) < 1, True)


def test_sigma95_3():
    arr2 = np.arange(1e4)*1.0
    np.random.shuffle(arr2)
    print 2 * pval.sigma_95(arr2)
    np.testing.assert_equal(np.abs(2 * pval.sigma_95(arr2) - 9540) < 1, True)


def test_delta_z():
    for N in [1, 100, 1000]:
        deltaz = pval.delta_z(np.arange(N), np.arange(N))
        np.testing.assert_equal(np.sum(np.abs(deltaz)), 0)


def test_outlierRate():
    #hard coded outlier rate
    arr1 = np.random.normal(size=1e6) * 0.15
    outRate = pval.outlierRate(arr1)
    np.testing.assert_equal(outRate + 0.68 > 0.99, True)


def test_outlierFraction():
    #hard coded outlier rate
    arr1 = np.random.normal(size=1e6) * 0.15
    outRate = pval.outlierFraction(arr1)
    np.testing.assert_equal(outRate + 68 > 99.7, True)


"""tests on pdf """

def test_normalize_pdf():
    N = 1e4
    pdf = np.random.dirichlet(np.arange(N))
    z = np.arange(N, dtype=float) / N
    normPdf = pval.normalize_pdf(pdf, z)
    integ = np.trapz(normPdf, x=z)
    np.testing.assert_equal(np.abs(integ - 1) < 0.001, True)


def create_data():
    #make things reproducible
    np.random.seed(0)
    import pandas as pd
    import os
    N = 1000
    #write hdf5 files
    df = pd.DataFrame()
    df['COADD_OBJECTS_ID'] = np.arange(N)

    for i, pdf in enumerate(['pdf_' + str(j) for j in range(50)]):
        df[pdf] = np.random.dirichlet(np.arange(N) + i)
    df['Z_SPEC'] = np.random.dirichlet(np.arange(N) + N)
    df['WEIGHT'] = np.random.dirichlet(np.arange(N) + N)
    df.to_hdf('data/validHDF', 'pdf')

    #write an invalid pdf
    #deliberate typo here
    df1 = pd.DataFrame()
    df1['COADDED_OBJECTS_ID'] = np.arange(N)
    for i, pdf in enumerate(['pdf_' + str(j) for j in range(3)]):
        df1[pdf] = np.random.dirichlet(np.arange(N) + i)
    df1['Z_SPEC'] = np.random.dirichlet(np.arange(N) + N)
    df1['WEIGHT'] = np.random.dirichlet(np.arange(N) + N)
    df1.to_hdf('data/invalidHDF', 'pdf')

    np.random.seed(0)
    #create the test fits files
    from astropy.table import Table
    d = {}
    d['Z_SPEC'] = np.random.dirichlet(np.arange(N) + N)
    d['COADD_OBJECTS_ID'] = np.arange(N)
    for i in ['MODE_Z', 'MEAN_Z', 'Z_MC']:
        d[i] = np.random.uniform(size=N)*2
    fit = Table(d)
    if os.path.exists('data/validPointPrediction.fits'):
        os.remove('data/validPointPrediction.fits')
    fit.write('data/validPointPrediction.fits', format='fits')

    d1 = {}
    d1['Z_SPEC'] = np.random.dirichlet(np.arange(N) + N)
    d1['COADDED_OBJECTS_ID'] = np.arange(N)
    for i in ['MODE_Z', 'Z_MC']:
        d1[i] = np.random.uniform(size=N)*2
    fit1 = Table(d1)
    if os.path.exists('data/invalidPointPrediction.fits'):
        os.remove('data/invalidPointPrediction.fits')

    fit1.write('data/invalidPointPrediction.fits', format='fits')



"""
create_data()
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