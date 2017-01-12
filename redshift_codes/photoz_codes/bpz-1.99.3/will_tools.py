import numpy as np
import astropy.io.fits as pyfits

# read 2D array using fits (topcat style, 1st extension)

def get_2Darray_fromfits(filename,cols='all',nrows='all',verbose='no'):
    data = pyfits.open(filename)[1].data
    if cols=='all':
        # just convert data to a 2Darray
        array_out = np.array(data)
    else:
        nc=len(cols)
        array_out = None
        for col in cols:
            array_in = data[col]
            array_out = np.hstack((array_out,array_in))
        array_out = array_out[1:].reshape((nc,len(array_in))).T
    return np.array(array_out,dtype=float)

def get_long_fromfits(filename,col='COADD_OBJECTS_ID',nrows='all',verbose='no'):
    data = pyfits.open(filename)[1].data
    try:
        array_in = data[col]
    except:
        print('Column not found:'+col)
    return np.array(array_in, dtype=long)
