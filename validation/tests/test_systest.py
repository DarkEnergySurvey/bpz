import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval
sys.path.append('../../systematic_tests')
import one_d_correlate as oneD


""" =============================
tests on hdf5 file format  ======
=================================
"""


def test_corr1():
    """test get correct error with non-existent hdf5 file """
    filename = 'data/__nonValidHDF__.hdf5'
    cols = ['COADD_OBJECTS_ID', 'Z_SPEC', 'pdf_0']
    err, mess = pval.valid_file(filename, cols)
    np.testing.assert_equal(err, False)
    np.testing.assert_equal(mess, 'file does not exist')
