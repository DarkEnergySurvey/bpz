#! /usr/bin/env python
"""This helper function should be run, once the redshift predictions have been made. It will massages data according to each specified science working group"""

def apply_cuts(d, sample):
    import numpy as np
    """d = astropy data table
    sample = science sample WL | LSS | GOLD
    """
    cols = {}
    if sample == 'LSS':
        cols['IN_SCIENCE_SAMPLE'] = np.array((d['MEAN_Z'] < 1.0) * (d['MEAN_Z'] > 0.6))
    elif sample == 'WL':
        cols['IN_SCIENCE_SAMPLE'] = np.ones(len(d), dtype=bool)
        cols['WL_WEIGHT'] = np.ones(len(d), dtype=bool)
    return cols

if __name__ == "__main__":
    import sys
    from astropy.io import fits as pyfits

    args = sys.argv[1:]
    sample = args[0]
    files = args[1:]

    for i in files:
        orig_table = pyfits.open(i)[1].data
        orig_cols = orig_table.columns
        cols = apply_cuts(orig_table, sample)
        new_cols = pyfits.ColDefs([pyfits.Column(name=col_name, array=cols[col_name], format='D') for col_name in cols.keys()])

        hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)

        fname = i.replace('.fits', sample + '.fits')
        hdu.writeto(fname)
        print ("Science sample {:} identified in {:}".format(sample, i))
        print ("Saving to {:} with new column names {:}".format(fname, cols.keys()))
