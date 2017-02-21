import healpy as hp
import math
import sys
from astropy.io import fits as pyfits
import os

"""
Written by ben hoyle.
Add IP_RING_XXX column to files
"""

def help():
    print ("python add_ip_ring.py Nside PathToFile[s]")

    print ("Where Nside should probaby==128 for our purposes.")

args = sys.argv[1:]

if len(args) < 2:
    print ('add_ip_ring.py Nside PathTo*Files')
    print ('generates a column IP_RING_XXX and saves to the fits file. It adds .IPXXX.fits to the file name')
    sys.exit()

Nside = args[1]
f = args[2:]

sample = str(Nside)

for i in f:
    orig_table = pyfits.open(i)[1].data
    orig_cols = orig_table.columns

    cols = {}
    cols['IP_RING_{:}'.format(Nside)] = hp.ang2pix(int(Nside), (90.0 - orig_table['DEC']) * math.pi / 180.0, orig_table['RA'] * math.pi / 180.0, nest=0)

    new_cols = pyfits.ColDefs([pyfits.Column(name=col_name, array=cols[col_name], format='D') for col_name in cols.keys()])

    hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)

    fname = i.replace('.fits', '.IP{:}.fits'.format(Nside))
    hdu.writeto(fname)
    print ("Saving to {:} with new column names {:}".format(i, cols.keys()))
    import os
    os.rename(fname, i)