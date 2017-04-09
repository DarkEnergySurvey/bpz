import healpy as hp
import math
import sys
import fitsio
from fitsio import FITS, FITSHDR
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
    fits = FITS(i, 'rw')

    ip = hp.ang2pix(int(Nside), (90.0 - np.array(fits[-1]['DEC'][:])) * math.pi / 180.0, np.array(fits[-1]['RA'][:]) * math.pi / 180.0, nest=0)

    col = 'IP_RING_{:}'.format(Nside)
    fits[-1].insert_column(col, ip)

    print ("Saving to {:} with new column names {:}".format(i, col))