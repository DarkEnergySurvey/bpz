import healpy as hp
import math
import sys
#from astropy.io import fits
from astropy.table import Table

"""
Written by ben hoyle.

Add IP_RING_XXX column to files

call like
python add_ip_ring.py Nside PathToFile[s]

Where Nside should probaby==128 for our purposes.

Notes:
Current outputs a new file, replicating the first, with the file ending .IPXXX.fits
To do
improve rather than output a new file, add a new (unique) column with the IPRING_XXX column to the old file.
"""

args = sys.argv[1:]

if len(args) < 2:
    print ('add_ip_ring.py Nside PathTo*Files')
    print ('generates a column IP_RING_XXX and saves to the fits file. It adds .IPXXX.fits to the file name')
    sys.exit()

Nside = args[1]
f = args[2:]
for i in f:
    d = Table.read(i)
    d['IP_RING_{:}'.format(Nside)] = hp.ang2pix(int(Nside), (90.0 - d['DEC']) * math.pi / 180.0, d['RA'] * math.pi / 180.0, nest=0)
    d.write(i + '.IP{:}.fits'.format(Nside))
