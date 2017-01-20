import numpy as np
import healpy as hp

gdmask = hp.read_map('y1a1_gold_1.0.2_wide_footprint_4096.fits.gz')
ipring, = np.where(gdmask >= 1)

msk = np.zeros_like(gdmask, dtype=bool)
msk[gdmask >= 1] = True

hp.write_map('y1a1_gold_1.0.2_wide_footprint_4096.mask.fits', msk)
