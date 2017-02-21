from __future__ import print_function
import numpy as np
import astropy.io.fits as pyfits 

"""
extensions:
0. spec
1. variance
2. sky
3. telcorr
4. table of values


"""
simspec = pyfits.open("simulated_bcc_vvds.fits")

wave = simspec[2].data['WAVE']
spec = simspec[1].data['SPEC_OBS']
spec[np.isnan(spec)] = 0.
print(spec.shape)
var = (simspec[1].data['SPEC_TRUE']/simspec[1].data['SPEC_SN'])**2
#var[np.isnan(var)] = -1.e-30

spec_len = len(spec[0])
n_spec = len(simspec[1].data['ID'])

#zerosky = np.zeros((n_spec, spec_len))
sky_data = pyfits.open('ESO_sky_VIMOS_LRred_notransm.fits')
zerosky = []
[zerosky.append(sky_data[1].data['sky']) for i in range(n_spec)]
zerosky = np.array(zerosky)
print(zerosky.shape)

tellcorr = np.zeros((spec_len))



# put in the header keywords and data
prihdr = pyfits.Header()
prihdr['OBJECT'] = 'SIMULATED VVDS DEEP'
prihdr['CRVAL1'] = simspec[2].data['WAVE'][0]
#prihdr['CRPIX1'] = 1
prihdr['CRPIX1'] = 1-350
prihdr['CDELT1'] = simspec[2].data['WAVE'][1]-simspec[2].data['WAVE'][0]
prihdu = pyfits.PrimaryHDU(spec[:300,350:950]*1e17,header=prihdr)

hdr1 = pyfits.Header()
hdr1['EXTNAME'] = 'VARIANCE'
hdu1 = pyfits.ImageHDU(data=var[:300,350:950]*1e20,header=hdr1)

hdr2 = pyfits.Header()
hdr2['EXTNAME'] = 'SKY'
#hdu2 = pyfits.ImageHDU(data=zerosky[:300,:],header=hdr1)
hdu2 = pyfits.ImageHDU(data=sky_data[1].data['sky'][350:950]*100.,header=hdr2)

hdr3 = pyfits.Header()
hdr3['EXTNAME'] = 'TELCORR'
#hdu3 = pyfits.ImageHDU(data=var[:300,:],header=hdr3)
hdu3 = pyfits.ImageHDU(data=tellcorr[350:950],header=hdr3)

c1 = pyfits.Column(name='ID', format='K', array=simspec[1].data['ID'])
c2 = pyfits.Column(name='Z', format='E', array=simspec[1].data['Z'])
c3 = pyfits.Column(name='OMAG_I', format='E', array=simspec[1].data['OMAG_I'])
CC = pyfits.ColDefs([c1, c2, c3])
hdu4 = pyfits.BinTableHDU.from_columns(CC)

hdulist = pyfits.HDUList([prihdu,hdu1,hdu2,hdu3,hdu4])
hdulist.writeto('sim_bcc_vvds_deep_ozdesform.fits', clobber=True)

simspec.close()
