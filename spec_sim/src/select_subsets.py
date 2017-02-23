from __future__ import print_function
import numpy as np
import astropy.io.fits as pyfits 
import random
import math as m
import sys
import dicts

"""
Want 10% overlap between subsamples, and a total of around 12000 objects in subsamples of 300 (i.e. 40 files).
There are dependent files. So we want this as a module per survey.

"""

def save_field(simspec, ids, outname, skyfile):
    wave = simspec[2].data['WAVE']
    # let's trim to 9800 - only for vimos, change to something more sensible in future.
    upper_lim = np.where(wave > 9800.)[0][0]
    spec = simspec[1].data['SPEC_OBS'][ids,:]
    spec[np.isnan(spec)] = 0.
    print(spec.shape)
    var = (simspec[1].data['SPEC_TRUE'][ids,:]/simspec[1].data['SPEC_SN'][ids,:])**2
    #var[np.isnan(var)] = -1.e-30

    spec_len = len(spec[0])
    n_spec = len(ids)

    ##zerosky = np.zeros((n_spec, spec_len))
    #sky_data = pyfits.open(skyfile)
    sky_data = np.loadtxt(sim_pars['sky_file'], comments='#')[:upper_lim,1]
    #zerosky = []
    #[zerosky.append(sky_data[1].data['sky']) for i in range(n_spec)]
    #zerosky = np.array(zerosky)
    #print(zerosky.shape)
    sky_data = np.array(sky_data)

    tellcorr = np.zeros((spec_len))

    # put in the header keywords and data
    prihdr = pyfits.Header()
    prihdr['OBJECT'] = 'SIMULATED SPECTRA'
    prihdr['CRVAL1'] = simspec[2].data['WAVE'][0]
    #prihdr['CRPIX1'] = 1
    prihdr['CRPIX1'] = 1-350 # ???? Oh yeah, I cut the end off the spectra, 350 off the front, 100 off the back.
    # do we still want to trim the VIMOS MR?? Or add more params? 350 is ok. Will need to check future surveys
    prihdr['CDELT1'] = simspec[2].data['WAVE'][1]-simspec[2].data['WAVE'][0]
    #prihdu = pyfits.PrimaryHDU(spec[:,350:-100]*1e17,header=prihdr)
    prihdu = pyfits.PrimaryHDU(spec[:,350:upper_lim]*1e17,header=prihdr)
    
    hdr1 = pyfits.Header()
    hdr1['EXTNAME'] = 'VARIANCE'
    hdu1 = pyfits.ImageHDU(data=var[:,350:upper_lim]*1e20,header=hdr1)

    hdr2 = pyfits.Header()
    hdr2['EXTNAME'] = 'SKY'
    #hdu2 = pyfits.ImageHDU(data=sky_data[1].data['sky'][350:950]*100.,header=hdr2)
    hdu2 = pyfits.ImageHDU(data=sky_data[350:upper_lim]*100.,header=hdr2)

    hdr3 = pyfits.Header()
    hdr3['EXTNAME'] = 'TELCORR'
    hdu3 = pyfits.ImageHDU(data=tellcorr[350:upper_lim],header=hdr3)

    c1 = pyfits.Column(name='ID', format='K', array=simspec[1].data['ID'][ids])
    c2 = pyfits.Column(name='Z', format='E', array=simspec[1].data['Z'][ids])
    c3 = pyfits.Column(name='OMAG_G', format='E', array=simspec[1].data['OMAG_G'][ids])
    c4 = pyfits.Column(name='OMAG_R', format='E', array=simspec[1].data['OMAG_R'][ids])
    c5 = pyfits.Column(name='OMAG_I', format='E', array=simspec[1].data['OMAG_I'][ids])
    c6 = pyfits.Column(name='OMAG_Z', format='E', array=simspec[1].data['OMAG_Z'][ids])
    c7 = pyfits.Column(name='OMAG_Y', format='E', array=simspec[1].data['OMAG_Y'][ids])
    CC = pyfits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
    hdu4 = pyfits.BinTableHDU.from_columns(CC)

    hdulist = pyfits.HDUList([prihdu,hdu1,hdu2,hdu3,hdu4])
    hdulist.writeto(outname, clobber=True)



# read in dict
sim_pars = dicts.dict_of_surveys[sys.argv[1]]

overall_spec_num = 75623 # number from all VIMOS surveys - constant across this set-up of surveys. 
input_file = sim_pars['outfile']
sky_file = sim_pars['sky_file']
seed = 100
random.seed(seed)

# read file
simspec = pyfits.open(input_file)

# work out full size of sample, and define fraction going into subsamples
spec_num = len(simspec[1].data['SPEC_OBS'])
print(spec_num)
target_num = 12000. * spec_num / overall_spec_num
n_bundles = int(m.ceil(target_num / 300.))
print(n_bundles)

# randomly select objects for subsamples and package into files
ids = []
for i in range(n_bundles):
    ids.append(random.sample(np.array(range(spec_num)), 300))

# how to ensure 10% overlap?
# for each file, we check how many objects are in any of the others.
# If not enough, we random sample from the ID list

if n_bundles > 1:
    for i in range(n_bundles):
        tmp = np.array(ids)[np.where(np.arange(len(ids))!=i)[0],:]
        n_repeat = len(np.intersect1d(ids[i],tmp)) # might need unique here, no seems not.
        print(n_repeat,len(np.unique(np.intersect1d(ids[i],tmp))))

        if n_repeat < 30:
            n_needed = int(m.ceil(30. - n_repeat))
            id_add = random.sample(np.reshape(tmp,(-1)), n_needed)
            #ids[i].append(id_add)
            ids_tmp = np.hstack((ids[i],id_add))
        else:
            ids_tmp = ids[i]

        # loop over bundles - make unique names 
        outname = sim_pars['outfile'].split('.')[0]+"_field"+str(i)+".fits"
        save_field(simspec, ids_tmp, outname, sim_pars['sky_file'])
        
else:
    outname = sim_pars['outfile'].split('.')[0]+"_field1.fits"
    save_field(simspec, ids[0], outname, sim_pars['sky_file'])
