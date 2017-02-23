import astropy.io.fits as pf
import numpy as np
import pylab as mplot
from scipy import interpolate
import sys
import dicts
import random
from tools_simspec import *

"""
Need to sort out number of sky pixels and star frations in dict file. (plus the annoying surveys and dependent files and culled spec samples).
"""

sim_pars = dicts.dict_of_surveys[sys.argv[1]]

#mag1 = 0
#mag2 = 24
N_spec = int(sim_pars['N_sim'] * sim_pars['BCC_frac'] * (1. - sim_pars['star_frac']))

h = 6.626 * 1e-34            # m^2 kg/s
c = 3 * 1e8                  # m/s

expt = 60.*sim_pars['exp_time']              # exposure time in sec
AA = sim_pars['mirror_area']            # collecting area in cm^2
slitloss = sim_pars['slit_loss']               # loss of object light due to slit
Npix = 10                    # number of pixels where object is extracted
sky_pix_fac = 4              # fraction of pixels where sky is estimated comapred to where object is extracted

# DES filters, this is only used for normalising the spectra
L_mean = np.array([ 4735.0, 6385.0, 7755.0, 9225.0, 9950.0 ])      # A
delta_L = np.array([ 1470.0, 1410.0, 1470.0, 1470.0, 500.0 ])      # A

seed = 100                   # seed for Poisson noise
plotting = 0                # weather to make plots later


truth = pf.open(sim_pars['cat_file'])[1].data # random BCC tile
template = pf.open('k_nmf_derived.newdefault.fits') # kcorrect templates
sky = np.loadtxt(sim_pars['sky_file'], comments='#')[:,1]*Npix*expt # sky from ESO
transmission = np.loadtxt(sim_pars['trans_file']) # transmission of VIMOS fom ESO
# ESO ETC website: http://www.eso.org/observing/etc/bin/gen/form?INS.NAME=VIMOS+INS.MODE=SPECTRO
wave = np.loadtxt(sim_pars['wave_array'])*10 # this is the default wavelength grid that we are going to work on

# read in templates
lam         = np.array(template[11].data)
tspec_v0        = np.array(template[1].data)       # no smoothing


# read in coefficients
IDs = truth['ID']
coeffs = truth['COEFFS']
g_mag = truth['TMAG'][:,0] 
r_mag = truth['TMAG'][:,1] 
i_mag = truth['TMAG'][:,2] 
z_mag = truth['TMAG'][:,3] 
y_mag = truth['TMAG'][:,4] 
g_omag = truth['OMAG'][:,0] 
r_omag = truth['OMAG'][:,1] 
i_omag = truth['OMAG'][:,2] 
z_omag = truth['OMAG'][:,3] 
y_omag = truth['OMAG'][:,4] 
z = truth['Z']

#mask = (i_mag<24) # this is the vvds selection cut!
#IDs = IDs[mask]
#coeffs = coeffs[mask]
#g_mag = g_mag[mask]
#r_mag = r_mag[mask]
#i_mag = i_mag[mask]
#z_mag = z_mag[mask]
#y_mag = y_mag[mask]
#g_omag = g_omag[mask]
#r_omag = r_omag[mask]
#i_omag = i_omag[mask]
#z_omag = z_omag[mask]
#y_omag = y_omag[mask]
#z = z[mask]

N = len(z)
random.seed(seed)
#ids = random.randint(N-1, size = N_spec)
ids = random.sample(np.array(range(N)), N_spec)
IDs = IDs[ids]
coeffs = coeffs[ids]
g_mag = g_mag[ids]
r_mag = r_mag[ids]
i_mag = i_mag[ids]
z_mag = z_mag[ids]
y_mag = y_mag[ids]
g_omag = g_omag[ids]
r_omag = r_omag[ids]
i_omag = i_omag[ids]
z_omag = z_omag[ids]
y_omag = y_omag[ids]
z = z[ids]


Spec_new = []
Spec_cal = []
SN = []

for id in range(N_spec):

    zs = z[id]
    mags = i_mag[id]

    # calculate true object spectrum
    spec, tspec = multiply_coeffs(tspec_v0, coeffs[id])
    lam_shift, spec_shift = redshift_spectrum(lam, spec, zs)

    # normalize object spectra
    scale = 10.**((mags - calculate_magnitude(lam_shift, spec_shift, 2))/2.5)
    spec_shift = spec_shift / scale
    
    # interpolate object spectra onto the same default wavelength grid
    f = interpolate.interp1d(lam_shift, spec_shift, bounds_error=False, fill_value=0)
    spec_new = f(wave)

    # convert object spectrum before and after transmission and slitloss into photons
    photons_true = calculate_photons(wave, spec_new, expt, AA)  
    photons_obj = slitloss*calculate_photons(wave, spec_new/100*transmission[:,1], expt, AA)

    # take care of the re-binning
    wave_ave = (wave[1:]+wave[:-1])/2
    trans_ave = ((transmission[:,1][1:]+transmission[:,1][:-1])/2/100)    
    photons_sky = (sky[1:]+sky[:-1])/2 # sky is already in photons, just re-grid
    
    # add Poisson noise on the final photons that reach the focal plane
    np.random.seed(seed)
    photons_poisson_obj = np.random.poisson(lam=photons_obj)
    photons_poisson_sky = np.random.poisson(lam=photons_sky)
    photons_poisson_sky_est = np.random.poisson(lam=photons_sky*sky_pix_fac)/sky_pix_fac
    noise = (photons_obj+photons_sky+photons_sky/sky_pix_fac)**0.5/trans_ave/slitloss 
    # need to check this! this is the noise in the final calibrated object spectra
    
    # reverse the sky-subtracted, calibrated object spectra from photons to F_lambda
    spec_cal = calculate_spec(wave_ave, (photons_poisson_obj+photons_poisson_sky-photons_poisson_sky_est)/trans_ave/slitloss, expt, AA)
    
    # interpolate again onto the same grid
    f = interpolate.interp1d((wave_ave[1:]+wave_ave[:-1])/2, spec_cal, bounds_error=False, fill_value=0)
    spec_cal = f(wave)
    f = interpolate.interp1d(wave_ave, photons_true/noise, bounds_error=False, fill_value=0)
    sn = f(wave)
    
    Spec_new.append(spec_new)
    Spec_cal.append(spec_cal)
    SN.append(sn)


prihdr = pf.Header()
prihdu = pf.PrimaryHDU(header=prihdr)

fmtstr = str(len(wave))+'E'

c1 = pf.Column(name='ID', format='K', array=IDs)
c2 = pf.Column(name='Z', format='E', array=z)
c3 = pf.Column(name='TMAG_G', format='E', array=g_mag) 
c4 = pf.Column(name='TMAG_R', format='E', array=r_mag) 
c5 = pf.Column(name='TMAG_I', format='E', array=i_mag) 
c6 = pf.Column(name='TMAG_Z', format='E', array=z_mag) 
c7 = pf.Column(name='TMAG_Y', format='E', array=y_mag) 
c8 = pf.Column(name='OMAG_G', format='E', array=g_omag) 
c9 = pf.Column(name='OMAG_R', format='E', array=r_omag) 
c10 = pf.Column(name='OMAG_I', format='E', array=i_omag) 
c11 = pf.Column(name='OMAG_Z', format='E', array=z_omag) 
c12 = pf.Column(name='OMAG_Y', format='E', array=y_omag) 
c13 = pf.Column(name='SPEC_TRUE', format=fmtstr, array=np.array(Spec_new)) 
c14 = pf.Column(name='SPEC_OBS', format=fmtstr, array=np.array(Spec_cal))
c15 = pf.Column(name='SPEC_SN', format=fmtstr, array=np.array(SN))

CC = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15]
hdu1 = pf.BinTableHDU.from_columns(CC, nrows=len(IDs))

c1 = pf.Column(name='WAVE', format='E', array=wave)
CC = [c1]
hdu2 = pf.BinTableHDU.from_columns(CC, nrows=len(wave))

hdulist = pf.HDUList([prihdu, hdu1, hdu2])
hdulist.writeto(sim_pars['outfile'], clobber=True)


