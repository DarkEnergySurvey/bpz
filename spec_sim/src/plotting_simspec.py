
    # plotting ##################
    if plotting==1:
   
        # plot the object spectra before and after transmission
        mplot.figure(figsize=(8,16))
        mplot.subplot(411)
        mplot.plot(wave, spec_new, label='galaxy spectra')   
        mplot.plot(wave, spec_new/100*transmission[:,1], label='galaxy + transmission')
        mplot.xlabel('A')
        mplot.ylabel('$F_{\lambda}$ (erg/s/cm^2/A)')
        mplot.xlim(5500, 10000)
        mplot.legend(loc='upper left')
        mplot.title('z='+str(zs)+', i_mag='+str(mags))

        # plot the measured spectra before and after adding Poisson noise
        mplot.subplot(412)
        mplot.plot(wave_ave, photons_poisson_obj+photons_poisson_sky, color='r', label='poisson noise')
        mplot.plot(wave_ave, photons_obj+photons_sky, color='k', label='mean')
        mplot.xlabel('A')
        mplot.ylabel('# of photons')
        mplot.xlim(5500, 10000)
        mplot.ylim(0,600000)
        mplot.legend(loc='upper left')
        mplot.title('Observed spectra in photons')

        # plot the calibrated, sky-subtracted spectra and the true spectra in F_lambda
        # this is what we use to measure redshifts 
        mplot.subplot(413)
        mplot.plot(wave, spec_new, color='k', label='true', ls='--')
        mplot.plot(wave, spec_cal, color='b', label='calibrated', alpha=0.4)
        mplot.xlabel('A')
        mplot.ylabel('$F_{\lambda}$ (erg/s/cm^2/A)')
        mplot.xlim(5500, 10000)
        mplot.ylim(-0.1*np.max(spec_new[(wave>5500)*(wave<10000)]), 1.1*np.max(spec_new[(wave>5500)*(wave<10000)]))
        mplot.legend(loc='upper left')
        mplot.title('Calibrated vs. true spectra in F_lambda')

        # plot the S/N
        mplot.subplot(414)
        mplot.plot(wave_ave, photons_true/noise, color='k', label='true')
        mplot.xlabel('A')
        mplot.ylabel('S/N')  
        mplot.xlim(5500, 10000)

        mplot.tight_layout()
        mplot.savefig('bcc_noise_spectra_'+str(id)+'.png')



A = pf.open('simulated_bcc_vvds.fits')[1].data
mplot.plot(wave, A['SPEC_TRUE'][100], color='k', ls='--')
mplot.plot(wave, A['SPEC_OBS'][100], color='b', alpha=0.4)
mplot.xlim(5500, 10000)
mplot.xlabel('A')
mplot.ylabel('$F_{\lambda}$ (erg/s/cm^2/A)')
mplot.ylim(-0.1*np.max(A['SPEC_TRUE'][100][(wave>5500)*(wave<10000)]), 1.1*np.max(A['SPEC_TRUE'][100][(wave>5500)*(wave<10000)]))
