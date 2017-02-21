# dictionary of params for each survey (object selection done separately)
#
# Param meanings:
#
# exp_time: exposure time in minutes
# N_sim: Number of objects lying in the DES footprint that were targetted by the survey
# sky_file: typical sky in photons/s/pixel (along slit). Includes full transmission function. 
# transmission_file: % transmission, fn of wavelength
# wave_array: wavelength array (per pixel)
# slit_loss: actually fraction *not* lost (current values are very rough!) 
# mirror_area: effective collecting area (if avail), other wise just area of primary mirror.
# outfile: output file...
# cat_file: Simulation file that contains the template co-effs, with photometric selection function and survey area already implemented.
# BCC_frac: the fraction of 'targetted' galaxies that could possibly be in cat_file (will be less than 1 if survey is deeper than simulation for instance).
# star_frac: fraction of targets that are likely to be stars 
#

# VVDS WIDE
dict_of_surveys = {'vvds_wide1_dict' : {'exp_time': 45,
                                        'N_sim': 13300,
                                        'sky_file': "ESO_sky_VIMOS_LRred_wtransm.txt",
                                        'trans_file': "VIMOS_LRred_transmission.txt",
                                        'wave_array': "VIMOS_LRred_wave.dat",
                                        'slit_loss': 0.7,
                                        'mirror_area': 51.86*100**2,
                                        'outfile': "vvds_wide1_sim.fits",
                                        'cat_file': "Buzzard_v1.1_truth_vvds_wide1.fit",
                                        'BCC_frac': 1.,
                                        'star_frac': 0.},

                   'vvds_wide2_dict' : {'exp_time': 45,
                                        'N_sim': 1750, # rough fraction of 6229 that is in footprint
                                        'sky_file': "ESO_sky_VIMOS_LRred_wtransm.txt",
                                        'trans_file': "VIMOS_LRred_transmission.txt",
                                        'wave_array': "VIMOS_LRred_wave.dat",
                                        'slit_loss': 0.7,
                                        'mirror_area': 51.86*100**2,
                                        'outfile': "vvds_wide2_sim.fits",
                                        'cat_file': "Buzzard_v1.1_truth_vvds_wide2.fit",
                                        'BCC_frac': 1.,
                                        'star_frac': 0.},

                   'vvds_wide3_dict' : {'exp_time': 45,
                                        'N_sim': 6500, # missing a fraction (of 7358), again rough.
                                        'sky_file': "ESO_sky_VIMOS_LRred_wtransm.txt",
                                        'trans_file': "VIMOS_LRred_transmission.txt",
                                        'wave_array': "VIMOS_LRred_wave.dat",
                                        'slit_loss': 0.7,
                                        'mirror_area': 51.86*100**2,
                                        'outfile': "vvds_wide3_sim.fits",
                                        'cat_file': "Buzzard_v1.1_truth_vvds_wide3.fit",
                                        'BCC_frac': 1.,
                                        'star_frac': 0.},

                   # VVDS DEEP
                   'vvds_deep1_dict' : {'exp_time': 270,
                                        'N_sim': 11360,
                                        'sky_file': "ESO_sky_VIMOS_LRred_wtransm.txt",
                                        'trans_file': "VIMOS_LRred_transmission.txt",
                                        'wave_array': "VIMOS_LRred_wave.dat",
                                        'slit_loss': 0.7,
                                        'mirror_area': 51.86*100**2,
                                        'outfile': "vvds_deep1_sim.fits",
                                        'cat_file': "Buzzard_v1.1_truth_vvds_deep1.fit",
                                        'BCC_frac': 1.,
                                        'star_frac': 0.},
                   
                   'vvds_deep2_dict' : {'exp_time': 270,
                                        'N_sim': 1572,
                                        'sky_file': "ESO_sky_VIMOS_LRred_wtransm.txt",
                                        'trans_file': "VIMOS_LRred_transmission.txt",
                                        'wave_array': "VIMOS_LRred_wave.dat",
                                        'slit_loss': 0.7,
                                        'mirror_area': 51.86*100**2,
                                        'outfile': "vvds_deep2_sim.fits",
                                        'cat_file': "Buzzard_v1.1_truth_vvds_deep2.fit",
                                        'BCC_frac': 1.,
                                        'star_frac': 0.},

                   # VIPERS
                   'vipers_dict' : {'exp_time': 45,
                                    'N_sim': 20452,
                                    'sky_file': "ESO_sky_VIMOS_LRred_wtransm.txt",
                                    'trans_file': "VIMOS_LRred_transmission.txt",
                                    'wave_array': "VIMOS_LRred_wave.dat",
                                    'slit_loss': 0.7,
                                    'mirror_area': 51.86*100**2,
                                    'outfile': "vipers_sim.fits",
                                    'cat_file': "Buzzard_v1.1_truth_vipers.fit",
                                    'BCC_frac': 1.,
                                    'star_frac': 0.0633},

                   # DEEP2
                   'deep2_dict' : {'exp_time': 60,
                                   'N_sim': 11986,
                                   'sky_file': "KECK_sky_DIEMOS_1200_wtransm.txt",
                                   'trans_file': "DIEMOS_1200_transmission.txt",
                                   'wave_array': "DIEMOS_1200_wave.dat",
                                   'slit_loss': 0.7,
                                   'mirror_area': 73.4*100**2,
                                   'outfile': "deep2_sim.fits",
                                   'cat_file': "Buzzard_v1.1_truth_deep2.fit",
                                   'BCC_frac': 1.,
                                   'star_frac': 0.},

                   # zCOSMOS
                   'zcosmos_dict' : {'exp_time': 60,
                                     'N_sim': 20689,
                                     'sky_file': "ESO_sky_VIMOS_MR_wtransm.txt",
                                     'trans_file': "VIMOS_MR_transmission.txt",
                                     'wave_array': "VIMOS_MR_wave.dat",
                                     'slit_loss': 0.7,
                                     'mirror_area': 51.86*100**2,
                                     'outfile': "zcosmos_sim.fits",
                                     'cat_file': "Buzzard_v1.1_truth_zcosmos.fit",
                                     'BCC_frac': 1.,
                                     'star_frac': 0.0488},

                   # WIGGLEZ
                   'wigglez_dict' : {'exp_time': 60,
                                     'N_sim': 20327,
                                     'sky_file': "AAT_sky_AAO_wtransm.txt",
                                     'trans_file': "AAO_transmission.txt",
                                     'wave_array': "AAO_wave.dat",
                                     'slit_loss': 1.0,
                                     'mirror_area': 12.*100**2,
                                     'outfile': "wigglez_sim.fits",
                                     'cat_file': "Buzzard_v1.1_truth_wigglez.fit",
                                     'BCC_frac': 1.,
                                     'star_frac': 0.},

                   # 3D-HST
                   '3dhst_dict' : {'exp_time': 0.,
                                     'N_sim': 0.,
                                     'sky_file': "0.",
                                     'trans_file': "0.",
                                     'wave_array': "0.",
                                     'slit_loss': 0.,
                                     'mirror_area': 0.,
                                     'outfile': "3dhst_sim.fits",
                                     'cat_file': "Buzzard_v1.1_truth_3dhst.fit",
                                     'BCC_frac': 1.,
                                     'star_frac': 0.},

                   # ACES
                   'aces_dict' : {'exp_time': 30, # check that they haven't been combined yet.
                                     'N_sim': 13963,
                                     'sky_file': "Magellan_sky_IMACS_wtransm.txt",
                                     'trans_file': "IMACS_transmission.txt",
                                     'wave_array': "IMACS_wave.dat",
                                     'slit_loss': 0.7,
                                     'mirror_area': 33.18*100**2,
                                     'outfile': "aces_sim.fits",
                                     'cat_file': "Buzzard_v1.1_truth_aces.fit",
                                     'BCC_frac': 1.,
                                     'star_frac': 0.},
                   
                   # DES_AAOMEGA
                   'des_aaom_dict' : {'exp_time': 0,
                                     'N_sim': 0,
                                     'sky_file': "AAT_sky_AAO_wtransm.txt",
                                     'trans_file': "AAO_transmission.txt",
                                     'wave_array': "AAO_wave.dat",
                                     'slit_loss': 1.,
                                     'mirror_area': 12.*100**2,
                                     'outfile': "des_aaomega_sim.fits",
                                     'cat_file': "Buzzard_v1.1_truth_des_aaomega.fit",
                                     'BCC_frac': 1.,
                                     'star_frac': 0.},

                   # EBOSS_DES_ELG
                   'eboss_dict' : {'exp_time': 0,
                                     'N_sim': 0,
                                     'sky_file': "Apache_sky_eboss_wtransm.txt",
                                     'trans_file': "eboss_transmission.txt",
                                     'wave_array': "eboss_wave.dat",
                                     'slit_loss': 1.,
                                     'mirror_area': 15.9*100**2,
                                     'outfile': "eboss_des_elg_sim.fits",
                                     'cat_file': "Buzzard_v1.1_truth_eboss_des_elg.fit",
                                     'BCC_frac': 1.,
                                     'star_frac': 0.}}
