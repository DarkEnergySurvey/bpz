#import lib.plots.plotLib as plotLib
from pylab import *
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-darkgrid')

import pyfits as pf
import glob
import yaml 
import pandas as pd
import bh_photo_z_validation_au as bh
from scipy.stats import gaussian_kde

'''
python make_real.py z_pht=MEAN_Z path=/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/train_1_noweight/ weight_s=WL_WEIGHTS
'''

"""
Plot many different n(z) together reading from a directory, where files are given in fits format following the DES photo-z wg standard in fits file

Author: Aurelio Carnero, Julia Gschwend

-input:
python plot_nz.py z_pht=MEAN_Z path='valid_1_noweight/' path_true='sims_spectra_representative_valid.WL_LSS_FLAGS.fits' [z_true=Z_TRUTH]

-outputs:
3 plots for now, the individual n(z) estimations and average, plus the average and error, and error.
"""
def load_yaml(filename):

    try:
        d = yaml.load(open(filename, 'r'))
        return d

    except:
        print "error loading yaml file " + filename
        print "check format here http://yaml-online-parser.appspot.com/"
        print "aborting"
        sys.exit()


args = sys.argv[1:]

inArgs = {}
for i in args:
    if '='in i:
        k, v = i.split('=')
        inArgs[k] = v

if 'path' not in inArgs.keys():
    print 'missing path for fits point statistics'
    print 'python plot_nz.py z_pht=MEAN_Z path=/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/train_1_noweight/valid_1_noweight/ path_true=sims_spectra_representative_valid.WL_LSS_FLAGS.fits [z_true=Z_TRUTH]'
    sys.exit('ERROR: path missing')

path = inArgs['path']

if 'z_pht' not in inArgs.keys():
    print 'missing photo-z estimation mode'
    print 'the options are MEAN_Z, MEDIAN_Z, MODE_Z or Z_MC'
    sys.exit('ERROR: missing photo-z estimation mode')
else:
    z_pht = inArgs['z_pht']

if 'z_true' in inArgs.keys():
    z_true = inArgs['z_true']
else:
    z_true = 'Z_TRUTH'

if 'weight_s' in inArgs.keys():
    weight_s = inArgs['weight_s']
    print weight_s
else:
    weight_s = False
    print 'WARNING no weights'



dim_lin = 6

bin_edge = np.linspace(0.4, 1., dim_lin+1)
bin_center = []
print bin_edge
for i in range(dim_lin):
	bin_center.append(bin_edge[i] + (bin_edge[1]-bin_edge[0])/2.)
print bin_center

pdf_results = []
results_stats = []

cla()
clf()
import random
from scipy import interpolate

x_plot = np.linspace(0., 2., 50.)
dim_algorithm = len(glob.glob(path + '*.hdf5'))

def random_distr(l,zmin,zmax):
    r = random.uniform(zmin, zmax)
    pd = l(r)
    random.seed()
    rr = random.uniform(0, 1)
    if rr <= pd:
	return r
    else:
	return 'non'
    

for i in glob.glob(path + '*.hdf5'):
	lab = i.split('/')[-1].split('_')[0]

	store = pd.HDFStore(i)
	df_orig = store['pdf']
	df = store['pdf']
	weight = df_orig[weight_s]
	z_spec = df_orig['Z_SPEC']
	store.close()
	centers = np.array([float(name[4:]) for name in df.columns if 'pdf' in name])

	y_pdf = np.array([name[:] for name in df.columns if 'pdf' in name])
	df = df[y_pdf]

	y_pdf = bh.normalisepdfs(df,centers)
	for index, row in y_pdf.iterrows():
		print len(row)
		print len(centers)
#		print df['Z_SPEC'][index]
		z = z_spec[index]
		print z
		pdf = []
		zmin = np.min(centers)
		zmax = np.max(centers)
		row_max = np.max(row)
		row = row/row_max
		print np.max(row)
#		for cc, rr in zip(centers, row):
#			pdf.append((cc, rr))
#		print pdf
		a = []
		f = interpolate.interp1d(centers, row)
		for jj in range(1000):
			kk = random_distr(f,zmin,zmax)
			
			if not isinstance(kk,str):
				print kk
				a.append(kk-z)
		print len(a)
#		plt.plot(centers,row)
		plt.hist(a,normed=1)
		plt.savefig('a.png')
		if index == 0:
			exit()
		
	exit()
	y_pdf = bh.stackpdfs(y_pdf)
	
	sums = np.sum(y_pdf)*(centers[1]-centers[0])

	y_pdf = y_pdf/sums


	pdf_results.append({'name':lab,'df':y_pdf,'centers':centers})

	res = bh.weighted_nz_distributions(df_orig, centers,
		weights=weight_s, 
		tomo_bins = bin_edge, 
		z_phot = np.array(df_orig[z_pht]), 
		n_resample=10)
	x = res['binning']
	spec = res['spec']
 	phot = res['phot']
	Delta_pdf = []
	pdf_photz = []
	pdf_specz = []
	pdf_photz_shifted = []
	pdf_spec_plot = []
	delta_pdf_shift = []
	for i in range(len(phot)):
		
		Delta_pdf.append(phot[i][0]-spec[i][0])

		phot_shift = bh.eval_pdf_point(phot[i][0], x-res['div_means'][i][0], x_plot)
		spec_shift = bh.eval_pdf_point(spec[i][0], x, x_plot)

		pdf_photz.append(phot[i][0])
		pdf_specz.append(spec[i][0])
		pdf_photz_shifted.append(phot_shift)
		pdf_spec_plot.append(spec_shift)
		delta_pdf_shift.append(phot_shift-spec_shift)

		
    	phot_means = res['phot_means']
    	spec_means = res['spec_means']
    	div_means = res['div_means']
		

	results_stats.append({'code':lab,'phot_means':phot_means,'spec_means':spec_means,'div_means':div_means,'zbins':bin_center,'x':x,'delta_pdf':Delta_pdf,'pdf_photz':pdf_photz, 'pdf_specz':pdf_specz,'x_plot':x_plot,'pdf_photz_shifted':pdf_photz_shifted,'pdf_spec_plot':pdf_spec_plot,'delta_pdf_shift':delta_pdf_shift})
	#results_stats.append({'code':lab,'phot_means':phot_means,'spec_means':spec_means,'div_means':div_means,'zbins':bin_center,'x':x,'delta_pdf':Delta_pdf,'pdf_photz':pdf_photz, 'pdf_specz':pdf_specz})

cla()
clf()



for i in range(dim_lin):
	fig, ax = plt.subplots()

	ax.annotate('z = %.3f' % res['spec_means'][i+1][0], (0.4, 0.9), textcoords='axes fraction', size=10)

	va = [[] for _ in xrange(len(x_plot))]

	mean_phot = []
	for j,res in enumerate(results_stats):
		x = res['x_plot']
		delta_pdf = res['delta_pdf_shift']
		print len(delta_pdf[i+1])
		print len(x)
		plt.plot(x-res['spec_means'][i+1][0],delta_pdf[i+1],label=res['code'])

		print res['spec_means'][i+1]

		etaeta = bh.eval_pdf_point(delta_pdf[i+1], x-res['spec_means'][i+1][0], x_plot-res['spec_means'][i+1][0])

		for j,dp in enumerate(etaeta):
			va[j].append(dp)
		mean_phot.append(res['spec_means'][i+1][0])

	mean_phot = np.mean(mean_phot)

	yyyy = []
	errrrr = []
	for yyy in va:
	
		yyyy.append(np.mean(yyy))
		errrrr.append(np.std(yyy))

	plt.errorbar(x_plot-mean_phot, yyyy, yerr=errrrr,fmt='-o',label='mean')
	plt.legend()
	plt.xlabel(r'$\Delta_z (z-z_{mean})$')
	plt.ylabel(r'$(dn/dz)_{phot}-(dn/dz)_{spec}$')
	plt.title(r'Difference between stacked pdfs (photo-z - spec-z) at %s < z < %s' % (str(bin_edge[i]),str(bin_edge[i+1])))
	
	plt.axvline(x=0.0, color='k')
	plt.xlim((-1,1))
	plt.ylim((-1.8,1.8))
	savefig('test_%s_.png' % str(i))
	cla()
	clf()
	
	fig, ax = plt.subplots()

	ax.annotate('z = %.3f' % res['spec_means'][i+1][0], (0.4, 0.9), textcoords='axes fraction', size=10)

	plt.plot(x_plot-mean_phot,errrrr,'o')
	plt.xlabel(r'$\Delta_z (z-z_{mean})$')
	plt.ylabel(r'$\sigma_{(dn/dz)_{phot}-(dn/dz)_{spec}}$')
	plt.title(r'Error in the difference stacked pdfs (photo-z - spec-z) at %s < z < %s' % (str(bin_edge[i]),str(bin_edge[i+1])))
	
	plt.axvline(x=0.0, color='k')
	plt.xlim((-1,1))
	plt.ylim((0.,1.))

	savefig('error_%s_.png' % str(i))
	cla()
        clf()
	fig, ax = plt.subplots()

	ax.annotate('z = %.3f' % res['spec_means'][i+1][0], (0.4, 0.9), textcoords='axes fraction', size=10)

	plt.errorbar(x_plot-mean_phot, yyyy, yerr=errrrr,fmt='-o',label='mean')
	plt.xlabel(r'$\Delta_z (z-z_{mean})$')
	plt.ylabel(r'$(dn/dz)_{phot}-(dn/dz)_{spec}$')
	plt.title(r'Difference between stacked pdfs (photo-z - spec-z) at %s < z < %s' % (str(bin_edge[i]),str(bin_edge[i+1])))
	
	plt.axvline(x=0.0, color='k')
	plt.xlim((-1,1))
	plt.ylim((-1.8,1.8))

	savefig('dif_%s_.png' % str(i))
	cla()
        clf()

