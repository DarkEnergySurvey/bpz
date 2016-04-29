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
python plot_nz.py z_pht=MEAN_Z path=/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/train_1_noweight/ path_true=/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/sims_spectra_representative_valid.WL_LSS_FLAGS.fits weight_s=WL_WEIGHTS
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

if 'path_true' not in inArgs.keys():
    print 'missing path for true fits file point statistics'
    print 'python plot_nz.py z_pht=MEAN_Z path=/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/train_1_noweight/valid_1_noweight/ path_true=sims_spectra_representative_valid.WL_LSS_FLAGS.fits [z_true=Z_TRUTH weight_s=WL_WEIGHTS]'
    sys.exit('ERROR: path for reference catalog missing')


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

path = inArgs['path']#'/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/train_1_noweight/valid_1_noweight/'
truefile = inArgs['path_true']
point_results = []
for i in glob.glob(path + '*_result_des.fits'):
	point_results.append(i)


#truefile = 'sims_spectra_representative_valid.WL_LSS_FLAGS.fits'
truedata = pf.open(truefile)[1]
z = truedata.data.field(z_true)

llist = []
truey,binEdges=histogram(z,bins=50,range=(0,2),normed=True)

bincenters = 0.5*(binEdges[1:]+binEdges[:-1])



plt.fill_between(bincenters, 0, truey, facecolor='gray')
labe = []
data_estimations = []

print np.sum(truey)


for res in point_results:

    data_1 = pf.open(res)[1].data
    phz = data_1[z_pht]
    tz = data_1['Z_SPEC']
    
    density = gaussian_kde(phz)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    ys = density(bincenters)

   

    llist.append(ys)
    lab = res.split('/')[-1].split('_')[0]
    print lab
    labe.append(lab)
    

    temp = plt.plot(bincenters, ys, antialiased=True, linewidth=2,label=lab)
    
    cc = temp[0].get_color()

    data_estimations.append({'code':lab,'x':bincenters,'y':ys,'color':cc})
#plotLib.plothistside(histo, hids, str(filter), '#',
#                'nc_'+str(filter)+'.png', mylog=True, step = False, showavg = True, lines=[mean_nc_maglim[filter]])

'''
plt.plot(bincenters,mean( array(llist), axis=0 ),'k--',lw=3,label='AVG',antialiased=True)
legend()
plt.ylabel('Density')
plt.xlabel('Redshift')
axes = plt.gca()
ylim = axes.get_ylim()
xlim = axes.get_xlim()
savefig('nz_sim.png')
cla()
clf()
plt.fill_between(bincenters, 0, truey, facecolor='gray',antialiased=True)
err = std( array(llist), axis=0 )
me = mean( array(llist), axis=0 )
plt.plot(bincenters,me,'k--',antialiased=True)
plt.fill_between(bincenters, me-err, me+err,antialiased=True)
plt.ylabel('Density')
plt.xlabel('Redshift')
axes = plt.gca()
axes.set_ylim(ylim)
savefig('nz_sim_average.png')
cla()
clf()
plot(bincenters,err)
plt.ylabel('RMS')
plt.xlabel('Redshift')

savefig(path+'error_nz_sim.png')
cla()
clf()
'''

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

x_plot = np.linspace(0., 2., 100.)
dim_algorithm = len(glob.glob(path + '*.hdf5'))
for i in glob.glob(path + '*.hdf5'):
	lab = i.split('/')[-1].split('_')[0]

	store = pd.HDFStore(i)
	df_orig = store['pdf']
	df = store['pdf']

	weight = df_orig[weight_s]
	store.close()
	centers = np.array([float(name[4:]) for name in df.columns if 'pdf' in name])

	y_pdf = np.array([name[:] for name in df.columns if 'pdf' in name])
	df = df[y_pdf]

	y_pdf = bh.normalisepdfs(df,centers)

	y_pdf = bh.stackpdfs(y_pdf)
	
	sums = np.sum(y_pdf)*(centers[1]-centers[0])

	y_pdf = y_pdf/sums

	cc = [d["color"] for d in data_estimations if d['code'] == lab][0]
	
#	plt.plot(centers, y_pdf, antialiased=True, linewidth=2,label=lab,color=cc)

	pdf_results.append({'name':lab,'df':y_pdf,'centers':centers})

	res = bh.weighted_nz_distributions(df_orig, centers,
		weights=weight_s, 
		tomo_bins = bin_edge, 
		z_phot = np.array(df_orig[z_pht]), 
		n_resample=20)
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


print len(x_plot)
'''
for i in range(dim_lin):
	va = [[] for _ in xrange(len(x_plot))]

	mean_phot = []
	for res in results_stats:
		x = res['x']
		delta_pdf = res['delta_pdf']
		plt.plot(x-res['phot_means'][i][0],delta_pdf[i])
		print res['phot_means'][i]
		etaeta = bh.eval_pdf_point(delta_pdf[i], x-res['phot_means'][i][0], x_plot-res['phot_means'][i][0])
		for j,dp in enumerate(etaeta):
			va[j].append(dp)
		mean_phot.append(res['phot_means'][i][0])

	mean_phot = np.mean(mean_phot)

	yyyy = []
	errrrr = []
	for yyy in va:
	
		yyyy.append(np.mean(yyy))
		errrrr.append(np.std(yyy))

	plt.errorbar(x_plot-mean_phot, yyyy, yerr=errrrr,fmt='-o')
#	plt.plot(x,yyyy,'o')

	plt.xlim((-1,1))
	savefig('test_%s.png' % str(i))
	cla()
	clf()
	plt.plot(x_plot-mean_phot,errrrr,'o')
	plt.xlim((-1,1))
	savefig('error_%s.png' % str(i))
	cla()
        clf()
	plt.plot(x_plot-mean_phot,yyyy,'o')
	plt.xlim((-1,1))
	savefig('dif_%s.png' % str(i))
	cla()
        clf()
'''
for i in range(dim_lin):
	va = [[] for _ in xrange(len(x_plot))]

	mean_phot = []
	for res in results_stats:
		x = res['x_plot']
		delta_pdf = res['delta_pdf_shift']
		print len(delta_pdf[i])
		print len(x)
		plt.plot(x-res['spec_means'][i][0],delta_pdf[i])

		print res['spec_means'][i]

		etaeta = bh.eval_pdf_point(delta_pdf[i], x-res['spec_means'][i][0], x_plot-res['spec_means'][i][0])

		for j,dp in enumerate(etaeta):
			va[j].append(dp)
		mean_phot.append(res['spec_means'][i][0])

	mean_phot = np.mean(mean_phot)

	yyyy = []
	errrrr = []
	for yyy in va:
	
		yyyy.append(np.mean(yyy))
		errrrr.append(np.std(yyy))

	plt.errorbar(x_plot-mean_phot, yyyy, yerr=errrrr,fmt='-o')
#	plt.plot(x,yyyy,'o')

	plt.xlim((-1,1))
	savefig('test_%s_.png' % str(i))
	cla()
	clf()
	plt.plot(x_plot-mean_phot,errrrr,'o')
	plt.xlim((-1,1))
	savefig('error_%s_.png' % str(i))
	cla()
        clf()
	plt.plot(x_plot-mean_phot,yyyy,'o')
	plt.xlim((-1,1))
	savefig('dif_%s_.png' % str(i))
	cla()
        clf()

exit()

'''
legend()
plt.ylabel('Density')
plt.xlabel('Redshift')
axes = plt.gca()
#ylim = axes.get_ylim()
axes = plt.gca()
axes.set_xlim(xlim)
axes.set_ylim(ylim)

plt.fill_between(bincenters, 0, truey, facecolor='gray',antialiased=True)
savefig('pdf_test_weight.png')

cla()
clf()
'''
print results_stats
phot_means = {}
spec_means = {}
div_means = {}
for i in range(dim_lin):
	phot_means[i] = []
	spec_means[i] = []
	div_means[i] = []

for res in results_stats:
	for i in range(dim_lin):
		print res['phot_means'][i]
		phot_means[i].append(res['phot_means'][i])
#		print res['spec_means'][i]
		spec_means[i].append(res['spec_means'][i])
#		print res['div_means'][i]
		div_means[i].append(res['div_means'][i])
		print

x = []
y = []
for bc in bin_center:
	x = x + [bc]*len(div_means[0])
for pm in div_means:
	for mp in div_means[pm]:
		print mp
		y.append(mp[0])

print y
y_mean, y_err = [], []
for i in range(dim_lin):
#	xx = []
	yy = []
	for j in range(dim_algorithm):
#		xx.append(x[j+i*5])
		yy.append(y[j+i*dim_algorithm])

	y_mean.append(np.mean(yy))
	y_err.append(np.std(yy))
plt.figure()
print bin_center, y_mean
plt.errorbar(bin_center, y_mean, yerr=y_err,fmt='-o')
plt.title("Difference between spec and photz using Y1 %s" % weight_s)
plt.ylabel(r'$\Delta_{photoz,spec}$')
plt.xlabel('Photo-z')
plt.savefig('div_mean_%s_%s.png' % (weight_s,str(dim_lin)))

plt.clf()
x = []
y = []
for bc in bin_center:
	x = x + [bc]*len(phot_means[0])
for pm in phot_means:
	for mp in phot_means[pm]:
		print mp
		y.append(mp[1])

print y
y_mean, y_err = [], []
for i in range(dim_lin):
#	xx = []
	yy = []
	for j in range(dim_algorithm):
#		xx.append(x[j+i*5])
		yy.append(y[j+i*dim_algorithm])

	y_mean.append(np.mean(yy))
	y_err.append(np.std(yy))
plt.figure()
print bin_center, y_mean
plt.errorbar(bin_center, y_mean, yerr=y_err,fmt='-o')
plt.title("Error in the mean using Y1 %s" % weight_s)
plt.ylabel(r'$\sigma_{mean}$')
plt.xlabel('Photo-z')
plt.savefig('err_mean_%s_%s.png' % (weight_s,str(dim_lin)))

plt.clf()
x = []
y = []
for bc in bin_center:
	x = x + [bc]*len(phot_means[0])
for pm in phot_means:
	for mp in phot_means[pm]:
		print mp
		y.append(mp[0])


print y
y_mean, y_err = [], []
for i in range(dim_lin):
#	xx = []
	yy = []
	for j in range(dim_algorithm):
#		xx.append(x[j+i*5])
		yy.append(y[j+i*dim_algorithm])

	y_mean.append(np.mean(yy))
	y_err.append(np.std(yy))
plt.figure()
print bin_center, y_mean
plt.errorbar(bin_center, y_mean, yerr=y_err,fmt='-o')
plt.title("Photoz mean using Y1 %s" % weight_s)
plt.ylabel(r'$photoz_{mean}$')
plt.xlabel('Photo-z')
plt.savefig('photoz_mean_%s_%s.png' % (weight_s,str(dim_lin)))

plt.clf()

fwhm = []

for i in glob.glob(path + '*.hdf5'):
	lab = i.split('/')[-1].split('_')[0]
	ff,topick, fff = bh.nz_test(i, lab, lab, save_plot=True, 
		weight_list = ['Y1_WEIGHTS'], 
		point_list =  [z_pht], 
		bin_list = [bin_edge],resample = 20 
		) 
	print 'FWHM'
	print '-----------'
	print fff
	print '----------'
	print
	fwhm.append(fff)

y = []
for ff in fwhm:
	yy = ff[0]	
        for rr in range(dim_lin):
                y.append(yy[rr+1])


y_mean, y_err = [], []

for i in range(dim_lin):
        yy = []
        for j in range(dim_algorithm):

                yy.append(y[j+i*dim_algorithm])

        y_mean.append(np.mean(yy))
        y_err.append(np.std(yy))
plt.figure()
print bin_center, y_mean
plt.errorbar(bin_center, y_mean, yerr=y_err,fmt='-o')
plt.title("pdf FWHM using Y1 %s" % weight_s)
plt.ylabel(r'fwhm')
plt.xlabel('Photo-z')
plt.savefig('fwhm_%s_%s.png' % (weight_s,str(dim_lin)))

plt.clf()

#print results_stats
print '----------------------------'
#print results_pdf


