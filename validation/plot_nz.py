#import lib.plots.plotLib as plotLib
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import glob
import yaml 


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
    print 'python plot_nz.py z_pht=MEAN_Z path=/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/train_1_noweight/valid_1_noweight/ path_true=sims_spectra_representative_valid.WL_LSS_FLAGS.fits [z_true=Z_TRUTH]'
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

path = inArgs['path']#'/home/carnero/Dropbox/DES_photoz_wg/project38/nz_codes/train_1_noweight/valid_1_noweight/'
truefile = inArgs['path_true']
results = []
for i in glob.glob(path + '*.fits'):
	results.append(i)


#truefile = 'sims_spectra_representative_valid.WL_LSS_FLAGS.fits'

truedata = pf.open(truefile)[1]
z = truedata.data.field(z_true)

llist = []
truey,binEdges=histogram(z,bins=50,range=(0,2),normed=True)

bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#plt.plot(bincenters,truey,'-')

nor = max(truey)
#truey,binEdges=histogram(z,bins=50,range=(0,2),normed=nor)

plt.fill_between(bincenters, 0, truey, facecolor='gray')
labe = []
data_estimations = []
for res in results:

    data_1 = pf.open(res)[1].data
    phz = data_1[z_pht]
    tz = data_1['Z_SPEC']
    data_estimations.append(phz)
#    phz,tz = loadtxt(res,usecols=(0,2),unpack=True)
#    photoz.append(phz)
    y,binEdges=histogram(phz,bins=50,range=(0,2),normed=nor)
#    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    llist.append(y)
    lab = res.split('_')[-1].split('.fits')[0]
    labe.append(lab)
    
    plot(bincenters,y,'-',label=lab)
#    semilogy(bincenters,y,'-')
    

#plotLib.plothistside(histo, hids, str(filter), '#',
#                'nc_'+str(filter)+'.png', mylog=True, step = False, showavg = True, lines=[mean_nc_maglim[filter]])
plot(bincenters,mean( array(llist), axis=0 ),'k--',lw=3,label='AVG')
legend()
plt.ylabel('Density')
plt.xlabel('Redshift')

savefig(path+'nz_sim.png')
cla()
clf()
plt.fill_between(bincenters, 0, truey, facecolor='gray')
err = std( array(llist), axis=0 )
me = mean( array(llist), axis=0 )
plt.plot(bincenters,me,'k--')
plt.fill_between(bincenters, me-err, me+err)
plt.ylabel('Density')
plt.xlabel('Redshift')

savefig(path+'nz_sim_average.png')
cla()
clf()
plot(bincenters,err)
plt.ylabel('RMS')
plt.xlabel('Redshift')

savefig(path+'error_nz_sim.png')
cla()
clf()


p = load_yaml('./testConfig/photoz.yaml')

bin_edge = p['point']['bins'][0]['MEAN_Z'].split('[')[1].split(']')[0].split(',')
bin_edge = map(str.strip,bin_edge)
bin_edge = map(float,bin_edge)
print bin_edge

for ii in range(len(bin_edge)-1):
    zmin=bin_edge[ii]
    zmax=bin_edge[ii+1]
    mak = (tz> zmin)*(tz<=zmax)
    llist = []
    for jj,de in enumerate(data_estimations):
	temp_de=de[mak]
	y,binEdges=histogram(temp_de,bins=50,range=(0,2),normed=nor)

	llist.append(y)
	plot(bincenters,y,'-',label=labe[jj])
   
    
    plot(bincenters,mean( array(llist), axis=0 ),'k--',lw=3,label='AVG')

    legend()
    plt.ylabel('Density')
    plt.xlabel('Redshift')
    plt.axvline(x=zmin, linewidth=2,color='k',linestyle='--')
    plt.axvline(x=zmax, linewidth=2,color='k',linestyle='--')
    savefig(path+'z_nz_bins_%s_%s.png' % (str(zmin),str(zmax)))
    cla()
    clf()

    err = std( array(llist), axis=0 )
    me = mean( array(llist), axis=0 )
    plt.plot(bincenters,me,'k--')
    plt.fill_between(bincenters, me-err, me+err)
    plt.ylabel('Density')
    plt.xlabel('Redshift')
    plt.axvline(x=zmin, linewidth=2, color='k',linestyle='--')
    plt.axvline(x=zmax, linewidth=2, color='k',linestyle='--')

    savefig(path+'z_nz_avg_%s_%s.png' % (str(zmin),str(zmax)))
    cla()
    clf()


