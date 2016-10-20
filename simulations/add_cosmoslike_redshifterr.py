import cPickle as p
import copy
import glob
from astropy.table import Table, vstack
import numpy as np

#add scatter to spec-z -> cosmos-z
i_err = pickle.load(open('data/CosmosAlhambra_phot_spec_z_rescaling.p', 'r'))
z_cos = np.zeros(len(d))

MAG_I = 'MAG_I'
TRUE_Z = 'Z'
for f in files:
    d = Table.read(f)
    for i in np.arange(len(i_err['i_mag'])-1):
        ind = np.array(d[MAG_I] > i_err['i_mag'][i]) * np.array(d[MAG_I] <= i_err['i_mag'][i+1])
        if np.sum(ind) > 0:
            z_tmp = np.array(d[TRUE_Z][ind]) + np.random.choice(i_err['cosmos']['{:}'.format(i_err['i_mag'][i])], size=np.sum(ind), replace=True)
            while np.sum(z_tmp < 0) > 0:
                ind1 = z_tmp < 0
                print i_err['i_mag'][i], np.sum(ind1)
                z_tmp[ind1] = np.array(d[TRUE_Z][ind])[ind1] + np.random.choice(i_err['cosmos']['{:}'.format(i_err['i_mag'][i])], size=np.sum(ind1), replace=True)
            z_cos[ind] = copy.copy(z_tmp)
    d['Z_COSMOS'] = z_cos
    d.write(f + '.zcosmos.fits')
