
from bpz_tools import *

from bpz_tools import prior
#don't even ask what's going on inside. It's a mess. Use as a black box.

import os
file_path = '/' + '/'.join([i for i in os.path.realpath(__file__).split('/')[0:-1]])

class BPZ_PRIOR():

    def __init__(self, z=np.arange(0, 2.5, 0.01), tipo_prior='cosmos_Laigle', ninterp=8, SED_LIST=file_path + '/SED/CWWSB4_6.list'):
        self.nt0 = len([i for i in open(SED_LIST, 'r')])
        self.z = z
        self.ninterp = ninterp
        self.tipo_prior = tipo_prior
        
    def evaluate(self, I_MAG):
        """evaluate prior at a particular magnitude"""
        """return, prior, redshift range"""
        p_i = prior(self.z, I_MAG, self.tipo_prior, self.nt0, self.ninterp)
        p_i = np.sum(p_i, axis=1)
        p_i = p_i / np.sum(p_i)
        return p_i, self.z

if __name__ == '__main__':
    c = BPZ_PRIOR()
    p, z = c.evaluate(20.5)
    print ('at Imag=', 20.5)
    print (z, p)