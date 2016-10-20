#!/usr/bin/env python
from __future__ import print_function, division
from glob import glob
from mpi4py import MPI
import numpy as np
import os
import sys
import yaml
import subprocess

if __name__ == '__main__':

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    cfgfile = sys.argv[1]
    with open(cfgfile, 'r') as fp:
        cfg = yaml.load(fp)

    catfiles = np.array(glob(cfg['FilePath']))

    if rank==0:
        try:
            os.mkdir(cfg['OPath'])
        except:
            pass

    for alg in cfg['Algorithms']:
        a = cfg['Algorithms'][alg]
        if alg=='BPZ':
            if rank==0:
                try:
                    os.mkdir("{0}/bpz/".format(cfg['OPath']))
                except:
                    pass

            for f in catfiles[rank::size]:
                fs = f.split('/')
                fsym = "{0}/bpz/{1}".format(cfg['OPath'], fs[-1])
                os.symlink(f, fsym)
                os.symlink(a['ColFile'], fsym[:-4]+'.columns')
                subprocess.call(['python', "{0}/bpz_des/bpz_will_mag.py".format(cfg['ExecPath']), fsym])
