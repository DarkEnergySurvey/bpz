#! /usr/bin/env python
import sys


def show_help(args):
    print ("""add_modest.py MOF|COADD Fits*.fits
        FitsFile. will add additional columns MODEST_CLASS and MODEST_CONT if MOF and
        MODES_CLASS +_1  MODEST_CLASS_CONT +_1 if COADD""")
    print (args)
    sys.exit()


def generate_modest_class(d, phot):

    if phot == 'MOF':
        MODEST_CONT = np.array(d['CM_T'] + d['CM_T_ERR'])
        d['MODEST_CONT'] = MODEST_CONT
        d['MODEST_CLASS'] = np.array(MODEST_CONT < 0.017, dtype=int)
    else:
        d['MODEST_CLASS'] = np.array(d['SPREAD_MODEL_I'] + (5.0/3.0) * d['SPREADERR_MODEL_I']) > 0.005
        d['MODEST_CLASS_1'] = np.abs(d['SPREAD_MODEL_I'] + (5.0/3.0) * d['SPREADERR_MODEL_I']) < 0.002
        d['MODEST_CONT'] = np.array(d['SPREAD_MODEL_I'] + (5.0/3.0) * d['SPREADERR_MODEL_I'])
        d['MODEST_CONT_1'] = np.abs(d['SPREAD_MODEL_I'] + (5.0/3.0) * d['SPREADERR_MODEL_I'])
    return d


def add_modest_save_file(ln):
    fil, phot = ln
    d = Table.read(fil)
    d = generate_modest_class(d, phot)
    d.write(fil.replace('.fits', '') + '.MODEST.fits')

if __name__ == "__main__":
    
    args = sys.argv[1:]

    if 'help' in args or len(args) < 2:
        show_help(args)

    phot = args[0]
    files = args[1:]

    from astropy.table import Table
    import numpy as np

    n_jobs = 5
    if n_jobs == 1:
        for fil in files:
            add_modest_save_file([fil, phot])
    else:

        lst = []
        for fil in files:
            lst.append([fil, phot])

        from bh_parallelise import Parrallelise

        res1 = Parrallelise(n_jobs=n_jobs, method=add_modest_save_file, loop=lst).run()
