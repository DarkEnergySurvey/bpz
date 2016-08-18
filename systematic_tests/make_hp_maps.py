#! /usr/bin/env python

"""make pixelised maps of column properties, for each file at a user specified resolution"""
"""Author: Ben Hoyle """
"""v0.1 date: 13th May 2016"""
"""Requires, scipy, numpy, healpy, astropy"""

import math
import numpy as np
import sys
#leave double import please!
import numpy as numpy
import healpy as hp
from astropy.io import fits
from scipy.stats import binned_statistic


def radec_2_pix(ra, dec, nside=512, nest=False):
    """converts from ra,dec to IP pixel number"""
    theta, phi = (90. - dec) * math.pi / 180., ra * math.pi / 180.
    pixel_id = hp.ang2pix(nside, theta, phi, nest=nest)
    return pixel_id


def gen_map(ip_, arr, nside=512, statistic=np.mean):
    """this function calculates 'statistic' on 'arr' for all 'ip_' which sit in each healpix bin, given by nside"""
    # + 1 because binned_statistic is [2, 3)
    bins = np.arange(hp.nside2npix(nside) + 1, dtype=int)
    #this function calculates 'statistic' on arr for all ip_ which sit in each bin bins
    map_ = binned_statistic(ip_, arr, statistic=statistic, bins=bins)[0]
    return map_


def err_message(args):
    print "make_hp_maps.py mapFile[s].fits  [columns=Cols,To,Extract ra=RA dec=DEC z=MEAN_Z zbins=[0,0.1,0.2,0.3,0.5,0.7,0.9,1.1] nside=512 statistic=numpy.mean|len|numpy.std] "
    print args

if __name__ == "__main__":

    args = sys.argv[1:]
    if 'help' in args or len(args) < 1:
        err_message(args)
        sys.exit()

    #get inputs:
    files = []
    inputs = {}
    for i in args:
        if '=' not in i:
            files.append(i)
        else:
            ky, val = i.split('=')
            inputs[ky] = val

    if len(files) == 0:
        print "no included files"
        err_message(args)
        sys.exit()

    #which cols to extract (or all of them)
    if 'columns' in inputs:
        columns = inputs['columns'].split(',')

    #set ra , dec
    if 'ra' not in inputs:
        inputs['ra'] = 'RA'

    if 'dec' not in inputs:
        inputs['dec'] = 'DEC'

    #should we bin a column, if so, which column
    if 'bin_col' not in inputs:
        inputs['bin_col'] = 'MEAN_Z'

    #how should we calculate the 'pixel averaged statistic'
    if 'statistic' not in inputs:
        inputs['statistic'] = numpy.mean

    #how should we bin the requested column?
    if 'col_bins' in inputs:
        inputs['col_bins'] = eval(inputs['col_bins'])
        if 'bin_col' not in inputs:
            print '....bining on ' + inputs['bin_col'] + '....'

    #what resolution of the map do we care about?
    if 'nside' not in inputs:
        inputs['nside'] = 512

    # for each of the files
    for file_ in files:
        hdulist = fits.open(file_)

        #get columns if not inputted
        if 'columns' not in inputs:
            columns = hdulist[1].columns.names

        #get pixel id's of data
        pixel_id = radec_2_pix(hdulist[1].data[inputs['ra']], hdulist[1].data[inputs['dec']], nside=inputs['nside'])

        #make fits files for each column
        for col in columns:

            file_name = file_ + col + '.ring_map.fits'

            hp_map = gen_map(pixel_id, hdulist[1].data[col], nside=inputs['nside'])
            hp.write_map(file_name, hp_map)

            #cut columns by bins if requested
            if 'col_bins' in inputs and col == inputs['bin_col']:
                for i in range(len(inputs['col_bins']) - 1):

                    ind = (hdulist[1].data[col] < inputs['col_bins'][i + 1]) * (hdulist[1].data[col] > inputs['col_bins'][i])

                    if np.sum(ind) > 0:
                        file_name = file_ + col + '.' + str(inputs['col_bins'][i]) + '_' + str(inputs['col_bins'][i + 1]) + '.ring_map.fits'

                        hp_map = gen_map(pixel_id[ind], hdulist[1].data[col][ind], nside=inputs['nside'])
                        hp.write_map(file_name, hp_map)
