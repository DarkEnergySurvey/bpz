#! /usr/bin/env python

def show_help(args):
    print ("""./identify_matching_columns.py
matches two files with the same column and then creates a matched CLASS 0 or 1 for use in a ML-type classifiucation problem.

inputs
file1=base_file.fits
file2=file_to_match.fits
file2_suffix=Suffix_For_Matched_column
out_file=matched_file_output.fits

optional inputs
verbose=1 -- show steps
match_column=COADD_OBJECTS_ID -- default COADD_OBJECTS_ID
join=all1|1and2|1or2|all1|all2|1not2|2not1|1xor2 -- default all1

e.g %> identify_matching_columns.py  file1=shuffled.NGMIX_MOF_Y1A1_WIDE_PHOT_001_RAND_4pcnt.LSS.WL.CLASS0_1.fits file2=test_LSS_SAMPLE_ID.v1.fits file2_suffix=LSS1 out_file=tet.lss.output.fits
""")
    print(inArgs)

import sys
args = sys.argv[1:]

inArgs = {}
for i in args:
    if '='in i:
        k, v = i.split('=')
        inArgs[k] = v

#protect the user, from themselves and us!
if 'out_file' not in inArgs:
    show_help(inArgs)
    sys.exit()
   
#print stuff to screen?
verbose = 'verbose' in inArgs
if verbose:
    verbose = inArgs['verbose']

if 'join' not in inArgs:
    inArgs['join'] = 'all1'
if 'match_column' not in inArgs:
    inArgs['match_column'] = 'COADD_OBJECTS_ID'

print inArgs

import os
import glob
from subprocess import call
from astropy.io import fits as pyfits
import numpy as np

cmd = "stilts tmatch2 matcher=exact in1={file1} in2={file2} values1={match_column} values2={match_column} suffix2={file2_suffix} join={join} out={out_file}".format(**inArgs)

cmd = cmd.split(" ")

call(cmd)

orig_table = pyfits.open(inArgs['out_file'])[1].data
orig_cols = orig_table.columns

CLASS = np.array(orig_table['COADD_OBJECTS_ID' + inArgs['file2_suffix']]>0).astype(int)

new_cols = pyfits.ColDefs([
                pyfits.Column(name=inArgs['file2_suffix'] + '_CLASS', format='D', array=CLASS)])

hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)

class_fname = inArgs['out_file'].split('.fits')[0] + '.CLASS.fits'
hdu.writeto(class_fname)

print ("file with classes {:} written to {:}".format(inArgs['file2_suffix'] + '_CLASS', class_fname))
