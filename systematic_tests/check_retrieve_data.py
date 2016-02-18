#! /usr/bin/env python
import sys
import os

#some example data products that exist in DESDM
availableZtables = ['hoyleb.PHOTOZ_Y1_RF_V01', 'hoyleb.PHOTOZ_HWE_Y1_V0_1', 'hoyleb.PHOTOZ_ADA_Z_Y1_V0_3', 'hoyeb.PHOT0Z_Y1_BPZ_V01', 'hoyleb.PHOTOZ_DNF_Y1_V01']
availableSamples = ['hoyleb.LSS_COADD_OBJECTS_ID_1M', 'hoyleb.WL_COADD_OBJECTS_ID_1M', 'hoyleb.REDMAG6411HIGHLUM1004_1M', 'hoyleb.REDMAG6411HIGHDENS_1M', 'hoyleb.MAIN_LSS_Y1V01', 'hoyleb.RED_LSS_Y1V01']
avaiableData = ['hoyleb.IM3SHAPE_Y1V1', 'NSEVILLA.Y1A1_GOLD_1_0_2']


def show_help():
    """Show the user how to use this file"""
    print "./check_retrieve_data -h or help to see the help file"
    print "./check_retrieve_data redshift-table=" + '|'.join(availableZtables) + " sample-table=" + '|'.join(availableSamples) + " data-table=" + '|'.join(avaiableData) + " [,max-rows=5]"
    print "You may use any other tables, as long as you can see them in DESDM"
    print "We first check for the existence for the file, then download if it is not there."
    sys.exit()

args = sys.argv[1:]

#show help if asked
if ('-h' in args) or ('help' in args):
    show_help()

inArgs = {}
for i in args:
    if '='in i:
        k, v = i.split('=')
        inArgs[k] = v

if (inArgs == {}):
    show_help()

fileToUse = ''
for i in inArgs:
    if 'table' in i:
        fileToUse += inArgs[i].split('.')[-1] + '.'

fileToUse += 'fits'


if os.path.isfile(fileToUse):
    print "FYI: you already have this file. We do not need to download it."
    sys.exit()

#load other stuff
import easyaccess as ea
import numpy as np
from astropy.table import Table
import string
import random

#conenct to "easy-access"
connection = ea.connect()

##create a cursor object to handle the DB
cursor = connection.cursor()

inArgs = {}

#should we only select some data?
maxRows = ''
if 'max-rows' in inArgs:
    maxRows = ' where ROWNUM<' + inArgs['max-rows'] + 1

#construct the select query from the input tables.
prt1 = 'SELECT '
prt2 = ' '
tbl = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
tnum = 0

#play some encoding tricks because we are using crappy oracle
for i in inArgs:
    if 'table' in i:
        prt1 += tbl[tnum] + '.*,'
        prt2 += ' JOIN ' + inArgs[i] + ' ' + tbl[tnum]
        if tnum > 0:
            prt2 += ' ON ' + tbl[0] + '.COADD_OBJECTS_ID=' + tbl[tnum]+'.COADD_OBJECTS_ID '
        tnum += 1

prt1 = prt1[0:-1] + ' '
prt2 = prt2[6:] + ' '

#build the query to run on DESDM
query = prt1 + ' FROM ' + prt2 + ' ' + maxRows

print "preparing to run this query on DESDM:"
print query

#get the query results [into memory? be careful!]
d = connection.query_to_pandas(query, prefetch='')

#if we have repeated columns, e.g. COADD_OBJECTS_ID, add some random letters to the repeats

kys = []
for i in d.keys():
    if i not in kys:
        kys.append(i)
    else:
        rnd = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(3))
        kys.append(i + '_' + rnd)

#turn pandas data frame into a fits, because.. um.. fits.
t = Table(rows=np.array(d), names=kys)
t.write(fileToUse, format='fits')

#profit.
