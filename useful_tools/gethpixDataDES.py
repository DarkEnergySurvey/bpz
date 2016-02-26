#! /usr/bin/env python
"""
Script to retrieve hpix for a catalog ingested in DESDM

version: 0.1.0

Author: Aurelio Carnero aurelio.crosell@gmail.com

dependencies: numpy, astropy, easyaccess [And you can login from command line]

call like:

./gethpixDataDES.py table=tableName output=tableName.hpix.csv [,nside=4096, supress_warning=1, verbose=1]

-h | help show the help message

tableName: name of the table in DES DM for which you want hpix

verbose: 1 | True to print lots fo stuff to screen as we go

nside: nside hpix to return back

output: csv file to export hpix

E.g:
./gethpixDataDES.py table=brportal.WL_REDUC_V2 output=wl_reduc_v2.hpix.csv nside=4096 supress_warning=1 verbose=1
"""
import sys
import os

args = sys.argv[1:]

inArgs = {}
for i in args:
    if '='in i:
        k, v = i.split('=')
        inArgs[k] = v

#protect the user, from themselves and us!
if '-h' or 'help' in args or 'supress_warning' not in inArgs:
    print "./gethpixDataDES.py table=brportal.WL_REDUC_V2 output=out.csv [,nside=4096 ,supress_warning=0 ,verbose=1]"
    print "BE WARNED! output table will be written to disk here!"
    print "I will not do anything until you add supress_warning=1"
    print "to input. You have been warned!"
    

if os.path.isfile(inArgs['output']):
    print "FYI: this file will be destroyed!: " + inArgs['table'] 

if inArgs['supress_warning'] == 0:
    print "Nothing happens until supress_warning=1"
    1/0

#print stuff to screen?
verbose = 'verbose' in inArgs
if verbose:
    verbose = inArgs['verbose']
if verbose:
    print inArgs

#load other stuff
from astropy.table import Table
import easyaccess as ea
import numpy as np


#read table using astropy
d = Table.read(inArgs['file'])


#decide, or be told, how many upload trips to make
if 'rows_per_chunk' not in inArgs:
    rows_per_chunk = int(4e7)
else:
    rows_per_chunk = int(inArgs['rows_per_chunk'])


#if the file is smaller than number rows
if rows_per_chunk > len(d):
    rows_per_chunk = len(d)


#determine number of chunks to split file into
n_chunks = int(len(d) / rows_per_chunk)
if verbose:
    print "file size", len(d)
    print "splitting file into", n_chunks," chunks"
    print "rows_per_chunk: ", rows_per_chunk


#conenct to "easy-access"
connection = ea.connect()

##create a cursor object to handle the DB
cursor = connection.cursor() 

#To do. this is a pain in ORACLE!

#check to see if table already exists:
#q = "select * from " + inArgs['table'] + " where rownum<2"
#res = cursor.execute(q).descrption



for i in range(n_chunks):
    #remove file if exists
    if os.path.isfile(inArgs['table']+'.fits'):
        os.remove(inArgs['table']+'.fits')

    #load each batch
    ind = np.arange(rows_per_chunk, dtype=int) + i * rows_per_chunk

    #finally get all remaining rows
    if i == n_chunks-1:
        ind = np.arange(len(d) - i* rows_per_chunk, dtype=int) + i * rows_per_chunk
    
    if verbose:
        print "this chunk", i, " total chunks", n_chunks
        print "n_rows this time", len(ind)
        print "min max index:", np.amin(ind), np.amax(ind)

    indr = np.array([False] * len(d))
    indr[ind] = True
    #write this chunk to disk
    d[indr].write(inArgs['table']+'.fits')

    #upload this chunk
    if i ==0:
        connection.load_table(inArgs['table']+'.fits', name=inArgs['table'])
    else:
        connection.append_table(inArgs['table']+'.fits', name=inArgs['table'])

if verbose:
    print "Action completed. Please check in DESDM. report bugs to Ben Hoyle!"

