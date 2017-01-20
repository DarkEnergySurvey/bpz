#! /usr/bin/env python
"""
Script to upload data from a file to DESDM

version: 0.1.0

Author: Ben Hoyle benhoyle1212@gmail.com

dependencies: numpy, astropy, easyaccess [And you can login from command line]

call like:

./uploadDataDES -h file=fileName table=tableName [,rows_per_chunk=3e7, supress_warning=0, /verbose=1

-h | help show the help message

tableName: name of the table in DES DM

verbose: 1 | True to print lots fo stuff to screen as we go

rows_per_chunk: how many rows to read in/ send to DESDM each setting.

supress_warning: This *must* be set to 1, this assumes you know that we are going to 
                  write over some files. Beware!

To do:
delete tableName if it exists? currently doesn't handle nicely. and is a pain in oracle.

E.g:
./uploadDataDES.py file=BH_ADA_STAR_GAL_GOLDv101.fits table=ADA_STAR_GAL_GOLDv101 supress_warning=1 verbose=1
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
    print "./uploadDataDES file=fileName table=temp1_bh.csv.fits  [,rows_per_chunk=3e7, supress_warning=0, /verbose=1]"
    print "BE WARNED! tablename.fits will be written to disk here!"
    print "I will not do anything until you add supress_warning=1"
    print "to input. You have been warned!"
    

if os.path.isfile(inArgs['table'] + '.fits'):
    print "FYI: this file will be destroyed!: " + inArgs['table'] + '.fits'

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

