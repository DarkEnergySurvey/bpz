#! /usr/bin/env python
"""
Script to retrieve hpix for a catalog ingested in DESDM

version: 0.1.0

Author: Aurelio Carnero aurelio.crosell@gmail.com

dependencies: easyaccess 

call like:

./gethpixDataDES.py table=tableName [output=tableName.hpix.csv ,nside=4096, verbose=1]

-h | help show the help message

tableName: name of the table in DES DM for which you want hpix

verbose: 1 | True to print lots fo stuff to screen as we go

nside: nside hpix to return back

output: csv file to export hpix

E.g:
./gethpixDataDES.py table=brportal.WL_REDUC_V2 output=wl_reduc_v2.hpix.csv nside=4096 verbose=1
"""
import sys
import os

args = sys.argv[1:]

inArgs = {}
for i in args:
    if '='in i:
        k, v = i.split('=')
        inArgs[k] = v

print args
#protect the user, from themselves and us!
if 'help' in args:
    print "./gethpixDataDES.py table=brportal.WL_REDUC_V2 [output=out.csv ,nside=4096 ,verbose=1]"
    print "BE WARNED! output table out.csv will be written to disk here!"
    print "You have been warned!"
    sys.exit(0)
if '-h' in args:
    print "./gethpixDataDES.py table=brportal.WL_REDUC_V2 [output=out.csv ,nside=4096 ,verbose=1]"
    print "BE WARNED! output table out.csv will be written to disk here!"
    print "You have been warned!"
    sys.exit(0)

if 'table' not in inArgs.keys():
	print 'missing table'
	print './gethpixDataDES.py table=brportal.WL_REDUC_V2 output=out.csv [,nside=4096 ,supress_warning=0 ,verbose=1]'
	sys.exit('ERROR: table name missing')
if 'nside' not in inArgs.keys():
	print 'nside default set to 4096'
	inArgs['nside'] = '4096'
if 'verbose' not in inArgs.keys():
	print 'verbose default to 1'
	inArgs['verbose'] = '1'

if 'output' not in inArgs.keys():
	print 'no output table name selected'
	print 'Output name will be %s.csv' % inArgs['table']
	inArgs['output'] = inArgs['table'] + '.csv'

if os.path.isfile(inArgs['output']):
    print "FYI: this file will be destroyed!: " + inArgs['output'] 

if inArgs['output'].endswith('.csv'):
	pass
else:
	print 'output table must end in csv'
	sys.exit('ERROR: output table must end in csv')

#print stuff to screen?
verbose = 'verbose' in inArgs
if verbose:
    verbose = inArgs['verbose']
if verbose:
    print inArgs

#load other stuff
import easyaccess as ea


#read table using astropy
table_name = inArgs['table']
nside = inArgs['nside'] 
#conenct to "easy-access"
connection = ea.connect()

##create a cursor object to handle the DB
cursor = connection.cursor() 

query = ("SELECT mcarras2.degrade(gold.HPIX, %s) as HPIX_%s, gold.coadd_objects_id "
	"FROM NSEVILLA.Y1A1_GOLD_1_0_2 gold, %s your where "
	"gold.coadd_objects_id = your.coadd_id" 
	% (nside, nside, table_name)
	) 
if verbose:
	print 'the query to save into %s' % inArgs['output']
	print query

if os.path.isfile(inArgs['output']):
	if verbose:
		print 'output exist, removing it...'
	os.remove(inArgs['output'])

if verbose:
	print "get hpix matching to gold table for ", table_name, "creating csv output called ", inArgs['output']

connection.query_and_save(query,inArgs['output'])

connection.close()

if verbose:
    print "Action completed. Please check in your directory. report bugs to Aurelio Carnero!"

