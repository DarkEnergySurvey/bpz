import despyastro
import numpy
import os
import sys
import subprocess
import time
import bpz.bpz_utils as bpz_utils

"""DB function helpers"""

# -------------------
# QUERY STRINGS
# -------------------

QUERY_COADD_CATS = """
SELECT
    fai.path as path,
    fai.filename as filename,
    fai.compression as compression,
    m.band as band
  FROM prod.proctag t, prod.catalog m, prod.file_archive_info fai
  WHERE t.tag='{tagname}'
        and t.pfw_attempt_id=m.pfw_attempt_id
        and m.tilename='{tilename}'
        and m.filetype='coadd_cat'
        and fai.filename=m.filename
        and fai.archive_name='desar2home' 
"""

def cmdline():

    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Get input files for DES Archive")
                                     # Inherit options from config_parser
                                     #parents=[conf_parser])
    parser.add_argument("tilename", action="store",default=None,
                        help="Name of Tilename to get")
    parser.add_argument("--tagname", action="store",default='Y3A1_COADD',
                        help="Tag Name")
    parser.add_argument("--section", action="store",default='db-desoper',
                        help="DB section to connect")
    parser.add_argument("--outpath", action="store",default=None,
                        help="Location of output files")
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Verbose?")
    parser.add_argument("--clobber", action="store_true", default=False,
                        help="Clobber files")
    args = parser.parse_args()

    if not args.outpath:
        args.outpath = "%s_bpz" % args.tilename

    return args


def get_dbh(db_section='db-desoper',verb=False):
    """ Get a DB handle"""
    from despydb import desdbi
    if verb: print "# Creating db-handle to section: %s" % db_section
    dbh = desdbi.DesDbi(section=db_section)
    return dbh

def get_root_archive(dbh, archive_name='desar2home',verb=False):
    """ Get the root-archive fron the database"""
    cur = dbh.cursor()
    # Get root_archive
    query = "select root from ops_archive where name='%s'" % archive_name
    if verb:
        print "# Getting the archive root name for section: %s" % archive_name
        print "# Will execute the SQL query:\n********\n** %s\n********" % query
    cur.execute(query)
    root_archive = cur.fetchone()[0]
    if verb: print "# root_archive: %s" % root_archive
    return root_archive

def get_root_https(dbh, archive_name='desar2home', logger=None):
    """ Get the root_https fron the database
    """
    if archive_name == 'desar2home':
        root_https = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive"
        return root_https

    cur = dbh.cursor()
    query = "SELECT val FROM ops_archive_val WHERE name='%s' AND key='root_https'" % archive_name
    if logger:
        logger.debug("Getting root_https for section: %s" % archive_name)
        logger.debug("Will execute the SQL query:\n********\n** %s\n********" % query)
    cur.execute(query)
    root_https = cur.fetchone()[0]
    if logger: logger.info("root_https:   %s" % root_https)
    cur.close()
    return root_https


def get_coadd_cats_from_db(dbh, tagname='Y3A1_COADD',**kwargs):
    """
    Execute database query to get coadd catalog files from the DB
    """

    tilename = kwargs.get('tilename')
    outpath  = kwargs.get('outpath')
    verbose  = kwargs.get('verbose')
    
    # Format and get the cat query
    if verbose:
        print "# Finding coadd catalogs for tilename:%s" % tilename
    query_coadd_cats = QUERY_COADD_CATS.format(tagname=tagname, **kwargs)
    cats = despyastro.query2dict_of_columns(query_coadd_cats, dbhandle=dbh)
    root_https   = get_root_https(dbh)

    # Here we fix 'COMPRESSION from None --> '' if present
    cats = fix_DB_None(cats)
    Nimages = len(cats['FILENAME'])
    cats['FILEPATH'] = [os.path.join(cats['PATH'][k],cats['FILENAME'][k]+cats['COMPRESSION'][k]) for k in range(Nimages)]
    cats['FILEPATH_HTTPS'] = [os.path.join(root_https,filepath) for filepath in cats['FILEPATH']]
    cats['FILEPATH_LOCAL'] = [os.path.join(outpath,filepath) for filepath in cats['FILENAME']]
    return cats

def fix_DB_None(rec):
    """
    Fix None --> '' in recordarray or dict returned by DB queries
    """

    if type(rec) == numpy.ndarray:
        if 'COMPRESSION' in rec.dtype.names:
            compression = [ '' if c is None else c for c in rec['COMPRESSION'] ]
            rec['COMPRESSION'] = numpy.array(compression)
    elif type(rec) == dict:
        if 'COMPRESSION' in rec.keys():
            rec['COMPRESSION'] = [ '' if c is None else c for c in rec['COMPRESSION'] ]
    return rec



def transfer_input_files(infodict, clobber=False, section='db-desoper', verbose=True):

    from despymisc import http_requests
    import tempfile
    
    """ Transfer the files contained in an info dictionary"""

    # Figure if we want to use wget or http_request
    url = "https://desar2.cosmology.illinois.edu/DESFiles/desarchive/ACT/check_connection.txt"
    localfile = os.path.join(tempfile.gettempdir(),"connect.txt")
    # Get a file using the $HOME/.desservices.ini credentials
    try:
        http_requests.download_file_des(url,localfile,section=section)
        use_wget = False
    except:
        print "WARNING: could not fetch file using Request class. Will try using old fashion wget now"
        use_wget = True
        
    # Now get the files via http
    Nfiles = len(infodict['FILEPATH_HTTPS'])
    for k in range(Nfiles):
           
        url       = infodict['FILEPATH_HTTPS'][k]
        localfile = infodict['FILEPATH_LOCAL'][k]
       
        # Make sure the file does not already exists exits
        if not os.path.exists(localfile) or clobber:
           
            dirname   = os.path.dirname(localfile)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
               
            sys.stdout.write("\rGetting:  %s (%s/%s)" % (url,k+1,Nfiles))
            sys.stdout.flush()

            if not use_wget:
                http_requests.download_file_des(url,localfile,section=section)
            else:
                status = get_file_des_wget(url,localfile,section=section,clobber=clobber)
                if status > 0:
                    raise RuntimeError("\n***\nERROR while fetching file: %s\n\n" % url)
        else:
            sys.stdout.write("\rSkipping: %s (%s/%s) -- file exists" % (url,k+1,Nfiles))
            sys.stdout.flush()

    print "\nDone"
       


def get_file_des_wget(url,localfile,section='http-desarchive',desfile=None,clobber=False):

    from despymisc import http_requests

    """
    A way to catch errors on http_requests.
    This whole fuction maybe you should go
    """

    # Read the credentials for the .desservices file
    USERNAME, PASSWORD, URLBASE = http_requests.get_credentials(desfile=desfile, section=section)
    WGET = "wget -q --user {user} --password {password} {url} -O {localfile}"
    kw = {'user':USERNAME, 'password':PASSWORD, 'url':url, 'localfile':localfile}
    cmd = WGET.format(**kw)
    if clobber and os.path.exists(localfile):
        os.remove(localfile)
   
    args = cmd.split()
    status = subprocess.call(args,env=os.environ.copy())
    return status



def main_transfer():

    t0 = time.time()
    args = cmdline()
    dbh = get_dbh(db_section=args.section,verb=args.verbose)

    # Make sure that outpath exists
    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)
    
    rec = get_coadd_cats_from_db(dbh,tilename=args.tilename,
                                 tagname=args.tagname,
                                 db_section=args.section,
                                 outpath=args.outpath,
                                 verbose=args.verbose)
    print "# Total DB Query time %s" % bpz_utils.elapsed_time(t0)
    if args.verbose: "# Will write files to %s" % args.outpath

    t1 = time.time()
    transfer_input_files(rec, clobber=args.clobber, section=args.section)
    print "# Total Transfer time %s" % bpz_utils.elapsed_time(t1)
    print "# Grand Total time %s" % bpz_utils.elapsed_time(t0)
    
