import despyastro
import numpy

"""DB function helpers"""

# -------------------
# QUERY STRINGS
# -------------------

QUERY_COADD_CATS = """
SELECT
    fai.path as path,
    fai.filename as filename,
    fai.compression as compression,
    m.band as band,
    m.pfw_attempt_id as pfw_attempt_id
  FROM prod.proctag t, prod.catalog m, prod.file_archive_info fai
  WHERE t.tag='{tagname}'
        and t.pfw_attempt_id=m.pfw_attempt_id
        and m.tilename='{tilename}'
        and m.filetype='coadd_cat'
        and fai.filename=m.filename
        and fai.archive_name='desar2home' 
"""

def get_dbh(db_section='desoper',verb=False):
    """ Get a DB handle"""
    from despydb import desdbi
    if verb: print "# Creating db-handle to section: %s" % db_section
    dbh = desdbi.DesDbi(section=db_section)
    return dbh

def get_root_archive(dbh, archive_name='prodbeta',verb=False):
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

def get_coadd_cats_from_db(dbh, **kwargs):

    """
    Execute database query to get coadd catalog files from the DB
    """

    dbh = get_dbh(db_section='desoper',verb=True)
    query_coadd_cats = QUERY_COADD_CATS.format(**kwargs)
    cats = despyastro.query2rec(query_coadd_cats, dbhandle=dbh)
    print cats
    return cats

def transfer_input_files(infodict, clobber, section, logger=None):

    from despymisc import http_requests
   
    """ Transfer the files contained in an info dictionary"""
   
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
               
            logger.info("Getting:  %s (%s/%s)" % (url,k+1,Nfiles))
            sys.stdout.flush()

            try:
                # Get a file using the $HOME/.desservices.ini credentials
                http_requests.download_file_des(url,localfile,section=section)
            except:
                warning = """WARNING: could not fetch file: %s using Request class. Will try using old fashion wget now"""  % url
                logger.info(warning)
                status = get_file_des_wget(url,localfile,section=section,clobber=clobber)
                if status > 0:
                    raise RuntimeError("\n***\nERROR while fetching file: %s\n\n" % url)

        else:
            logger.info("Skipping: %s (%s/%s) -- file exists" % (url,k+1,Nfiles))
       


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
    if clobber:
        os.remove(localfile)
   
    status = work_subprocess(cmd)
    return status
