import bpz.bpz_utils as bpz_utils
import logging

# Create a logger for all functions
LOGGER = bpz_utils.create_logger(level=logging.NOTSET,name='BPZ')


QUERY_SINGLE = """
SELECT
   coadd_object_id,
   tilename,
   RA,DEC,
   CM_FLUX_CORRECTED_G,
   CM_FLUXERR_CORRECTED_G,
   CM_FLUX_CORRECTED_R,
   CM_FLUXERR_CORRECTED_R,
   CM_FLUX_CORRECTED_I,
   CM_FLUXERR_CORRECTED_I,
   CM_FLUX_CORRECTED_Z,
   CM_FLUXERR_CORRECTED_Z,
   CM_MAG_CORRECTED_I
 FROM {tablename}, GTT_STR
 WHERE tilename=GTT_STR.STR
 ORDER BY coadd_object_id
"""

QUERY_MOF = """
SELECT
   coadd_object_id,
   tilename,
   cm_flux_g, 
   cm_flux_g/cm_flux_s2n_g as cm_fluxerr_g,
   cm_flux_r,
   cm_flux_r/cm_flux_s2n_r as cm_fluxerr_r,
   cm_flux_i,
   cm_flux_i/cm_flux_s2n_i as cm_fluxerr_i,
   cm_flux_z,
   cm_flux_z/cm_flux_s2n_z as cm_fluxerr_z
 FROM {tablename}, GTT_STR
 WHERE tilename=GTT_STR.STR
 ORDER BY coadd_object_id
"""

QUERY_SEX = """
select
  COADD_OBJECT_ID,
  EBV_SFD98,
  TILENAME,
  ALPHAWIN_J2000,
  DELTAWIN_J2000,

  FLUX_AUTO_G,
  FLUX_AUTO_R,
  FLUX_AUTO_I,
  FLUX_AUTO_Z,
  FLUX_AUTO_Y,
  FLUXERR_AUTO_G,
  FLUXERR_AUTO_R,
  FLUXERR_AUTO_I,
  FLUXERR_AUTO_Z,
  FLUXERR_AUTO_Y,

  MAG_AUTO_G,
  MAG_AUTO_R,
  MAG_AUTO_I,
  MAG_AUTO_Z,
  MAG_AUTO_Y,
  MAGERR_AUTO_G,
  MAGERR_AUTO_R,
  MAGERR_AUTO_I,
  MAGERR_AUTO_Z,
  MAGERR_AUTO_Y

  FROM {tablename}, GTT_STR
  where TILENAME=GTT_STR.STR
  ORDER BY coadd_object_id

"""


def get_dbh(db_section='db-desoper',verb=False):
    """ Get a DB handle"""
    from despydb import desdbi
    if verb: LOGGER.info("Creating db-handle to section: %s" % db_section)
    dbh = desdbi.DesDbi(section=db_section)
    return dbh

def get_root_archive(dbh, archive_name='desar2home',verb=False):
    """ Get the root-archive fron the database"""
    cur = dbh.cursor()
    # Get root_archive
    query = "select root from ops_archive where name='%s'" % archive_name
    if verb:
        logger.info("Getting the archive root name for section: %s" % archive_name)
        logger.info("Will execute the SQL query:\n********\n** %s\n********" % query)
    cur.execute(query)
    root_archive = cur.fetchone()[0]
    if verb: LOGGER.info("root_archive: %s" % root_archive)
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

