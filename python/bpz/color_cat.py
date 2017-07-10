"""
Utilities to assemble the input 'color' catalog needed by BPZ
"""

import time
import fitsio
import numpy
import despyastro
import bpz.bpz_utils as bpz_utils
import bpz.db_utils as db_utils
import extinction 
import logging

# Create a logger for all functions
LOGGER = bpz_utils.create_logger(level=logging.NOTSET,name='BPZ')

def cmdline_coadd_files():

    import argparse
    import yaml

    # 1. We make a proto-parse use to read in the default yaml
    # configuration file, Turn off help, so we print all options in response to -h
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("-c", "--config",help="BPZ yaml config file")
    args, remaining_argv = conf_parser.parse_known_args()
    conf_defaults = yaml.load(open(args.config))

    # 2. This is the main parser
    parser = argparse.ArgumentParser(description="Build the BPZ input catalog",
                                     # Inherit options from config_parser
                                     parents=[conf_parser])
    parser.add_argument("--incats", action="store",default=None,nargs='+',required=True,
                        help="Name of input fits catalog(s)")
    parser.add_argument("--outcat", action="store",default=None,required=True,
                        help="Name of output 'color' fits catalog needed by BPZ")
    parser.add_argument("--PRIOR_MAGNITUDE", action="store",default=None,
                        help="Column Identifier for PRIOR_MAGNITUDE")
    parser.add_argument("--INPUT_MAGS",action="store_true", default=False,
                        help="Work with MAGS [True] or FLUX [False]")
    parser.add_argument("--ID",action="store", default=None,
                        help="Column name for ID (i.e. NUMBER)")
    parser.add_argument("--RA",action="store", default=None,
                        help="Column name for RA (i.e. ALPHAWIN_J2000)")
    parser.add_argument("--DEC",action="store", default=None,
                        help="Column name for Decl. (i.e. DELTAWIN_J2000)")
    parser.add_argument("--no-extinction",action="store_false", dest='extinction', 
                        help="Do not perform extiction correction")
    parser.add_argument("--extinction",action="store_true", dest='extinction', default=True,
                        help="Perform extiction correction")
    # Set the defaults of argparse using the values in the yaml config file
    parser.set_defaults(**conf_defaults)
    args = parser.parse_args()
    
    # Return the dictionary 'filters' from the yaml file, but keyed to 'band' instead of response name
    args.filters = {}
    for filter_name in conf_defaults['filters'].keys():
        band = conf_defaults['filters'][filter_name]['band']
        args.filters[band] = conf_defaults['filters'][filter_name]

    # The BAND -- not very elegant
    args.PRIOR_MAGNITUDE_BAND = args.PRIOR_MAGNITUDE.split('_')[-1]
    args.PRIOR_MAGNITUDE_NAME = args.PRIOR_MAGNITUDE.split('_'+args.PRIOR_MAGNITUDE_BAND)[0]

    # Get the band's list
    args.bands = [args.filters[fname]['band'] for fname in args.filters]

    #print "Will use:"
    #for k, v in vars(args).iteritems():
    #    print "%s: %s" % (k, v)

    return args

def extinction_correct(args,data_out):

    # get the E(B-V) for each position
    ra  = data_out[args.RA]
    dec = data_out[args.DEC]
    ebv,l,b = extinction.get_EBV_SFD98(ra,dec)

    for BAND in args.bands:
        key = "%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND)
        if args.INPUT_MAGS:
            data_out[key] = data_out[key] - args.efactor[BAND]*ebv
        else:
            dflux = 10**(0.4*args.efactor[BAND]*ebv)
            data_out[key] = data_out[key]*dflux
        LOGGER.info("Correcting %s" % key)
    return data_out


def write_colorcat(args,data_in):

    # Define dtypes and record array for ID, RA and DEC
    dtypes = [(args.ID,'i8'),
              (args.RA,'f4'),
              (args.DEC,'f4')]
        
    for BAND in args.bands:
        dtypes.append(("%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND),'f4'))
        dtypes.append(("%s_%s" % (args.filters[BAND]['ERR'],BAND),'f4'))

    # Add the prior magnitude, but make sure we do not duplicate
    colnames = [ colname for (colname, format) in dtypes]
    if args.PRIOR_MAGNITUDE not in colnames:
        dtypes.append((args.PRIOR_MAGNITUDE,'f4'))

    nrows = len(data_in[args.bands[0]])
    data_out = numpy.zeros(nrows, dtype=dtypes)

    # Now we populate the rec array over all bands
    for BAND in args.bands:
        # Loop over all dtypes and insert
        for key in data_in[BAND].dtype.names:
            if key in [args.ID, args.RA, args.DEC]:
                outkey = key
            else:
                outkey = "%s_%s" % (key,BAND)
            data_out[outkey] = data_in[BAND][key] 

    # Add the prior outside the loop
    if args.PRIOR_MAGNITUDE not in colnames:
        BAND = args.PRIOR_MAGNITUDE_BAND
        key  = args.PRIOR_MAGNITUDE_NAME
        outkey = "%s_%s" % (key,BAND)
        data_out[outkey] = data_in[BAND][key] 

    # If we want extinction corrected fluxes/mags
    if args.extinction:
        t0 = time.time()
        LOGGER.info("Computing E(B-V) for all objects")
        data_out = extinction_correct(args,data_out)
        LOGGER.info("E(B-V) Compute Time %s" % bpz_utils.elapsed_time(t0))
    else:
        LOGGER.info("Skipping extinction correction")

    # Now we write the file
    fitsio.write(args.outcat, data_out, extname='OBJECTS', clobber=True)
    return

def read_cats(args):
    
    # Read in all of the input catalogs, and extract only the columns we care
    id=False
    data_in = {}
    for incat in args.incats:
        
        tab =  fitsio.FITS(incat)
        header = tab[0].read_header()
        BAND = header['BAND'].strip()
        
        if BAND not in args.bands:
            print "%s not in bands -- ignoring" % BAND
            continue

        # Define the input colums per band
        col_names = [args.filters[BAND]['MAG_OR_FLUX'],
                     args.filters[BAND]['ERR']]
        if not id:
            col_names.insert(0,args.ID)  # i.e 'NUMBER'
            col_names.insert(1,args.RA)  # i.e.'ALPHAWIN_J2000'
            col_names.insert(2,args.DEC) # i.e.'DELTAWIN_J2000'
            id = True

        # Check if we need to add the PRIOR_MAGNITUDE to the read columns
        if BAND == args.PRIOR_MAGNITUDE_BAND and args.PRIOR_MAGNITUDE_NAME not in col_names:
            col_names.append(args.PRIOR_MAGNITUDE_NAME)
            
        data_in[BAND] = tab['OBJECTS'].read(columns=col_names)
        format = "%s,"* len(data_in[BAND].dtype.names)
        msg = "Reading "+ format[:-1] % data_in[BAND].dtype.names
        msg = msg + " from: %s" % incat
        LOGGER.info(msg)

    return data_in

def build_from_coadd_files():

    t0 = time.time()
    args = cmdline_coadd_files()
    data_in = read_cats(args)
    write_colorcat(args,data_in)
    LOGGER.info("Wrote ColorCat %s in %s" % (args.outcat,bpz_utils.elapsed_time(t0)))



# --- SQL Functions ------

def cmdline_sql():

    import argparse
    import yaml

    # 1. We make a proto-parse use to read in the default yaml
    # configuration file, Turn off help, so we print all options in response to -h
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("-c", "--config",help="BPZ yaml config file")
    args, remaining_argv = conf_parser.parse_known_args()
    conf_defaults = yaml.load(open(args.config))

    # 2. This is the main parser
    parser = argparse.ArgumentParser(description="Build the BPZ input catalog from SQL",
                                     # Inherit options from config_parser
                                     parents=[conf_parser])
    parser.add_argument("--tilename", action="store",default=None,required=True,
                        help="Name of Tilename to get")
    parser.add_argument("--outcat", action="store",default=None,required=True,
                        help="Name of output 'color' fits catalog needed by BPZ")
    parser.add_argument("--section", action="store",default='db-dessci',
                        help="DB section to connect")
    parser.add_argument("--PRIOR_MAGNITUDE", action="store",default=None,
                        help="Column Identifier for PRIOR_MAGNITUDE")
    parser.add_argument("--INPUT_MAGS",action="store_true", default=False,
                        help="Work with MAGS [True] or FLUX [False]")
    parser.add_argument("--ID",action="store", default=None,
                        help="Column name for ID (i.e. NUMBER)")
    parser.add_argument("--no-extinction",action="store_false", dest='extinction', 
                        help="Do not perform extiction correction")
    parser.add_argument("--extinction",action="store_true", dest='extinction', default=True,
                        help="Perform extiction correction")
    # Advanced options
    parser.add_argument("--QUERY_MOF", action="store",default=db_utils.QUERY_MOF,
                        help="Optional string with template for MOF SQL query")
    parser.add_argument("--QUERY_SEX", action="store",default=db_utils.QUERY_SEX,
                        help="Optional string with template for SExtractor SQL query")
    parser.add_argument("--PHOTO_TYPE", action="store",default='MOF_MIX',
                        help="Type of Photometry we want: MOF_MIX/MOF_ONLY/SEX_ONLY")
    # Set the defaults of argparse using the values in the yaml config file
    parser.set_defaults(**conf_defaults)
    args = parser.parse_args()
    
    # Return the dictionary 'filters' from the yaml file, but keyed to 'band' instead of response name
    args.filters = {}
    for filter_name in conf_defaults['filters'].keys():
        band = conf_defaults['filters'][filter_name]['band']
        args.filters[band] = conf_defaults['filters'][filter_name]

    # The BAND -- not very elegant
    args.PRIOR_MAGNITUDE_BAND = args.PRIOR_MAGNITUDE.split('_')[-1]
    args.PRIOR_MAGNITUDE_NAME = args.PRIOR_MAGNITUDE.split('_'+args.PRIOR_MAGNITUDE_BAND)[0]

    # Get the band's list
    args.bands = [args.filters[fname]['band'] for fname in args.filters]

    #print "Will use:"
    #for k, v in vars(args).iteritems():
    #    print "%s: %s" % (k, v)

    return args


def get_cat(args,query_type, dbh=None):
    """ Generic query catalog function """
    t0 = time.time()
    # Get db connection
    if not dbh:
        dbh = db_utils.get_dbh(db_section=args.section,verb=True)
    # string query
    if query_type == 'MOF':
        str_query = args.QUERY_MOF.format(tilename=args.tilename)
    elif query_type == 'SEX':
        str_query = args.QUERY_SEX.format(tilename=args.tilename)
    else:
        exit("ERROR: Query type %s not supported" % query_type)
    if args.verbose: LOGGER.info("Will execute %s query: %s" % (query_type, str_query))
    reccat = despyastro.query2rec(str_query, dbhandle=dbh)
    LOGGER.info("Total %s sql Time:%s" % (str_query,bpz_utils.elapsed_time(t0)))
    return reccat

def get_cats_from_sql(args):

    t0 = time.time()
    # Get db connection
    dbh = db_utils.get_dbh(db_section=args.section,verb=True)
    if args.PHOTO_TYPE=='MOF_MIX':
        mofcat = get_cat(args,'MOF',dbh=dbh)
        sexcat = get_cat(args,'SEX',dbh=dbh)
    elif args.PHOTO_TYPE=='MOF_ONLY':
        mofcat = get_cat(args,'MOF',dbh=dbh)
    elif args.PHOTO_TYPE=='SEX_ONLY':
        sexcat = get_cat(args,'SEX',dbh=dbh)
    else:
        exit("ERROR: PHOTO_TYPE:%s not defined" % args.PHOTO_TYPE) 
    LOGGER.info("Total sql Time:%s" % bpz_utils.elapsed_time(t0))

    MOF_FILTERS = [f.upper() for f in args.mof_filters]
    
    # Find the -9999 objects in MOF catalog, we can do it in any band
    FILTER = MOF_FILTERS[0]
    no_mof_idx = numpy.where(mofcat['CM_FLUX_%s' % FILTER] == -9999)

    # Here we pre-make the output record array.
    # It should contain all the output colums we want
    # Define dtypes and record array for ID, RA and DEC
    dtypes = [('COADD_OBJECT_ID','i8'),
              ('EBV_SFD98','f4')]

    for BAND in args.bands:
        dtypes.append(("%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND),'f4'))
        dtypes.append(("%s_%s" % (args.filters[BAND]['ERR'],BAND),'f4'))
    nrows = len(mofcat['COADD_OBJECT_ID'])
    data_out = numpy.zeros(nrows, dtype=dtypes)
    print dtypes

    # Case 1, we want a mix of Sextractor+mof
    for BAND in args.bands:
        FILTER = BAND.upper()
        print BAND, FILTER
        if BAND in args.mof_filters.keys():
            data_out["FLUX_MOFMIX_%s" % BAND]    = mofcat['CM_FLUX_%s' % FILTER]
            data_out["FLUXERR_MOFMIX_%s" % BAND] = mofcat['CM_FLUXERR_%s' % FILTER]
            # Replace the -9999 with SExtractor values
            data_out["FLUX_MOFMIX_%s" % BAND][no_mof_idx]    = sexcat['FLUX_AUTO_%s' % FILTER][no_mof_idx]
            data_out["FLUXERR_MOFMIX_%s" % BAND][no_mof_idx] = sexcat['FLUXERR_AUTO_%s' % FILTER][no_mof_idx]
        elif BAND in args.sex_filters.keys():
            data_out["FLUX_MOFMIX_%s" % BAND]    = sexcat['FLUX_AUTO_%s' % FILTER]
            data_out["FLUXERR_MOFMIX_%s" % BAND] = sexcat['FLUXERR_AUTO_%s' % FILTER]
        else:
            exit("ERROR: filters %s not in input filters" % BAND)

    # Now we write the file
    fitsio.write('mofmix.fits', data_out, extname='OBJECTS', clobber=True)

            
    #print mofcat['COADD_OBJECT_ID']
    #print sexcat['COADD_OBJECT_ID']
    return

def build_from_sql():

    t0 = time.time()
    args = cmdline_sql()
    get_cats_from_sql(args)

  
    
#if __name__ == '__main__':
#    build_from_coadd_files()


