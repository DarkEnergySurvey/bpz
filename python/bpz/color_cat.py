"""
Utilities to assemble the input 'color' catalog needed by BPZ
"""

import os
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

def extinction_correct(args,data_out,ebv_key='EBV_SFD98'):

    # get the E(B-V) for each position

    # Check if we already have e(B-V) in the record array
    if ebv_key in data_out.dtype.names:
        ebv = data_out[ebv_key]
    else:
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

def build_color_cat_from_coadd_files():

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
    parser.add_argument("--tilename", action="store",nargs='+', default=None,required=True,
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
    parser.add_argument("--extinction",action="store_true", dest='extinction', default=False,
                        help="Perform extiction correction")
    parser.add_argument("--no-extinction",action="store_false", dest='extinction', default=True,
                        help="Do not perform extiction correction")
    # Advanced options
    parser.add_argument("--QUERY_MOF", action="store",default=db_utils.QUERY_MOF,
                        help="Optional string with template for MOF SQL query")
    parser.add_argument("--QUERY_SEX", action="store",default=db_utils.QUERY_SEX,
                        help="Optional string with template for SExtractor SQL query")
    parser.add_argument("--QUERY_SINGLE", action="store",default=db_utils.QUERY_SINGLE,
                        help="Optional string with template for SINGLE table SQL query")
    parser.add_argument("--PHOTO_MODE", action="store",default='MOF_MIX',
                        help="Type of Photometry we want: MOF_MIX/MOF_ONLY/SEX_ONLY")
    parser.add_argument("--table_SINGLE", action="store",default='MCARRAS2.Y3A2_MOF_CM_CORRECTED',
                        help="Modify the name of the table that contains the SExtractor COADD photometry")
    parser.add_argument("--table_MOF", action="store",default='y3a2_mof',
                        help="Modify the name of the table that contains the MOF photometry")
    parser.add_argument("--table_SEX", action="store",default='Y3A2_COADD_OBJECT_SUMMARY',
                        help="Modify the name of the table that contains the SExtractor COADD photometry")

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

    # Final sanity check on the QUERY strings in case empty in the conf file
    if not args.QUERY_MOF:
        args.QUERY_MOF = db_utils.QUERY_MOF
    if not args.QUERY_SEX:
        args.QUERY_SEX = db_utils.QUERY_SEX
    if not args.QUERY_SINGLE:
        args.QUERY_SINGLE = db_utils.QUERY_SINGLE

    if len(args.tilename) and os.path.exists(args.tilename[0]):
        args.tilelist = args.tilename[0]
        args.tilename = []
        for line in open(args.tilelist):
            args.tilename.append(line.split()[0])

    # Split them if they are comma separated and format for SQL
    args.tilename = parse_comma_separated_list(args.tilename)
    # Keep them as a list strcture for later
    args.tilename_list = args.tilename
    args.tilename = ["'%s'" % tilename for tilename in args.tilename]
    args.tilename = "(%s)" % ','.join(args.tilename)

    #print "Will use:"
    #for k, v in vars(args).iteritems():
    #    print "%s: %s" % (k, v)
    return args

def make_gtt_table(args,dbh=None):
    t0 = time.time()
    # Get db connection
    if not dbh:
        dbh = db_utils.get_dbh(db_section=args.section,verb=True)
    list = [ [tilename] for tilename in args.tilename_list]
    cur=dbh.cursor()
    LOGGER.info("Cleaning up GTT_STR")
    cur.execute('delete from GTT_STR')
    LOGGER.info("Inserting tilelist to GTT_STR.STR")
    dbh.insert_many('GTT_STR',['STR'],list)

def get_phot_catalog(args,query_type, dbh=None):
    """ Generic query catalog function """
    t0 = time.time()
    # Get db connection
    if not dbh:
        dbh = db_utils.get_dbh(db_section=args.section,verb=True)
    # string query
    if query_type == 'MOF':
        str_query = args.QUERY_MOF.format(tilename=args.tilename, tablename=args.table_MOF)
    elif query_type == 'SEX':
        str_query = args.QUERY_SEX.format(tilename=args.tilename, tablename=args.table_SEX)
    elif query_type == 'SINGLE':
        str_query = args.QUERY_SINGLE.format(tilename=args.tilename, tablename=args.table_SINGLE)
    else:
        exit("ERROR: Query type %s not supported" % query_type)
    if args.verbose: LOGGER.info("Will execute %s query: %s" % (query_type, str_query))
    reccat = despyastro.query2rec(str_query, dbhandle=dbh, verb=args.verbose)
    if args.verbose: LOGGER.info("Total query time: %s" % bpz_utils.elapsed_time(t0))

    if reccat is False:
        exit("ERROR: Query returned False")

    # Return the record-array catalog from the query
    return reccat

def get_MIX(args,data_in,data_out):

    """ Populate the data_out record array with the proper mix of MOF and SEXtractor quantities"""

    # Find the -9999 objects in MOF catalog, we can do it in any band
    BAND = args.mof_filters.keys()[0]
    no_mof_idx = numpy.where(data_in['MOF'][args.mof_filters[BAND]['MAG_OR_FLUX']] == -9999)

    # Fill in the mask
    data_out['MOF_MASK'][no_mof_idx] = 1

    # Case we want a mix of MOF+SExtractor
    for BAND in args.bands:
        # Short-cuts to outkey value/err
        outkey_val = "%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND)
        outkey_err = "%s_%s" % (args.filters[BAND]['ERR'],BAND) 
        sexkey_val = args.sex_filters[BAND]['MAG_OR_FLUX']
        sexkey_err = args.sex_filters[BAND]['ERR']
        if BAND in args.mof_filters.keys():
            mofkey_val = args.mof_filters[BAND]['MAG_OR_FLUX']
            mofkey_err = args.mof_filters[BAND]['ERR']
            data_out[outkey_val] = data_in['MOF'][mofkey_val]
            data_out[outkey_err] = data_in['MOF'][mofkey_err]
            # Replace the -9999 with SExtractor values
            data_out[outkey_val][no_mof_idx] = data_in['SEX'][sexkey_val][no_mof_idx]
            data_out[outkey_err][no_mof_idx] = data_in['SEX'][sexkey_err][no_mof_idx]
        elif BAND in args.sex_filters.keys():
            data_out[outkey_val] = data_in['SEX'][sexkey_val]
            data_out[outkey_err] = data_in['SEX'][sexkey_err]
        else:
            exit("ERROR: filter %s not in input filters" % BAND)
    return data_out

def get_ONLY(args,data_in, data_out):

    """ Get either the MOF or SExtractor catalog """

    # Define the inputs we want
    if args.PHOTO_MODE == 'MOF_ONLY':
        gen_filters = args.mof_filters
        gencat = data_in['MOF']
    elif args.PHOTO_MODE  == 'SEX_ONLY':
        gen_filters = args.sex_filters
        gencat = data_in['SEX']
    else:
        exit("ERROR: Query type: %s is not defined" % args.PHOTO_MODE)

    for BAND in args.bands:
        # Short-cuts to outkey value/err
        outkey_val = "%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND)
        outkey_err = "%s_%s" % (args.filters[BAND]['ERR'],BAND) 
        if BAND in gen_filters.keys():
            genkey_val = gen_filters[BAND]['MAG_OR_FLUX']
            genkey_err = gen_filters[BAND]['ERR']
            data_out[outkey_val] = gencat[genkey_val]
            data_out[outkey_err] = gencat[genkey_err]
        else:
            exit("ERROR: filter %s not in input filters" % BAND)

    return data_out

def get_SINGLE(args,data_in,data_out):

    """ Get the data_out catalog filled"""

    gencat = data_in[args.PHOTO_MODE]
    for BAND in args.bands:
        # Short-cuts to outkey value/err
        outkey_val = "%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND)
        outkey_err = "%s_%s" % (args.filters[BAND]['ERR'],BAND) 
        genkey_val = "%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND.upper())
        genkey_err = "%s_%s" % (args.filters[BAND]['ERR'],BAND.upper())
        data_out[outkey_val] = gencat[genkey_val]
        data_out[outkey_err] = gencat[genkey_err]

    return data_out

def read_catalogs_sql(args):

    t0 = time.time()
    # Get db connection
    dbh = db_utils.get_dbh(db_section=args.section,verb=True)

    # Make the GTT_STR table
    make_gtt_table(args,dbh=dbh)

    data_in = {}
    # Get the input catalogs using queries
    if args.PHOTO_MODE == 'SEX_ONLY':
        data_in['SEX'] = get_phot_catalog(args,'SEX',dbh=dbh)
    elif args.PHOTO_MODE == 'SINGLE':
        data_in['SINGLE'] = get_phot_catalog(args,'SINGLE',dbh=dbh)
    else:
        data_in['SEX'] = get_phot_catalog(args,'SEX',dbh=dbh)
        data_in['MOF'] = get_phot_catalog(args,'MOF',dbh=dbh)

        # Make sure the have the same lenght
        N_MOF = len(data_in['MOF']['COADD_OBJECT_ID'])
        N_SEX = len(data_in['SEX']['COADD_OBJECT_ID'])
        if N_MOF != N_SEX:
            LOGGER.info("ERROR: MOF and SExtractor queries returned different catalogs")
            LOGGER.info("Will attemp to examine the arrays before quitting")
            for tilename in args.tilename_list:
                id_MOF = numpy.where(data_in['MOF']['TILENAME'] == tilename)[0]
                id_SEX = numpy.where(data_in['SEX']['TILENAME'] == tilename)[0]
                if len(id_MOF) != len(id_SEX):
                    LOGGER.info("%s -- N(MOF)=%d, N(SEx)=%d" % (tilename, len(id_MOF), len(id_SEX)))
            exit("ERROR: MOF and SExtractor catalogs don't have the same number of objects")

        # Make sure that the two catalogs have the same COADD_OBJECT_ID
        equal_cats = (data_in['MOF']['COADD_OBJECT_ID']==data_in['SEX']['COADD_OBJECT_ID']).all()
        if not equal_cats:
            exit("ERROR: MOF and SExtractor catalogs don't have the same objects")
        
    LOGGER.info("Total sql Time:%s\n" % bpz_utils.elapsed_time(t0))
    return data_in

def write_colorcat_sql(args,data_in):

    # Here we pre-make the output record array.
    # It should contain all the output colums we want

    MODE=args.PHOTO_MODE

    if args.PHOTO_MODE == 'MOF_MIX' or args.PHOTO_MODE=='SEX':
        MODE='SEX'
    else:
        MODE=args.PHOTO_MODE

    # Figure out the dtypes for additional columns
    dtypes = []
    for colname in args.ADDITIONAL_OUTPUT_COLUMNS:
        dtype = ("%s" % data_in[MODE][colname].dtype)
        dtypes.append((colname,dtype))

    for BAND in args.bands:
        dtypes.append(("%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND),'f4'))
        dtypes.append(("%s_%s" % (args.filters[BAND]['ERR'],BAND),'f4'))
    # Add the prior magnitude, but make sure we do not duplicate
    colnames = [ colname for (colname, format) in dtypes]
    if args.PRIOR_MAGNITUDE not in colnames:
        dtypes.append((args.PRIOR_MAGNITUDE,'f4'))

    if args.PHOTO_MODE=='MOF_MIX':
        dtypes.append(('MOF_MASK','i4'))

    nrows = len(data_in[MODE][:])
    data_out = numpy.zeros(nrows, dtype=dtypes)

    # Populate the basic information
    for key in args.ADDITIONAL_OUTPUT_COLUMNS:
        data_out[key] = data_in[MODE][key]

    if args.PHOTO_MODE=='MOF_MIX':
        data_out = get_MIX(args,data_in,data_out)
    elif args.PHOTO_MODE=='SINGLE':
        data_out = get_SINGLE(args,data_in,data_out)
    else:
        data_out = get_ONLY(args,data_in,data_out)

    # Add the prior outside the loop
    if args.PRIOR_MAGNITUDE not in colnames:
        data_out[args.PRIOR_MAGNITUDE] = data_in[MODE][args.PRIOR_MAGNITUDE.upper()] 

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
    LOGGER.info("Wrote ColorCat: %s" % args.outcat)
    return 

def parse_comma_separated_list(inputlist):
    if inputlist[0].find(',') >= 0:
        return inputlist[0].split(',')
    else:
        return inputlist

def build_color_cat_from_sql():

    t0 = time.time()
    args = cmdline_sql()
    data_in = read_catalogs_sql(args)
    write_colorcat_sql(args,data_in)
    LOGGER.info("Total time: %s" % (bpz_utils.elapsed_time(t0)))

