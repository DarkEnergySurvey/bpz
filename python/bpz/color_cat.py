"""
Utilities to assemble the input 'color' catalog needed by BPZ
"""

import time
import fitsio
import numpy
import bpz.bpz_utils as bpz_utils

def cmdline(parent_parser=None):

    import argparse
    import yaml

    parser = argparse.ArgumentParser(description="Build the BPZ input catalog",
                                     parents=[conf_parser])
    parser.add_argument("--incats", action="store",default=None,nargs='+',required=True,
                        help="Name of input fits catalog(s)")
    parser.add_argument("-c", "--config",help="BPZ yaml config file")
    parser.add_argument("--outcat", action="store",default=None,required=True,
                        help="Name of output 'color' fits catalog needed by BPZ")
    args = parser.parse_args()

    # Read in the full yaml config file, we only care about the 'filters' section, we'll ignore the rest
    defaults = yaml.load(open(args.config))
    
    # Return the dictionary 'filters' from the yaml file, but keyed to 'band' instead of response name
    args.filters = {}
    for filter_name in defaults['filters'].keys():
        band = defaults['filters'][filter_name]['band']
        args.filters[band] = defaults['filters'][filter_name]

    # Get the band's list
    args.bands = [args.filters[fname]['band'] for fname in args.filters]
    return args

def write_colorcat(args,data_in):

    # Define dtypes and record array
    dtypes = [('NUMBER','i8')]
    for BAND in args.bands:
        dtypes.append(("%s_%s" % (args.filters[BAND]['MAG_OR_FLUX'],BAND),'f8'))
        dtypes.append(("%s_%s" % (args.filters[BAND]['ERR'],BAND),'f8'))

    nrows = len(data_in[args.bands[0]])
    data_out = numpy.zeros(nrows, dtype=dtypes)

    # Now we populate the rec array over all bands
    for BAND in args.bands:
        # Loop over all dtypes and insert
        for key in data_in[BAND].dtype.names:
            if key == 'NUMBER':
                outkey = key
            else:
                outkey = "%s_%s" % (key,BAND)
            data_out[outkey] = data_in[BAND][key] 

    # Now we write the file
    fitsio.write(args.outcat, data_out, extname='OBJECTS', clobber=True)
    return

def read_cats(args):
    
    # Read in all of the input catalogs, and extract only the columns we care
    number=False
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
        if not number:
            col_names.insert(0,'NUMBER')
            number = True

        data_in[BAND] = tab['OBJECTS'].read(columns=col_names)
        format = "%s,"* len(data_in[BAND].dtype.names)
        msg = "# Reading "+ format[:-1] % data_in[BAND].dtype.names
        msg = msg + " from: %s" % incat
        print msg

    return data_in


def build():

    t0 = time.time()
    args = cmdline()
    data_in = read_cats(args)
    write_colorcat(args,data_in)
    print "# Wrote ColorCat %s in %s" % (args.outcat,bpz_utils.elapsed_time(t0))
    
if __name__ == '__main__':
    build()
