#!/usr/bin/env python

import sys

def show_help(args):
    print ("""cut_data.py dataFile column numpy.greater|less_equal|greater_equal|equal|not_equal VALUE""")
    print(args)
    sys.exit()

args = sys.argv[1:]

if len(args) != 4:
    show_help(args)

import numpy as np
import numpy as numpy
from astropy.table import Table

def get_function(function_string):
    import importlib
    module, function = function_string.rsplit('.', 1)
    module = importlib.import_module(module)
    function = getattr(module, function)
    return function

print args, eval(args[3])

crit = get_function(args[2])
d = Table.read(args[0])
print ('before {:}'.format(args[1:]), len(d))

try: 
    print ('evaluating {:}'.format(args[3]))
    d = d[crit(d[args[1]], eval(args[3]))]
except:
    print ('not evaluating {:}'.format(args[3]))
    d = d[crit(d[args[1]], args[3])]

print ('after {:}'.format(len(d)))

d.write(args[0].replace('.fits','') + '.cut.fits')
