This directory contains an updated version of BPZ.
Authors: Ben Hoyle

1) About

1.1) Requirements

2) Example

3) Output 

3.1) point prediction files
3.2) pdf prediction files

4) comparison with old BPZv1.99.3

5) Saving and recovering the best fitting templates, and synthetic fluxes

6) To do

================================
1) About
================================
--- What this code does: ---
This is a re-writing of BPZ. We make the code a lot cleaner, so that it could be parralised or cythonised in the future.

The prior functions has re-written the priors into a nice class that can be callable from anywhere, not just BPZ

It is also trivial to change the priors. Simply look at sed_prior_file.py for an example of how to change the priors.

--- What this code doesn't do:---
BPZ does lots of error checking under the hood and deals with missing / unseen data nicely. This is not currently implemented here.

BPZ allows you to understand if you photometry is offset, and how to fix this. This is also not implemented.

1.1) Requirements

import sys
from joblib import Parallel, delayed
import numpy as np
import inspect
import random as rdm
import yaml
import os
import copy
from astropy.io import fits as pyfits
import copy
import time
import pandas as pd

-in the photoz-wg /validation/ directory
import bh_photo_z_validation as pval

-- in this directory
from galaxy_type_prior import GALAXYTYPE_PRIOR

================================
2) Example usage. 
================================
Use a text editor to look at bpzConfig.yaml. There you can make all the adjustments, set directory paths that you used to make within the BPZ code.

To generate an exmaple yaml file run
%>./bpzv1.py 

-- this write exampleBPZConfig.yaml to disk (if it doesn't already exist) 

To run the code

%>./bpzv1.py PathToConfig.yaml PathToListofFitsFiles.fits

PathToListofFitsFiles.fits can contain a * so that all files will be processsed separately


where PathToConfig.yaml is a configuration file, that defines how to run the code. See bpzConfig.yaml for an example

================================
3) Output
================================
3.1) point predictions: the output file has the same name as the input fits file, with .BPZ.fits at the end.

The code generates a output fits files, which contain 'MEAN_Z'  'Z_SIGMA' 'MEDIAN_Z': 'Z_MC': 'Z_SIGMA68': and any additional columns e.g. REDSHIFT or MAG_I or COADDED_OBJECTS_ID that you asked for within the config.yaml file.

Additional outputs are: 
KL_POST_PRIOR - the information gain (Kullbeck Leibler divergence) between prior and posterior
chi2 - the minimum \Chi^2 to the closest template synthetic example
TEMPLATE_TYPE - the maximum posterior template type, as defined as a linear combination of input templates using a floating point notation [as per standard BPZ]! 
    #each sed (in order) is a number between 1 -> Num seds.
    #interpolatated sed are fractional quantites between, 
        e.g. TEMPLATE_TYPE = 1.4  meaning 0.6 of sed1 and 0.4 of sed2 
        e.g. TEMPLATE_TYPE = 3.1  meaning 0.9 of sed3 and 0.1 of sed4

TEMPLATE_ID - the SED ID that corresponds to the best fitting template.

3.2) pdf prediction files: the output file has the same name as the input fits file, with .h5 at the end

Various extensions:
/point_predictions/ table extension
    -- all info that is stored in .fits files, can be found here
    -- the ID specified in the config file

/info/ table extension
    -- z_bin_centers: bin centers (also written as the column names. This is an easy access)

/point_predictions/ table extension
    -- the pdfs labelled like "pdf_{:0.4}.format(bin_center)"
    -- the ID specified in the config file


================================
4) Comparison with older BPZ
================================

In the test/ directory there is .cat file and a .fits file of the same data. 

The .cat file works with older version of bpz

#compare the output of original BPZ and this version
cd ../bpz-1.99.3/
%> python bpz_tcorrv1.py ../bpzv1/test/WL_CLASS.METACAL.rescaled.slr.cosmos.v2._96_200_sampled.fits.cat  -COLUMNS columns/y1a1_final_spec_valid.columns  -INTERP 8 -VERBOSE 0

then back in the ../bpzv1/test/
%>markus_bpz_out_fits.py WL_CLASS.METACAL.rescaled.slr.cosmos.v2._96_200_sampled.fits.bpz

And for comparison the newer version 
%>./bpzv1.py test/bpzConfig.yaml test/WL_CLASS.METACAL.rescaled.slr.cosmos.v2._96_200_sampled.fits

--- results ---
Compare the outputs. The priors are identical [if old BPZ is correctly re-normalised -- it currently is not!], and the likelihoods are identical. Some small differences creep
in due to different interpolation algorithms. This is a smaller than 2% effect.

output:
array:
delta_z = z_max_post - REDSHIFT
metrics:
('median(delta_z), mean(delta_z), std(delta_z), outFrac(delta/(1+z) > 0.15), len(delta_z)', 

old-bpz
0.0275658443 0.003453 0.5619  30.46609 4999
new-bpz
0.027212330, -0.0069863, 0.56566, 31.1062, 4999)

-- old template order
0.02721233069248824, -0.0069863751258707587, 0.56566552010489746, 31.10622124424885, 

--new template order
 0.011479710711922331, -0.0032816341776811184, 0.56635402992871497, 29.90598119623924

--Standard template list: 8 interps 
0.011479710711922331, -0.0032816341776811184, 0.56635402992871497, 29.90598, 4999

intpolate x8 all template types with each other
0.019759900170284506, -0.017462470344914566, 0.557489, 29.76595, 4999

intpolate x3 all template types with each other
0.018568740987925159, -0.021051188088463273, 0.56483379518312071, 30.0260, 4999)

intpolate x5 all template types with each other
0.019447923890760777, -0.018850748000445672, 0.55952095548652392, 30.00600

================================
6) Saving and recovering the best fitting templates, and synthetic fluxes
================================

6.1 -- first run bpzv1.py with a .yaml file with the path to the SED_DIR set (default used), and provide an output pickle file name
e.g. 
output_sed_lookup_file: BPZ_templates.p
SED_DIR: %s../../templates/SED/

-- to save trouble, bpz will stop if the output file already exists!

6.2) read in resulting .fits files (or hdf5) and pickle file: 
e.g.
%>cd photoz-wg/redshift_codes/photoz_codes/bpzv1/test/
In python 2.7 

import cPickle as pickle
from astropy.table import Table

sed_file = pickle.load(open('BPZ_templates.p', 'r'))
print sed_file.keys()
data = Table.read('WL_CLASS.METACAL.rescaled.slr.cosmos.v2._96_200_sampled.BPZ.fits')
print data.keys()

#now get the template for the first galaxy:
template_id = data['TEMPLATE_ID'][0]

#access the lambda wavelength range, and get the SED of that (composite galaxy)
lamb, sed = sed_file['SED'][template_id]


#FYI the other keys in the output pickle file may also be useful.
#['template_type', 'flux_per_z_template_band', 'filters_dict', 'SED', 'z_bins', 'filter_order']
flux_per_z_template_band: the flux array of the synthetic "data" array(redshiftbins, template_type, band)
filter_order: describing the order of the template_type
filters_dict: describing what mags/fluxes/errors and offsets were applied to each filter
template_type: what fraction of one galaxy type or another the template consists of

================================
5) To do
================================
- deal with unseen / missing data
- cythonise code
- parrallelise code with MPI? -- currently with joblib


-- under construction ---
#first compile the cython code that allows a speed up of BPZ
%>python setup.py build_ext --inplace
