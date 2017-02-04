This directory contains an updated version of BPZ.


Authors: Ben Hoyle

It is usable, but still under construction.

example usage. 
%>./bpzv1.py PathToConfig.yaml PathToListofFitsFiles.fits

where PathToConfig.yaml is a configuration file, that defines how to run the code. See bpzConfig.py for an example




#compare the output of original BPZ and this version

cd ../bpz-1.99.3/
%> python bpz_tcorrv1.py ../bpzv1/test/WL_CLASS.METACAL.rescaled.slr.cosmos.v2._96_200_sampled.fits.cat  -COLUMNS columns/y1a1_final_spec_valid.columns  -INTERP 8 -VERBOSE 0

with 
%>./bpzv1.py test/bpzConfig.yaml test/WL_CLASS.METACAL.rescaled.slr.cosmos.v2._96_200_sampled.fits



#first compile the cython code that allows a speed up of BPZ

%>python setup.py build_ext --inplace



