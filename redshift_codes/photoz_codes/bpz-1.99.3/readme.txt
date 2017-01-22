Directory Contents

--- readme.txt --
This file describing the other files

-- example.txt -- 
A worked example with a small example BPZ ready file

--- markus_bpz_out_fits.py ---
Executable file to convert current BPZ output to a redshift-validation ready point prediction fits files.

%>python markus_pbz_out_fits.py Paths*ToBPZPredictions.bpz


--- bpz_prior_wrapper.py ---
A handy wrapper to extract the prior used for BPZ. 
Call like:

%>python
from bpz_prior_wrapper import BPZ_PRIOR
c = BPZ_PRIOR()
#get the prior for a galaxy with magnitude 20.5
p, z = c.evaluate(20.5)
print ('at Imag=', 20.5)
print (z, p)


--- bpz_tcorr.py  --
BPZ main code. Edited / some improvements by Will.

Call like:
%>python bpz_tcorr.py  test_DES/exampleBPZfile.cat -COLUMNS columns/y1a1_final_spec_valid.columns  -INTERP 8
