BPZ:
1. The invocation script is bpz_tcorr_fits.py if your catalogue is in fits format; bpz_tcorr.py for the traditional ascii input. Both output ascii results files (it'd be nice to output hdf5 directly for the probabilities). Each script needs a path to the ABcorr directory setting. Alternatively you could revert the line where that is mentioned (L98) and simlink ABcorr to AB.

2. For fits input, see WL_CLASS.rescaled.columns

3. One of the resampled validation files Ben made is included as an example to run on.
Invoke as,
>python bpz_tcorr_fits.py WL_CLASS.rescaled.fits -INTERP 8
>python post_process_bpz.py WL_CLASS.rescaled.bpz

*** YOU MAY NEED TO ALTER THE PYFITS IMPORT STATEMENT IN BOTH SCRIPTS ***
