This readme describes the tools/scripts found in this directory.

==============
Overview:
==============
1) check_retrieve_data.py  -- file to obtain data from DESDM

2) one_d_correlate.py   - this file correlates every column of a file against every other.

3) make_hp_maps.py - this file extracts each column from a fits file, and uses Ra/DEc to make healpix maps

4) two_d_correlate.py - this file correlates two healpix files against each other (using PolSpice) and takes into account masks.

==============
Script Details
==============
1.1)
This is a generic script to match any tables in DESDM and return all of the column results as a fits files. If column names are replicated a random suffix is added to the returned rows.

This will connect to DESDM (you need a username, password and easy access all set up for this to work). You can join any tables that you can "see" if you were to log in.

--Usage--
To get help information about the script, and about a selected list of available tables.
%>./check_retrieve_data.py -h  

As an example;  to download data from DESDM, getting the WL sample of objects with BPZ redshift predictions and WL measured quantites, such as shear etc

%>./check_retrieve_data.py redshift-table=hoyleb.PHOT0Z_Y1_BPZ_V01 sample-table=hoyleb.WL_COADD_OBJECTS_ID_1M data-table=hoyleb.IM3SHAPE_Y1V1

If you only want 50 rows (in CO_ADD_ID order) add max-rows=50 to the command.

If the file already exists in directory you are calling from, then the file will *not* be re-downloaded.

SELECT A.*,B.*,C.*  FROM  hoyleb.WL_COADD_OBJECTS_ID_1M A JOIN hoyleb.PHOT0Z_Y1_BPZ_V01 B ON A.COADD_OBJECTS_ID=B.COADD_OBJECTS_ID  JOIN hoyleb.IM3SHAPE_Y1V1 C ON A.COADD_OBJECTS_ID=C.COADD_OBJECTS_ID where ROWNUM<5;   

==============
Script Details
==============

2.1)
This script accepts an input file (e.g. from part 1) and calculate the mutual information (MI, which is like a generalised Pearsons correlation coefficient) between each column and every other column.
For more details about MI see 
https://en.wikipedia.org/wiki/Mutual_information

THe script also calculates the KS-test statistic and the Pearson's correlation coefficient

Dependencies:
import bh_photo_z_validation as pval
import minepy as MINE

e.g.
pip install minepy
photoz-wg/validation/ to system path
import sys
sys.path.append('PathToValidationScripts')


For help see:
%> ./one_d_correlate.py -h  

Call the script like:
./one_d_correlate.py File1.fits [ file2.fits file..N.fits]
or
./one_d_correlate.py File*.fits 

output:
    CorrelationResults_1d_[FileName.fits].yaml

Next open the ipython notebook and load
photoz-wg/validation/Visualise_systematics_1d_correlation_output.ipynb

Run the code [changing file paths!] and examine the plots!

==============
Script Details
==============
3.1)
This script loads in RA/Dec + every other column "Col1" of a fits file in turn, and generates a healpix map of Col1. You can decided how to compress the data in each pixel. the default is to take the numpy.mean vaule of all Col1 data in each pixel, but you may want to count the number of data in each pixel, so use 'len' or identify the stdev of Col1 in each pixel, so use numpy.std , or your own statistic.

Example:

make_hp_maps.py mapFile[s].fits  [columns=Cols,To,Extract ra=RA dec=DEC z=Z_MEAN zbins=[0,0.3,0.50.7,0.9] nside=512 statistic=numpy.mean|len|numpy.std] 

You may also bin the data along Bin_Column, and make maps for all Col1 data in each bin_col bin 


==============
Script Details
==============
4.1) 
This script loads in a mask file (total_mask.fits) and calcalates the cross correlation between the two maps, for example those created by make_hp_maps.py 

Example
two_d_correlate.py pathTohealPix/FitFiles*.fits mask=maskFile, output_file=fileName [thetamax=thetamax] 
thetamax is the apodadisation scale, for details see ftp://ftp.iap.fr/pub/from_users/hivon/PolSpice/latest/README 
""-thetamax [dfloat|NO](NO)
   maximum value of angle theta used in the integrals to compute the
   power spectra from the correlation functions.""

