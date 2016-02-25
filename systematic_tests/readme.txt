This readme describes the tools/scripts found in this directory.

==============
Overview:
==============
1) check_retrieve_data.py  -- file to obtain data from DESDM

2) one_d_correlate.py   - this file correlates every column of a file against every other.

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
