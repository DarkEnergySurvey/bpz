add_modest.py
-- function to add the Starg-Galaxy MODEST continous value and discreet classification to COADD and MOF photometry

cut_data.py
-- a general routine to make a selection on a column, based on >=< values 

identify_matching_columns.py
-- A wrapper to topcat "stilts" to match two files on a given column (COADDED_OBJECTS_ID)

markus_bpz_out_fits.py
--- A function to read the BPZ output prediction pdfs, and extract point predictions of interese

apply_science_sample_cuts.py
-- a script that performs post-processing to extract WL, LSS, YNGold science samples *after* photo-z predictions have been estimated.

lima_knn_weights.py
-- a function to perform python only data weighting based on lima et al. It's general, allowing as may columns as input "features" as required

nd_abundance_matching.py
-- pretty code to perfrom abundance matching between two data files

cosmos_science_sample_mapping.py
-- under development. Wrapper to nd_abundance_matching.py

gethpixDataDES.py
-- add a HEALPIX pixel ID to each column in the data

uploadDataDES.py
-- this used to upload data to des. It's been depreciated by the internal DESDM easyaccess function load_table

pdf_recalibration.py 
-- use this function to rescaled redshift PDFs (in the agreed .h5 photo-z wg format) such that they have more statiscally meaning as a pdf. Details in the Y1 validation photo-z paper.

Call like:
>pdf_recalibration.py PathToPDF.h5 out=OutputRSfile.p true_z=XXX bin_col=YYY

if true_z and bin_col are given, PathToPDF.h5  must contain /point_predictions/ with true_z e.g. true_z=Z_SPEC|REDSHIFT
or
>pdf_recalibration.py InputRSfile.p PathToPDF*.h5
--will load InputRSfile.p and apply the pdf rescaling to PathToPDF*.hdf5 and produce a new .h5 files with PathToPDF*.ReScaled.h5