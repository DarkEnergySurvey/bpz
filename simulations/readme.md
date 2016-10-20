


1) Knn matching of dataset_1 (e.g. real adata) and datasets_2 (e.g simulations) 

call like:
python match_training_data_sims.py config/match_training_data_sims.yaml

inside config/match_training_data_sims.yaml 
you can define which columns map from one file to another.

The output is a fits file of KNN (1) matches, with the matched columns dataset_2 and the ID from  dataset_1

Notes:
Ensure ID in dataset_1 is a different field name than ID in dataset_2 [this could be fixed if required.]
Assumes that dataset_1 is one file
Assumes that datasets_2 is a list of files. Use the wildcard *.
Assumes that all datasets already have healpix IPs as columns in the file. [use add_ip_ring.py before if not]


2) Add IP_RING_XXX column to files
call like
python add_ip_ring.py Nside PathToFile[s]

Where Nside should probaby==128 for our purposes.

Notes:
Current outputs a new file, replicating the first, with the file ending .IPXXX.fits
To do
improve rather than output a new file, add a new (unique) column with the IPRING_XXX column to the old file


3)Machine learning training 

Removes all the hassle of building a machine learning experiment. Please consult benhoyle1212@gmail.com before using for purposes outside of DES photo-z

%>standardised_ml_training.py config/exampleMLConfig.yaml

If called without a config file, an example will be written to disk

Possible ToDo:
extend from non DecisionTree based methods


4)DES ML predictions
Using a system trained in 3, make predictions on new data
./des_ml_output.py data/trainedMachine.p PathToFile[*s] COLUMS,TO,KEEP,COMMA,SEP

Removes all the hassle of obtaining pdfs from a trained machine (using standardised_ml_training.py with tree methods to generate DES formatted data predictions)

%>des_ml_output.py MLA.p files.* COLUMS,TO,KEEP ID,REDSHIFT

ToDo:
currently hard coded the parrelisation, using
n_jobs = 3

"""
