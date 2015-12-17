PhotoZ validation code

1. About

2. Calling the code

3. data formats

4. More information


1. ===== About this doc and the code =====
This doc describes how to call and use the photo-z validation code, photoz_metrics.py

Authors Ben Hoyle, Chris Bonnet.

Dependencies, 
numpy version 1.10.X
scipy version 0.16.1
pandas version 0.17.1
astropy version 1.0.7

If you cannot control your own python libraries, I recommend installing and using Canopy Enthought or Anaconda.


2.  === How to call the code ====
There are three ways of calling the code:

1)
%>./photoz_metrics.py

This is write an exampleValidation.yaml "test" file to the directory, describing how to define your own tests.
Then you can define the tests, and files, and run the code like:

%>./photoz_metrics.py exampleValidation.yaml


2) %>./photoz_metrics.py  PathToPointPredictions1.fits PathToPointPredictions2.fits ...

This will run all of the battery of tests, which are related to point predictions, and are found in the testConfig/ directory on each of the input files. It will also check that the fits files are formatted correctly, with the column names required by the tests.


3) %>./photoz_metrics.py  PathToPointPredictions1.fits PathToPointPredictions2.fits ...

For a detailed walkthrough of the code, run the ipython notebook
%>ipython notebook

and navigate to/ load the file
../notebooks/ValidationScriptExample.ipynb 


3. === Format of the files ===

Point predictions file must be a fits file, and it must have the columns:

[upper-case]
COADD_OBJECTS_ID  
Z_MC
MEDIAN_Z
MODE_Z
MEAN_Z
Z_SPEC -- spec-z of galaxy

[lower-case] (sorry for temporary mixing of types)
weights_valid

It must also have any extra columns that you use in the test file, in the correct case.


4. === More information ===
For an interactive demo open the ipython notebook

%>cd notebook/
%>ipython notebook

And then load/click on the ValidationScriptExample.ipynb

contact benhoyle1212@gmail.com, or see me in photo-z hipchat.

