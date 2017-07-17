# BPZ v1

Description
-----------

This is the DESDM version of the bpzv1 code that will be used to estimate photo-z for DES during production as well as development/validation. Although the code can be run on any catalog size, the idea is that BPZ will be run per tile using as inputs coadd catalogs from the multi-epoch pipeline, or SQL queries to the database/archive, for a given TILENAME and TAGNAME (i.e. release)

In order to run BPZ on a given the following steps need to the followed:
1. Fetch the coadd input catalogs for a given TILENAME/TAGNAME from the DES archive.
2. Combine the input catalogs fetched above to create a 'color' catalog per FILENAME. This is the input file for BPZ that contain the MAGS/FLUXES.
3. Run BPZ

In the case when we get the catalogs using a SQL query to the DESDM Database, step 1 and 2 are combined.


Installation and Setup
----------------------

The recommended method of installation is via the DESDM EUPS software distribution. This method (although slow) will install all pre-requisites automatically. 

To install:

```bash
eups distrib install bpz 1.3.2+1 —nolocks
```

To setup (i.e. load) the code:
```bash
setup -v bpz 1.3.2+1
```

The quick guide to EUPS can be found [here](https://opensource.ncsa.illinois.edu/confluence/display/DESDM/The+Impatient%27s+Guide+to+DESDM+EUPS+installation)


Example 1:
----------

In this example we fetch the information needed to build the input 'color' catalog using SQL and a connection to the DESDM database for tilename DES2246-4457. We use the example yaml configuration file located in bpz/etc/bpz-comfig-example.yaml. Nearly all of the option in the YAML config file can be overridden using command-line arguments.


1. Create the input color catalog using a mix of MOF and Extractor catalogs

```bash
build_colorcat_sql -c bpz-comfig-example.yaml --tilename DES2246-4457  --outcat DES2246-4457_cats/DES2246-4457_r2583p01_color.fits 
```

This will download the coadd catalogs fits files and store them in the directory: DES2246-4457_cats

Note that to connect to the archive, you'll need a $HOME/.desservices.ini file. Here is an example, please update your credentials.

```
#
# DES services configuration
# Please modify the passwords accordingly
#
[db-desoper]
user = your-user-name
passwd = your-DESDM-Database-passwd
name = desoper
server = leovip148.ncsa.uiuc.edu
port = 1521

[db-dessci]
user = your-user-name
passwd = your-DESDM-Database-passwd
name = dessci
server = desdb.ncsa.illinois.edu
port = 1521
```

2. Run BPZ using the same config file and the 'color' catalog built in step 1.

```bash
bpzv1  -c bpz-comfig-example.yaml --incat DES2246-4457_cats/DES2246-4457_r2583p01_color.fits --n_jobs 6  --outbpz DES2246-4457_cats/DES2246-4457_r2583p01_bpz.fits
```




Example 2:
----------

In this example we will run BPZ over tilename DES2246-4457 and tag name Y3A1_COADD. We use the example yaml configuration file located in bpz/etc/bpz-comfig-example.yaml. Nearly all of the option in the YAML config file can be overridden using command-line arguments.

1. We fetch the files automatically from the DES Archive. The code needs to read the login credential from a 'DES Services' files. This file is need to live in $HOME/.desservices.ini  and it usually created by easyaccess.

```bash
fetch_catalogs  DES2246-4457  --tagname Y3A1_COADD --outpath DES2246-4457_cats --clobber --verbose 
```

2. Using the above inputs, construct a `color` fits catalog. This is the input to run BPZ. By default the fluxes/mags will be extinction corrected using E(B-V) user-supplied factor per band. To avoid extinction correction use the option --no-extinction.
```bash
build_colorcat -c bpz-comfig-example.yaml --incats DES2246-4457_cats/DES2246-4457_r2583p01_*cat.fits  --outcat DES2246-4457_cats/DES2246-4457_r2583p01_color.fits
```

3. Run BPZ using the same config file and the 'color' catalog built in step 2.

```bash
bpzv1  -c bpz-comfig-example.yaml --incat DES2246-4457_cats/DES2246-4457_r2583p01_color.fits --n_jobs 6  --outbpz DES2246-4457_cats/DES2246-4457_r2583p01_bpz.fits
```


