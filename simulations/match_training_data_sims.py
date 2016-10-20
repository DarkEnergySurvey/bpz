import yaml
import glob
from astropy.table import Table, vstack
import numpy as np
import lib
from sklearn.neighbors import NearestNeighbors
import sys
import os
import copy

"""
Written by ben hoyle.

Knn matching of dataset_1 (e.g. real adata) and datasets_2 (e.g simulations)

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
"""


def help_message():
    print ('./match_training_data_sims.py [,-h, -help, PathToYamlFile.yaml')
    print ('e.g. PathToYamlFile.yaml = config/match_training_data_sims.yaml contains information about how the match is to be performed')
    print ('if it is not passed in, the routine will look for it in config/match_training_data_sims.yaml')
    print ('otherwise print this help message')
    sys.exit()


args = sys.argv[1:]
if '-h' in args:
    help_message()

if __name__ == "__main__":

    #get the location of this script.
    script_dir = '/'.join([i for i in os.path.realpath(__file__).split('/')[0:-1]])

    if len(args) == 0 and os.path.isfile(script_dir + '/config/match_training_data_sims.yaml') is False:
        help_message()

    path_to_config = script_dir + '/config/match_training_data_sims.yaml'
    if len(args) == 1:
        path_to_config = args[1]

    try:
        print ('reading {:}'.format(path_to_config))
        config = yaml.load(open(path_to_config, 'r'))
    except:
        print ('error reading {:}\n'.format(path_to_config))
        help_message()

    if lib.key_exist_not_None(config, 'nside_ring_match') is False:
        config['nside_ring_match'] = 128

    #get the sims columns and the data columns for the NN matching. Ensure the order is perserved of key/value pairs
    sim_cols = []
    data_cols = []
    for c1 in config['sims_to_data_feature_mapping']:
        sim_cols.append(c1)
        data_cols.append(config['sims_to_data_feature_mapping'][c1])

    #load training data and get RA/DEC
    inpD, outData = lib.dataSet(config['data_path']['DATA'], data_cols, [config['ID']['DATA'], config['ip_ring']['DATA']]).loadData()

    #config['sims_to_data_coord_mapping']['DATA'] == [RA, DEC]
    ip_ring_128_data = outData[config['ip_ring']['DATA']]

    uniq_ip_ring_128_data = np.unique(ip_ring_128_data)

    #load all sim data that sits within uniq_ip_ring_128_data
    sim_files = glob.glob(config['data_path']['SIMS'])

    inpS, outSims = lib.get_data_ip(uniq_ip_ring_128_data, sim_cols, sim_files, config['ID']['SIMS'], config['ip_ring']['SIMS'])

    #results table.
    res = Table()

    #perform analysis for each super-pixel
    for ip in uniq_ip_ring_128_data:

        #construct the matching routine
        knn = NearestNeighbors(n_neighbors= 1, algorithm='ball_tree')

           #extract data only for these IPs
        ins_ip = outSims[config['ip_ring']['SIMS']] == ip
        ind_ip = outData[config['ip_ring']['DATA']] == ip

        #if we have found some data, continue
        if np.sum(ins_ip) > 0 and np.sum(ind_ip) > 0:
            knn.fit(inpS[ins_ip])

            dist, indices = knn.kneighbors(inpD[ind_ip], return_distance=True)
            dist = np.ravel(dist)
            indices = np.ravel(indices)

            #construct a ip-only results table
            temp = Table()

            #store pertinent SIM information
            for i in outSims:
                temp[i] = outSims[i][indices]

            #store the columns we care about from the SIMS
            for i, col in enumerate(sim_cols):
                temp[col] = inpS[indices, i]

            #store the original DATA ID
            temp[config['ID']['DATA']] = outData[config['ID']['DATA']][ind_ip]

            #store the matchign distance for prosperity
            temp['MATCHING_DISTANCE'] = dist

            #if this is the first iteration
            if len(res) == 0:
                res = copy.deepcopy(temp)
            else:
                #stack results
                res = vstack([res, temp])
    if len(res) > 0:
        res.write(config['output_path'])