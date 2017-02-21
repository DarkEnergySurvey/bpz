
"""
Here we try to pull together all the parts for this project:
1. for each survey, train from the sims and human flags
   - output human-flag incompleteness  
2. loop through an array of flag-limits, for each flag-limit, 
   calculate incompleteness and mean redshift bias of weighted 
   redshift distributions
3. do 2. N times to get error bars associated with random forest
4. for each survey the output is two plots of:
   - incompleteness as a function of flag-limit
   - mean redshift bias as a function of flag-limit
"""

from __future__ import print_function
import numpy as np
import calculate_weights
import astropy.io.fits as pf
import sys
from sklearn.ensemble import RandomForestRegressor

survey = sys.argv[1]
nid = int(sys.argv[2])

flaglim_array = np.arange(40)*1.0/10 + 0.1
file_names = ['vvds_deep1_sim_field4_feat.fits']
full_sim_name = 'vvds_deep1_sim._feature.fits'
data_dir = '/Users/chihwaychang/Desktop/Work/spectra_sims_new/spectra_sims/data_files/'

# first collect all info in all the human-shifted files
Flags = np.array([])
Features = np.zeros((0,29))
Z = np.array([])

for i in range(len(file_names)):
	specfile = pf.open(data_dir + file_names[i])[1].data
	Flags = np.concatenate((Flags, specfile['QOP']), axis=0)
	features = np.hstack((specfile['line_SN'],specfile['break_SN'].reshape((len(specfile['break_SN']),1))))
	features[features>1000.] = 0.
	features[np.isnan(features)] = 0.
	Features = np.concatenate((Features, features), axis=0)
	Z = np.concatenate((Z, specfile['Z']), axis=0)

# flag 6 is wierd, we'll call it 2.5 for now
Flags[Flags==6] = 2.5

# define half of it to be the training sample and half to be the validation sample 
ids_all = np.random.choice(np.arange(len(Z)), size=len(Z), replace=False)
N = len(Z)/2
ids_train = ids_all[:N]
ids_target = ids_all[N:]

train_small = Features[ids_train]
train_small[np.isnan(train_small)] = 0.
train_small[np.isinf(train_small)] = 0.
train_flag = Flags[ids_train]

target_small = Features[ids_target]
target_small[np.isnan(target_small)] = 0.
target_small[np.isinf(target_small)] = 0.
target_flag = Flags[ids_target]

# run random forest
rf = RandomForestRegressor(n_estimators=2000, max_leaf_nodes=20, max_features=12, min_samples_leaf=2, warm_start=True, random_state=nid+100)
rf.fit(train_small, train_flag)

# predict larger sim set
full_sim = pf.open(data_dir+full_sim_name)[1].data
full_sim_features = np.hstack((full_sim['line_SN'],full_sim['break_SN'].reshape((len(full_sim['break_SN']),1))))
full_sim_small = full_sim_features.copy()
full_sim_small[np.isnan(full_sim_small)] = 0.
full_sim_small[np.isinf(full_sim_small)] = 0.
predict_flag = rf.predict(full_sim_small)

gr = full_sim['OMAG_G']-full_sim['OMAG_R']
ri = full_sim['OMAG_R']-full_sim['OMAG_I']
iz = full_sim['OMAG_I']-full_sim['OMAG_Z']
targ_data = np.vstack((full_sim['OMAG_I'],gr,ri,iz))

zmean_true = []
zmean_incomplete_raw = []
zmean_incomplete_weighted = []
incompleteness = []

# loop through the flag limits
for i in range(len(flaglim_array)):
	flaglim = flaglim_array[i]
	incompleteness.append(len(predict_flag[predict_flag>flaglim])*1.0/len(predict_flag))
	zmean_true.append(np.mean(full_sim['Z']))
	zmean_incomplete_raw.append(np.mean(full_sim['Z'][predict_flag>flaglim]))
	id_incomplete = np.arange(len(full_sim))[predict_flag>flaglim]

	if len(id_incomplete)>10:
		print(len(id_incomplete))
		# with spectra selection
		spec_data = np.vstack((full_sim['OMAG_I'][id_incomplete],gr[id_incomplete],ri[id_incomplete],iz[id_incomplete]))
		train_weights = calculate_weights.best_lima_knn_weights(spec_data.T, targ_data.T)
		zmean_incomplete_weighted.append(np.sum(full_sim['Z'][predict_flag>flaglim]*train_weights)/np.sum(train_weights))
	else: 
		zmean_incomplete_weighted.append(-1.0)

np.savez('spectra_incompleteness_'+str(survey)+'_'+str(nid)+'.npz', flaglim=flaglim_array, 
	incompleteness=np.array(incompleteness), zmean_true=np.array(zmean_true),
	zmean_incomplete_raw=np.array(zmean_incomplete_raw), zmean_incomplete_weighted=np.array(zmean_incomplete_weighted))



