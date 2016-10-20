#! /usr/bin/env python
import numpy as np
import sys
import os
import yaml
import sml_mla as ml
import lib
import cPickle as pickle
from joblib import Parallel, delayed
import copy
import textwrap
import inspect

"""
Author Ben Hoyle

Removes all the hassle of building a machine learning experiment. Please consult benhoyle1212@gmail.com before using for purposes outside of DES photo-z

%>standardised_ml_training.py config/exampleMLConfig.yaml

If called without a config file, an example will be written to disk

Possible ToDo:
extend from non DecisionTree based methods

"""

#loop over all trees to get a redshift dist
def allPredictions(cl, inp):
    pred = cl.predict(inp)

    pdf = np.zeros((len(pred), len(cl.estimators_)))
    if hasattr(cl, 'estimator_weights_'):
        weights = cl.estimator_weights_ / np.sum(cl.estimator_weights_)
        weightsInt = np.array(len(weights) * 3 / np.sum(weights) * weights, dtype=int)
        wpdf = np.zeros((len(pred), np.sum(weightsInt)))
        del pdf
    sm = 0
    for i, DT in enumerate(cl.estimators_):
        if isinstance(DT, (list, tuple, np.ndarray)):
            DT = DT[0]
        DTz = DT.predict(inp)

        if hasattr(cl, 'estimator_weights_'):
            #make a weighted pdf
            for j in range(weightsInt[i]):
                wpdf[:, sm] = DTz
                sm += 1
        else:
            pdf[:, i] = DTz

    if hasattr(cl, 'estimator_weights_'):
        return np.std(wpdf[:, 0:sm], axis=1)
    return np.std(pdf, axis=1)

# parralelisation function
# beware some global variables like score_metric is in here
def xvFitPredict(_clf, _score_metric, _indho):

    # train on non held out
    _intr = np.setdiff1d(np.arange(len(inputs)), _indho, assume_unique=True)

    inp_ = copy.deepcopy(inputs[_intr])
    outp_ = copy.deepcopy(outputs[_intr])

    #if we have training weights, then bootstrap resample from weights
    if train_weights:
        wt_ = np.array(train_weights_arr[_intr], dtype=float)
        wt_ = wt_ / np.sum(wt_)
        indWtrain = np.random.choice(np.arange(len(_intr), dtype=int), size=len(_intr), p=wt_, replace=True)
        inp_, outp_ = inp_[indWtrain], outp_[indWtrain]

    _clf.fit(inp_, outp_)

    dictXv = {}
    # calculate stats on held out fraction
    dictXv['predict'] = _clf.predict(inputs[_indho])

    if hasattr(_clf, 'estimators_'):
        dictXv['sigma'] = allPredictions(_clf, inputs[_indho])
    if hasattr(_clf, 'predict_proba'):
        dictXv['predict_proba'] = _clf.predict_proba(inputs[_indho])

    #we can process which objects get a prediction?

    # calcaulate score
    dictXv['scoreXv'] = _score_metric(outputs[_indho], dictXv['predict'], sample_weight=train_weights_arr[_indho])
    dictXv['outputs'] = outputs[_indho]
    dictXv['index'] = _indho
    dictXv['id'] = ID[_indho]
    return dictXv


def writeExampleConfig():
    import os.path
    import os
    path = os.getcwd() + '/'
    if os.path.isfile('exampleMLConfig.yaml') is False:

        f = open('exampleMLConfig.yaml', 'w')
        txt = textwrap.dedent("""
#add paths to data
trainDataPath: {0}""".format(path) + """data/IPRING.4096.CosmosXDES._409_500_.txt.Buzzard.sims.fits.clean.fits

#if blank use cross validation, else use this.
validDataPath:

#ID column. This is required
ID: ID

#add feature list
features: {
   input: [MAG_AUTO_G,MAG_AUTO_R,MAG_AUTO_I,MAG_AUTO_Z,MAG_AUTO_G-MAG_AUTO_R,MAG_AUTO_G-MAG_AUTO_I,MAG_AUTO_G-MAG_AUTO_Z,MAG_AUTO_R-MAG_AUTO_I,MAG_AUTO_R-MAG_AUTO_Z,MAG_AUTO_I-MAG_AUTO_Z],
   output: Z
}
machineLearningArchitecture: [sklearn.ensemble.RandomForestRegressor, sklearn.ensemble.AdaBoostRegressor, sklearn.ensemble.GradientBoostingRegressor]

base_estimator: sklearn.tree.DecisionTreeRegressor

#your outputs will be stored here. Classifier and relevant cross validation info
outPutFileBase: """ + """{0}data/outxv20.p


#if validDataPath empty then populate this cross validation
numberCrossValidation: 20

#should we have a small number of base learners only for the first N rounds?
small_iters_with_low_n_estimators: 4

#how many machines should we make, we save only the best. Empty=1
iterations: 150

#verbose output?
verbose: False

#should we add weights during validation?
valid_weights:

#should we add weights during training?
train_weights:

#number of cores (for adaboost in xval)
n_jobs: 6

#maximimum number of MLAs [useful if standarising how many machines are explored]
maximum_NMLAs: 150

#which score metric e.g. sklearn.metrics.f1_score | sklearn.metrics.r2_score | ml_metrics.log_loss | ml_metrics.Deltaz_sigma_68_sigma_95 <-
score_metric: lib.Deltaz_sigma_68_sigma_95

        """.format(path))
        f.write(txt)
        f.close()

if __name__ == "__main__":
    # get config file name
    args = sys.argv[1:]

    inArgs = {}
    for i in args:
        if '=' in i:
            k, v = i.split('=')
            inArgs[k] = v
        elif '.yaml' in i:
            inArgs['config'] = i

    print 'input arguments', inArgs

    if ('config' not in inArgs.keys()):
        print "please supply config=PathToConfigFile | or PathToConfigFile"
        print "a exampleConfig.yaml has written to this directory"
        writeExampleConfig()
        sys.exit()

    # load config file
    try:
        config = yaml.load(open(inArgs['config'], 'r'))
    except:
        print "error loading config file"
        sys.exit()

    print 'config', config

    verbose = False
    if 'verbose' in config:
        verbose = config['verbose']

    outPutClassifier = config['outPutFileBase']

    #enable very fast loading
    outPutClassifier_temp = outPutClassifier + '.temp.p'

    # load the features lists
    features = config['features']

    if isinstance(features['output'], (list, tuple, np.ndarray)):
        features['output'] = features['output'][0]

    if config['ID'] in features['output']:
        print "ID cannot be the output feature"

    outp_id = [features['output'], config['ID']]

    #Do we have weights during training?
    train_weights = False
    if 'train_weights' in config:
        if config['train_weights'] is not None:
            outp_id.append(config['train_weights'])
            train_weights = True

    if verbose:
        print 'loading training data'
    # load the data
    inputs, outputs = ml.dataSet(config['trainDataPath'], features['input'], outp_id).loadData()
    if verbose:
        print 'trainign data loaded', np.shape(inputs)

    #are there NaN's to remove
    ind = np.array([True] * len(inputs))
    for i in range(len(inputs[0])):
        ind *= np.isfinite(inputs[:, i])
    if np.sum(ind) != len(ind):
        print 'removing some NaNs in training file.'
        print 'this many:', np.sum(ind == False)
        inputs = inputs[ind]
        for i in outputs:
            outputs[i] = outputs[i][ind]

    ID = copy.deepcopy(outputs[config['ID']])
    if train_weights:
        train_weights_arr = outputs[config['train_weights']]
    else:
        train_weights_arr = np.ones(len(ID)) * 1.0

    train_weights_arr = np.array(train_weights_arr, dtype=float) / np.sum(train_weights_arr, dtype=float)

    outputs = outputs[features['output']]

    Ntr = len(outputs)

    # we do k-fold cross validation
    numXVal = config['numberCrossValidation']

    # number of iterations to explore
    iterations = 1
    if 'iterations' in config:
        if config['iterations'] is not None:
            iterations = config['iterations']

    # load Machine learning architecture
    MLA = [lib.get_function(i)() for i in config['machineLearningArchitecture']]

    if config['base_estimator'] is not None:
        for clf in MLA:
            p = clf.get_params()
            if 'base_estimator' in p:
                p['base_estimator'] = lib.get_function(config['base_estimator'])()
                clf.set_params(**p)

    score_metric = lib.get_function(config['score_metric'])

    res_best = {'score': -9999, 'N_MLAs': -1}

    if verbose:
        print ('here: outPutsystem_temp', outPutClassifier_temp)
    if os.path.exists(outPutClassifier_temp) is False:
        pickle.dump(res_best, open(outPutClassifier_temp, "w"))

    imax = 0
    if 'small_iters_with_low_n_estimators' in config:
        if config['small_iters_with_low_n_estimators'] is not None:
            imax = int(config['small_iters_with_low_n_estimators'])

    # keep a running record of how many MLAs have been explored.
    N_MLAs = 0

    for i in range(iterations):
        params = lib.getRandomParams()

        for j, clf in enumerate(MLA):

            # get parameters from this MLA
            curr_params = clf.get_params()

            predictions = {}

            if i < imax:
                params['n_estimators'] = np.random.randint(1, 15)

            # combine these with those specified
            curr_params = lib.combineParameters(curr_params, params)

            curr_params = lib.tree_loss(curr_params, clf)

            clf.set_params(**curr_params)

            if verbose:
                print clf
                print clf.get_params()

            results = {}

            #define held out (almost) equal size continuous arrays
            indHo = np.array_split(np.arange(Ntr), numXVal)

            # now parralise this step with joblib
            resXv = Parallel(n_jobs=config['n_jobs'])(
                    delayed(xvFitPredict)(clf, score_metric, indHold) for indHold in indHo)

            scoreXv = np.zeros(numXVal).astype(float)
            for k, dictXv in enumerate(resXv):
                scoreXv[k] = dictXv['scoreXv']
            score = np.mean(scoreXv)
            results['resXv'] = resXv

            # load temp file  this is the best current solution
            res_best = pickle.load(open(outPutClassifier_temp, "r"))
            if 'N_MLAs' in res_best:
                if res_best['N_MLAs'] > N_MLAs:
                    N_MLAs = res_best['N_MLAs']
            N_MLAs = N_MLAs + 1
            res_best['N_MLAs'] = N_MLAs

            #keep track of how many MLAs have been explored.
            pickle.dump(res_best, open(outPutClassifier_temp, "w"))

            print ('this score, best score;', score, res_best['score'])
            if (res_best['score'] < score):
                if verbose:
                    print ('new winner', params)
                    print ('best score', score)
                    print ('N_MLAs', N_MLAs)
                res_best['score'] = score
                res_best['features'] = features
                res_best['results'] = results
                res_best['config'] = config
                res_best['N_MLAs'] = N_MLAs

                clf.fit(inputs, outputs)

                # save the best clf
                res_best['clf'] = clf
                print (clf)
                del clf

                # save this system for future use
                #let's check again to see if we have learnt a better model in the mean time?
                res_best1 = pickle.load(open(outPutClassifier_temp, "r"))
                #do we still have a better fit?
                if res_best1['score'] < res_best['score']:
                    pickle.dump(res_best, open(outPutClassifier, "w"))
                    res_best['clf_params'] = res_best['clf'].get_params(deep=True)
                    del res_best['clf']
                    if 'results' in res_best:
                        del res_best['results']
                    pickle.dump(res_best, open(outPutClassifier_temp, "w"))

            #have we trained enough?
            if 'maximum_NMLAs' in config:
                if config['maximum_NMLAs'] is not None:
                    if res_best['N_MLAs'] > config['maximum_NMLAs']:
                        sys.exit()
