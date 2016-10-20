import numpy as np
from astropy.table import Table
import copy
import sys
from joblib import Parallel, delayed

#helper functions
def get_function(function_string):
    import importlib
    module, function = function_string.rsplit('.', 1)
    module = importlib.import_module(module)
    function = getattr(module, function)
    return function

def key_exist_not_None(dct, ky):
    if ky not in dct:
        return False
    return dct[ky] is not None


#functions to load datasets
class dataSet:

    """The Class that loads a data set with specified features

    Parameters
    ----------
    dataFile:    str, or pyTable
            This string defines the path to the file
            Or is a pre-loaded PyTable
    featuresIn: str list/array
        This array list all the features. output an N-D array with the same order
    featuresOut: str list/array
        The output features as a dictionary
    Returns
    -------
    method() to load the data
    input data array
    output array

    """

    def __init__(self, dataFile, featuresIn, featuresOut):

        if isinstance(dataFile, str):
            self.mustloadData = True
            #guess the fileTime!
            self.fileType = None
            if (dataFile[-4:] == '.fit') or (dataFile[-7:] == '.fit.gz') or (dataFile[-8:] == '.fits.gz') or (dataFile[-5:] == '.fits'):
                self.fileType = 'fits'

            if (dataFile[-4:] == '.csv') or (dataFile[-7:] == '.csv.gz') or (dataFile[-8:] == '.txt.gz') or (dataFile[-5:] == '.txt'):
                self.fileType = 'ascii.csv'
        else:
            self.mustloadData = False

        self.loadData = self.loadFitsData
        self.dataFile = dataFile
        self.featuresOut = featuresOut
        self.featuresIn = featuresIn

    def load_data(self):
        if self.fileType is not None:
            return Table.read(self.dataFile, format=self.fileType)
        else:
            #guess the file type. This *may* work
            return Table.read(self.dataFile)

    def check_valid_file(self, d_, cols_):
        kys = d_.keys()
        for i in cols_:
            no_split = True
            for sp in ['-', '+', '*', '/']:
                if sp in i:
                    no_split = False
                    c1, c2 = i.split(sp)
                    if c1 not in kys:
                        print c1, ' not in file'
                        sys.exit()
                    if c2 not in kys:
                        print c2, ' not in file'
                        sys.exit()
            if no_split and (i not in kys):
                print i, ' not in file'
                sys.exit()

        return True

    def loadFitsData(self):

        if self.mustloadData:
            d = self.load_data()
        else:
            d = self.dataFile

        self.check_valid_file(d, self.featuresIn + self.featuresOut)

        ourRes = {}
        if (type(self.featuresOut) == list):
            for i, itm in enumerate(self.featuresOut):
                if '-' in itm:
                    mag1, mag2 = itm.split('-')
                    v = np.array(d[mag1].data - d[mag2].data)
                elif '+' in itm:
                    mag1, mag2 = itm.split('+')
                    v = np.array(d[mag1].data + d[mag2].data)
                elif '*' in itm:
                    mag1, mag2 = itm.split('*')
                    v = np.array(d[mag1].data * d[mag2].data)
                elif '/' in itm:
                    mag1, mag2 = itm.split('/')
                    v = np.array(d[mag1].data).astype(float) / np.array(d[mag2].data).astype(float)
                else:
                    v = np.array(d[itm].data)
                ourRes[itm] = v
            N = len(ourRes[itm])
        else:
            ourRes = d[self.featuresOut]
            N = len(ourRes)

        arr = np.zeros((N, len(self.featuresIn)))
        for i, itm in enumerate(self.featuresIn):
            if '-' in itm:
                mag1, mag2 = itm.split('-')
                v = np.array(d[mag1].data - d[mag2].data)
            elif '+' in itm:
                mag1, mag2 = itm.split('+')
                v = np.array(d[mag1].data + d[mag2].data)
            elif '*' in itm:
                mag1, mag2 = itm.split('*')
                v = np.array(d[mag1].data * d[mag2].data)
            elif '/' in itm:
                mag1, mag2 = itm.split('/')
                v = np.array(d[mag1].data).astype(float) / np.array(d[mag2].data).astype(float)
            else:
                v = np.array(d[itm].data)
            arr[:, i] = v

        return arr, ourRes

        """ To do: add support for .txt /.csv/.dat """

    def loadcsvfile(self):
        d = np.genfromtxt(self.dataFile)
        arr = 0
        ourRes = 0
        return arr, ourRes

def get_data_ip(ips, cols, files, ID_, IP_):
    """ get_data_ip(ips, cols, files, ID_, IP_)
        load "cols" data from within "files" and keep ID "ID_" and IP_RING "IP_"
        only keep data if IP_RING is in list of ips
    """

    #identify which data sit in all relevant pixels, by looping over data sets
    inpS = np.zeros(0)
    for i, fi in enumerate(files):
        #load cols to match from each file, and also the IP_RING and ID
        inpS_, outSims_ = dataSet(fi, cols, [ID_, IP_]).loadData()

        #determine which healpix IPs to keep, remove any NaNs
        keep = np.array([ip in ips for ip in outSims_[IP_]])
        for j in range(len(inpS_[0])):
            keep *= np.isfinite(inpS_[:, j])

        inpS_ = inpS_[keep]
        for j in outSims_:
            outSims_[j] = outSims_[j][keep]
        del keep

        #start keeping a set of all data in all allowed IPs
        if len(inpS_) > 0:
            if len(inpS) == 0:
                inpS = copy.copy(inpS_)
                outSims = copy.deepcopy(outSims_)
            else:
                inpS = np.append(inpS, inpS_, axis=0)
                for j in outSims:
                    outSims[j] = np.append(outSims[j], outSims_[j])

    #return data array for use in matching, and dictionary of ID's IP_RING values
    return inpS, outSims


#functions to calculate metrics
def getStats(diff_, sample_weight=None):
    if sample_weight is not None:
        prob = np.array(sample_weight, dtype=float)
        prob /= np.sum(prob, dtype=float)
        diff = np.random.choice(diff_, replace=True, size=1000000, p=prob)
    else:
        diff = copy.copy(diff_)

    diff = np.sort(diff)
    sig = np.std(diff)
    avg = np.mean(diff)
    median = np.median(diff)
    sig68 = (diff[int(len(diff) * 0.84075)] - diff[int(len(diff) * 0.15825)]) / 2.0
    sig95 = (diff[int(len(diff) * 0.977)] - diff[int(len(diff) * 0.023)]) / 2.0
    errsig68 = sig68 / np.sqrt(2 * len(diff))
    fracOut = np.sum(np.abs(diff) > 0.15) * 100.0 / len(diff)
    frac2s = np.sum(np.abs(diff) > 2 * sig68) * 1.0 / len(diff)
    frac3s = np.sum(np.abs(diff) > 3 * sig68) * 1.0 / len(diff)
    res = {"median": median, "outlierFrac": fracOut,
           "sigma_68": sig68, "errsig68": errsig68,
           "sigma_95": sig95, "sigma": sig, "mean": avg,
           "frac_2sig68": frac2s,
           "frac_3sig68": frac3s,
           }
    return res


def HarmonicMean(v1, v2):
    return v1 * v2 / (v1 + v2)


def DeltaZ(arr1, arr2):
    diff = (arr1 - arr2) / (1 + arr1)
    return diff


#wrapper for ML cost functions
def Deltaz_sigma_68_sigma_95(arr1, arr2, sample_weight=None):
    diff = DeltaZ(arr1, arr2)
    stats = getStats(diff, sample_weight=sample_weight)
    return -1 * HarmonicMean(stats['sigma_68'], stats['sigma_95'])


#ML helper functions
def getRandomParams(p=None):
    """ A master list of all possible parameters
    which have been modified in standard ML jobs
    Later steps remove those which are unneccesary
    """

    params = {}
    params['min_samples_leaf'] = np.random.randint(1, 100)
    params['max_depth'] = np.random.randint(1, 100)
    params['max_features'] = np.random.uniform()
    params['subsample'] = np.random.uniform()
    params['n_estimators'] = np.random.randint(10, 250)
    params['learning_rate'] = np.power(10, -2.5 * np.random.random())

    params['base_estimator__min_samples_leaf'] = np.random.randint(1, 100)
    params['base_estimator__max_depth'] = np.random.randint(1, 100)
    params['base_estimator__max_features'] = np.random.uniform()
    params['base_estimator__subsample'] = np.random.uniform()

    #this will use all available cores!!!
    params['n_jobs'] = -1
    params['loss'] = np.random.choice(['linear', 'square', 'exponential'])
    params['class_weight'] = np.random.choice(['auto', 'subsample', None])

    return params


def combineParameters(pCLF, pInBase):
    """combine two dictionaries, only keeping keys in dict pCLF"""
    pIn = copy.deepcopy(pInBase)
    for i in pIn:
        if i in pCLF.keys():
            pCLF[i] = pIn[i]
    return pCLF


def tree_loss(params_, clf_):
    """Fixes the sci-kit learn loss function errors, on a MLA by MLA basis
    By determining if the loss function is in the list of available functions
    If not set to one which is allowed"""
    from sklearn.ensemble import GradientBoostingRegressor as gr
    from sklearn.ensemble import GradientBoostingClassifier as gc
    from sklearn.ensemble import AdaBoostRegressor as ar
    from sklearn.ensemble import BaggingRegressor as br
    from sklearn.ensemble import BaggingClassifier as bc
    from sklearn.linear_model import RidgeClassifier as rgdc
    from sklearn.linear_model import Ridge as rgd
    from sklearn.linear_model import SGDRegressor as sgd
    from sklearn.linear_model import SGDClassifier as sgdc
    from sklearn.linear_model import LogisticRegression as lr
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as ldra
    from sklearn.tree import DecisionTreeClassifier as dtc
    #from sklearn.svm import SVR

    fix_items = {
        'loss': {
            gr: ['ls', 'lad', 'huber', 'quantile'],
            gc: ['deviance', 'exponential'],
            ar: ['linear', 'square', 'exponential'],
            sgd: ['squared_loss', 'huber', 'epsilon_insensitive', 'squared_epsilon_insensitive'],
            sgdc: ['hinge', 'log', 'modified_huber', 'squared_hinge', 'perceptron', 'squared_loss', 'huber', 'epsilon_insensitive', 'squared_epsilon_insensitive']
        },
        'max_samples': {
            br: [1.0],
            bc: [1.0]
        },
        'solver': {
            rgd: ['auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag'],
            rgdc: ['auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag'],
            lr: ['newton-cg', 'lbfgs', 'liblinear', 'sag'],
            ldra: ['svd', 'lsqr', 'eigen']
        },
        'class_weight': {
            ar: ['auto'],
            gc: ['auto'],
            gr: ['auto'],
            sgdc: ['auto', None],
            sgd: ['auto', None],
            rgdc: ['balanced', None],
            dtc: ['balanced', None]
        },
        'learning_rate': {
            sgdc: ['optimal'],
            sgd: ['invscaling']
        },
        'max_features': {
            bc: [1.0],
            br: [1.0]
        },
        'max_iter': {
            rgdc: [10000]
        }
    }

    for fix_param in fix_items.keys():
        for cl, allowed_param_vals in fix_items[fix_param].iteritems():
            if isinstance(clf_, cl):
                if fix_param in params_:
                    if params_[fix_param] not in allowed_param_vals:
                        params_[fix_param] = np.random.choice(allowed_param_vals)

    return params_


class Parrallelise:
    """ Creates a generic method to perform
    trivially Parrallel operations
    loop is a tuple of things that one wants to use in the declared function
    call like:
    from bh_parallelise import Parrallelise
    p = Parrallelise(n_jobs=2, method=myFunction, loop=[Tuple,Tuple,..]).run()
    Tuple is passed to myFunction and containes all the info required by the function

    Exampe call:

    #load all the files
    files = glob.glob('/Volumes/Untitled/DES/Y1_GOLD_V1_0_1/*.csv')

    arr = []
    for n, f in enumerate(files):
        arr.append(['ALL_DES_SPECTRA.fits', f])

    res1 = Parrallelise(n_jobs=5, method=stiltsMatch, loop=arr).run()

    """

    def __init__(self, n_jobs=-1, method=None, loop=None):
        self._n_jobs = n_jobs
        self._method = method
        self._loop = loop

    def run(self):
        results = Parallel(n_jobs=self._n_jobs)(delayed(self._method)(l) for l in self._loop)
        return results
