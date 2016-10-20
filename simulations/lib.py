import numpy as np
from astropy.table import Table
import copy
import sys

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


def Deltaz_sigma_68_sigma_95(arr1, arr2, sample_weight=None):
    diff = DeltaZ(arr1, arr2)
    stats = getStats(diff, sample_weight=sample_weight)
    return -1 * HarmonicMean(stats['sigma_68'], stats['sigma_95'])
