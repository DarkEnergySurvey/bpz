import numpy as np
import yaml

class metricBins:

    def __init__(self, group=None):
        self._group = None
        if group == 'photoz':
            self.group = group
            """Photo z bins go here
            """
            self.tomographic_bins = [0.0, 3.0]

        if group == 'weak_lensing':
            self._group = group
            """Weak lensing metrics go here

            """
            self.tomographic_bins = [0.0, 0.3, 0.6, 0.9, 1.3, 2.0]

        #otherwise _group is not set

    def get_param(self, ky):
        if ky in self:
            return self[ky]
        return None

    def set_param(self, param):
        for i in param:
            self[i] = param[i]

def photoz():
    return metricBins('photoz')

def weak_lensing():
    return metricBins('weak_lensing')

