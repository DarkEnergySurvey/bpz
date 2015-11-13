import numpy as np


class metricBins:

    def __init__(self, group=None):
        self._group = None
        if group == 'photoz':
            self._group = group
            """Photo z bins go here
            """
            self._tomographic_bins = [0.0, 3.0]

        if group == 'weak_lensing':
            self._group = group
            """Weak lensing metrics go here

            """
            self._tomographic_bins = [0.0, 0.3, 0.6, 0.9, 1.3, 2.0]

        #otherwise _group is not set

    def tomographic_bins(self):
        return self._tomographic_bins

    def group(self):
        return self.group
    