import numpy as np
import sys
sys.path.append('../')
import bh_photo_z_validation as pval


ngals = 100
sigs = np.random.uniform(size=ngals) * 0.1 + 0.02
x = np.arange(401)/400.0 * 2
xcentr = x[1:] - (x[1] - x[0]) / 2.0
arr = np.zeros((ngals, len(xcentr)))

for i in range(ngals):
    h = np.histogram(np.random.normal(size=1e5) * sigs[i] + 2 * sigs[i], bins=x)
    arr[i, :] = h[0] * 1.0

cum = np.cumsum(arr, axis=1)
cum = (cum.T / cum[:, -1]).T

for point in [0.1, 0.5, 0.8, 0.9]:
    #determine which index is just <= each point value]
    print [np.amax(xcentr[c <= point]) for c in cum]
