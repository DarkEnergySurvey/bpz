
"""Match two distriutions in N-d using abundance (percentile) matching"""
import numpy as np

def get_one_d_percentiles(arr, num_percentiles=1000):
    """order data by arr and get interger percentiles"""
    indo = np.argsort(arr)
    index = np.zeros(len(arr), dtype=int)
    for i, indx in enumerate(np.array_split(indo, num_percentiles)):
        index[indx] = i
    return index

def interger_order(nddata, num_percentiles=100):
    """get an n-data array of intergers describing where each data sits in each cloud"""
    int_order = np.zeros_like(nddata, dtype=int)
    for i in range(len(nddata[0])):
        int_order[:, i] = get_one_d_percentiles(nddata[:, i], num_percentiles=num_percentiles)
    return int_order


def match_data_clouds(data1, data2, num_percentiles=100):
    """get index of data2 that is abundance matched to data1"""
    int_ord_d1 = interger_order(data1, num_percentiles=num_percentiles)
    int_ord_d2 = interger_order(data2, num_percentiles=num_percentiles)

    #match id's in cloud 1 to an object in cloud 2

    from sklearn.neighbors import NearestNeighbors
    neighd2 = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(int_ord_d2)

    split = np.arange(len(data1), dtype=int)

    #50k mean split the file
    Ns = 50000
    if len(split) < Ns:
        split = [split]
    else:
        split = np.array_split(split, int(len(split) / Ns) + 1)

    match_d1_d2 = np.zeros(len(data1), dtype=int)

    for i in split:
        dist, id2 = neighd2.kneighbors(int_ord_d1[i], n_neighbors=1)
        id2 = np.ravel(id2)
        match_d1_d2[i] = id2
    return match_d1_d2


def test_nd_abundance_matching1():
    """Test in N-d this code uncorrelated gaussians"""

    import matplotlib.pyplot as plt
    import copy

    Nd = 3
    sig1, sig2 = np.random.uniform(size=Nd)*3, np.random.uniform(size=Nd)*3
    mu1, mu2 = np.random.uniform(size=Nd)*10, np.random.uniform(size=Nd)*10
    l1, l2 = 4e3, 4e5
    d1 = np.zeros((l1, Nd))
    d2 = np.zeros((l2, Nd))

    for i in range(Nd):
        d1[:, i] = np.random.normal(size=l1) * sig1[i] + mu1[i]
        d2[:, i] = np.random.normal(size=l2) * sig2[i] + mu2[i]


    ind1_2 = match_data_clouds(d1, d2, num_percentiles=100)

    d12 = copy.copy(d1)

    almost_black = '#262626'
    plt.rcParams['figure.figsize'] = (24, 16)
    plt.rcParams.update({'font.size': 12,
                         'axes.linewidth': 5,
                        'text.color': almost_black,
                        'xtick.major.size': 4,
                        'ytick.major.size': 4,
                        'legend.fancybox': True,
                        'figure.dpi': 300,
                        'legend.fontsize': 14,
                        'legend.framealpha': 0.8,
                        'legend.shadow': True,
                        'xtick.labelsize': 32,
                        'ytick.labelsize': 32})

    for i in range(Nd):
        d12[:, i] += (d2[ind1_2, i] - d1[:, i])

    f, axarr = plt.subplots(Nd, Nd)
    for i in range(Nd):
        for j in range(Nd):
            if j > i:
                axarr[i, j].plot(d1[:, i], d1[:, j], ',', alpha=0.6, rasterized=True, label='data1')
                axarr[i, j].plot(d2[:, i], d2[:, j], ',', alpha=0.6, rasterized=True, label='data2')
                axarr[i, j].plot(d12[:, i], d12[:, j], '*', alpha=0.6, rasterized=True, label='data1->2')
    plt.show()


def test_nd_abundance_matching2():
    """Test in N-d this code correlated gaussians"""
    Nd = 2

    ln =  [4000, 4000]

    for i in range(Nd):
        sig1 = [np.random.uniform()*3, np.random.uniform()*7]
        xx = np.array([np.random.uniform()*10, 30+np.random.uniform()*10])
        yy =  np.array([np.random.uniform()*30, 60+np.random.uniform()*30])
        means = [xx.mean(), yy.mean()]  
        stds = sig1
        corr = 0.2         # correlation
        covs = [[stds[0]**2          , stds[0]*stds[1]*corr], 
            [stds[0]*stds[1]*corr,           stds[1]**2]] 

        if i == 0:
            d1 = np.random.multivariate_normal(means, covs, ln[i])
        else:
            d2 = np.random.multivariate_normal(means, covs, ln[i])


    ind1_2 = match_data_clouds(d1, d2, num_percentiles=100)

    almost_black = '#262626'
    plt.rcParams['figure.figsize'] = (24, 16)
    plt.rcParams.update({'font.size': 12,
                         'axes.linewidth': 5,
                        'text.color': almost_black,
                        'xtick.major.size': 4,
                        'ytick.major.size': 4,
                        'legend.fancybox': True,
                        'figure.dpi': 300,
                        'legend.fontsize': 14,
                        'legend.framealpha': 0.8,
                        'legend.shadow': True,
                        'xtick.labelsize': 32,
                        'ytick.labelsize': 32})

    for i in range(Nd):
        d12[:, i] += (d2[ind1_2, i] - d1[:, i])

    f, axarr = plt.subplots(Nd, Nd)
    for i in range(Nd):
        for j in range(Nd):
            if j>=i:
                axarr[i, j].plot(d1[:, i], d1[:, j], ',', alpha=0.6, rasterized=True, label='data1')
                axarr[i, j].plot(d2[:, i], d2[:, j], ',', alpha=0.6, rasterized=True, label='data2')
                axarr[i, j].plot(d12[:, i], d12[:, j], '*', alpha=0.6, rasterized=True, label='data1-2')
    plt.legend()
    plt.show()
