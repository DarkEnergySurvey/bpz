import numpy as np

"""
Author: Ben Hoyle
"""

def lima_knn_weights(dtrain, dtest, n_neighbors=None):
    """Determine weights for training data to make it resemble test data
    ala Lima et al, using K-NN
    #equ. 24 of http://www.fma.if.usp.br/~mlima/papers/LimCunOyaFriLinShe08.pdf
    dtrain = np.array( Ngals, Ndimensions)
    dtest = np.array( Ngals, Ndimensions)

    -- actually the code now more closely resembles that found in the des photoz-wg/weights/ directory.
    """
    #perfrom Knn in input features to be near output features.
    from sklearn.neighbors import NearestNeighbors

    #set radius if not set
    if n_neighbors is None:
        n_neighbors = 5

    neighd2 = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(dtest)
    #determine KNN ball around training, and then test data
    neighd1 = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(dtrain)

    #what is the density of training data around 1) train and 2) test data
    rad_tr = np.zeros(len(dtrain))
    nn_tr = np.zeros(len(dtrain))
    ind = np.arange(len(dtrain))
    split = np.array_split(ind, 50)

    for i in range(len(split)):
        dist_d2, _ = neighd1.kneighbors(dtrain[split[i]], n_neighbors=n_neighbors, return_distance=True)
        rad_tr[split[i]] = np.array([np.amax(i1) if len(i1) > 0 else 0 for i1 in dist_d2])
        del dist_d2, _

        nn_tr[split[i]] = [len(neighd2.radius_neighbors(dtrain[j].reshape(1, -1), radius=rad_tr[j], return_distance=False)[0]) for j in split[i]]

    #determine weights, ala Lima et al:
    #equ. 24 of http://www.fma.if.usp.br/~mlima/papers/LimCunOyaFriLinShe08.pdf
    #---now changed to match the code in photoz-wg/weights/weights_code/calcWeights.cpp
    ratio_train = (nn_tr / n_neighbors)

    #set any infs to 0
    ratio_train[np.isfinite(ratio_train) == False] = 0

    #if this is a terrible fit, then don't weight at all.
    if np.all(ratio_train == 0):
        ratio_train = np.ones(len(ratio_train), dtype=float)

    ratio_train = ratio_train / np.sum(ratio_train)

    return ratio_train


def best_lima_knn_weights(dtrain, dtest, n_neighbors_array=None, verbose=False, return_all_weights=False):
    """Iteratively determine the best radius to use in the K-NN for lima weights
    by selecting on the average K.L. of the distributions

    dtrain = np.array( Ngals, Ndimensions)
    dtest = np.array( Ngals, Ndimensions)
    """

    from weighted_kde import gaussian_kde as gss_kde
    from scipy.stats import entropy
    if n_neighbors_array is None:
        n_neighbors_array = np.arange(8) * 2 + 1

    #short cuts int array
    Nd_ = range(len(dtest[0]))

    #determine the axis of the *test* data. Will use this in KL test
    axis = [np.linspace(np.amin(np.append(dtrain[:, i], dtest[:, i])), np.amax(np.append(dtrain[:, i], dtest[:, i])), 100) for i in Nd_]

    #deteremine KDE of the test data along each axis
    test_distNd = [gss_kde(dtest[:, i]).evaluate(axis[i]) for i in Nd_]
    indN0 = [test_distNd[i] > 1e-99 for i in Nd_]
    #set up weights and best KL values
    training_weights = None
    best_KL = 99.0

    weights_ = {}
    #for each radius
    for n_neighbors in n_neighbors_array:

        #determine weights for this radius kNN
        weights = lima_knn_weights(dtrain, dtest, n_neighbors=n_neighbors)

        #normalise for KDE
        p = weights / np.sum(weights)

        #DKE new weighted training distribution
        train_distNd = [gss_kde(dtrain[:, i], weights=p).evaluate(axis[i]) for i in Nd_]
        #determine KL along each dimension
        #print weights, test_distNd[i], train_distNd[i], dtrain[:, i], dtest[:, i]
        KL = [entropy(train_distNd[i][indN0[i]], test_distNd[i][indN0[i]]) for i in Nd_]

        if verbose:
            print 'radius:', n_neighbors, '<KL>:', np.mean(KL), 'MI:', KL
            print 'weights minmax:', np.amin(p), np.amax(p)
        #if lower average KL than before, use these weights
        if np.isfinite(np.mean(KL)) and np.mean(KL) < best_KL:
            best_KL = np.mean(KL)
            training_weights = weights
        if return_all_weights:
            weights_[n_neighbors] = weights[:]
    
    if return_all_weights:
        return training_weights, weights_
    
    return training_weights


def test_lima_knn_weights():
    """Visual unit test of lima_knn_weights"""
    
    from weighted_kde import gaussian_kde as gss_kde
    import matplotlib.pyplot as plt
    Nd = 3
    Nsamp = 3e3
    a = np.random.normal(size=Nd * Nsamp * 1).reshape((Nsamp * 1, Nd)) * 4 + np.random.uniform() * 5
    b = np.random.normal(size=Nd * Nsamp * 10).reshape((Nsamp * 10, Nd)) * 2 + np.random.uniform() * 4

    X_plot = np.linspace(np.amin([np.amin(a), np.amin(b)]), np.amax([np.amax(a), np.amax(b)]), 1000)
    weights = best_lima_knn_weights(a, b)
    p = weights / np.sum(weights)

    for i in range(Nd):
        a1_ = gss_kde(a[:, i]).evaluate(X_plot)
        a_lim = gss_kde(a[:, i], weights=p).evaluate(X_plot)
        b_ = gss_kde(b[:, i]).evaluate(X_plot)

        f = plt.figure()
        plt.plot(X_plot, a1_, label='un weighted-sampled')
        plt.plot(X_plot, a_lim, label='lima weights')
        plt.plot(X_plot, b_, label='match to this')
        plt.legend()
    
    plt.show()