
import numpy as np
import pandas as pd

@profile
def main():
    pdf_file = 'test.h5'

    verbose = True

    z_bins = np.arange(0.1, 3.5, 0.1)
    n_gals = 10000000

    z_max_post = np.arange(n_gals) + 1.0
    mean = np.arange(n_gals) + 1.0
    sigma = np.arange(n_gals) + 1.0
    median = np.arange(n_gals) + 1.0
    mc = np.arange(n_gals) + 1.0
    sig68 = np.arange(n_gals) + 1.0
    KL_post_prior = np.arange(n_gals) + 1.0
    pdfs_ = np.random.uniform(size=(n_gals, len(z_bins))) + 1.0


    if verbose:
        print 'entering pdf'

    inds = np.array_split(np.arange(n_gals), 50)
    for ind in inds:

        cols = {'MEAN_Z': mean[ind], 'Z_SIGMA': sigma[ind], 'MEDIAN_Z': median[ind],
                    'Z_MC': mc[ind], 'Z_SIGMA68': sig68[ind], 'KL_POST_PRIOR': KL_post_prior[ind]}

        if verbose:
            print 'making data frame'

        df2 = pd.DataFrame(cols)

        if verbose:
            print 'saving data frame'

        df2.to_hdf(pdf_file, key='point_predictions', format='table', append=True, complevel=5, complib='blosc')
        del df2
        del cols
        
        if verbose:
            print 'saved data frame'


        post_dict = {'KL_POST_PRIOR': KL_post_prior[ind], 'MEAN_Z': mean[ind]}
        for ii in np.arange(len(z_bins)):
            post_dict['pdf_{:0.4}'.format(z_bins[ii] + (z_bins[0]+z_bins[1])/2.0)] = pdfs_[ind, ii]
        if verbose:
            print 'generating DataFrame'

        df2 = pd.DataFrame(post_dict)

        if verbose:
            print 'writing pdf'

        df2.to_hdf(pdf_file, key='pdf_predictions', format='table', append=True,complevel=5, complib='blosc')
        del df2
    del pdfs_
    #free memory
    
    del post_dict
    if verbose:
        print 'leaving pdf'

main()