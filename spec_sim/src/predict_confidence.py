from sklearn.ensemble import RandomForestClassifier
#from numpy import genfromtxt, savetxt
import numpy as np
import astropy.io.fits as pyfits

def main():
    #dataset = genfromtxt(open('Data/train.csv','r'), delimiter=',', dtype='f8')[1:]    

    # read data
    classified = pyfits.open('sim_bcc_vvds_deep_features_classified.fits')[1].data
    test_data = pyfits.open('sim_bcc_vvds_deep_features.fits')[1].data
        
    target  = classified['QOP']
    train = np.hstack((classified['line_SN'],classified['break_SN'].reshape((len(classified['break_SN']),1))))

    test = np.hstack((test_data['line_SN'],test_data['break_SN'].reshape((len(test_data['break_SN']),1))))

    # clean bad values ***
    train[np.isnan(train)] = 0.
    test[np.isnan(test)] = 0.
    train[train>1000.] = 0.
    test[test>1000.] = 0.

    # construct features as a rank-ordered list of S/N in features?


    # add feature importance?


    # create and train the random forest
    rf = RandomForestClassifier(n_estimators=100, n_jobs=2)
    rf.fit(train, target)

    np.savetxt('results.csv', rf.predict(test), delimiter=',', fmt='%f')

if __name__=="__main__":
    main()

