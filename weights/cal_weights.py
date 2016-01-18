"""
Python script to calualate relative weighting

"""                                           

### TODO, should we add error space in the weights space ? ###

import pandas as pd
import os
import fitsio
from fitsio import FITS,FITSHDR


def read_spectra():
    data = fitsio.read('Y1A1_GOLD101_Y1A1trainValid_14.12.2015.fits')
    df_data = pd.DataFrame(data.byteswap().newbyteorder())
    df_data['g-r'] = df.MAG_AUTO_G - df.MAG_AUTO_R 
    df_data['r-i'] = df.MAG_AUTO_R - df.MAG_AUTO_I 
    df_data['i-z'] = df.MAG_AUTO_I - df.MAG_AUTO_Z 
    df_data['z-y'] = df.MAG_AUTO_Z - df.MAG_AUTO_Y 
    return df_data


def read_WL_catalogue():
    store = pd.HDFStore('/Volumes/Data/bonnett/DES/Y1/gold/Y1_gold_WL_v1.h5')
    df = store['Y1_gold_WL_v1']
    df_data['g-r'] = df.mag_auto_g - df.mag_auto_r 
    df_data['r-i'] = df.mag_auto_r - df.mag_auto_i 
    df_data['i-z'] = df.mag_auto_i - df.mag_auto_z 
    df_data['z-y'] = df.mag_auto_z - df.mag_auto_y  
    
    return df                                                        


def read_LSS_catalogue():
    store = pd.HDFStore('/Volumes/Data/bonnett/DES/Y1/gold/Y1_gold_LSS_v0.0.1.h5')
    df = store['table']
    df_data['g-r'] = df.mag_auto_g - df.mag_auto_r 
    df_data['r-i'] = df.mag_auto_r - df.mag_auto_i 
    df_data['i-z'] = df.mag_auto_i - df.mag_auto_z 
    df_data['z-y'] = df.mag_auto_z - df.mag_auto_y  
    
    return df


def get_weights_c(df1,df2):

    df_w1 = df1[['Z_SPEC','COADD_OBJECTS_ID','COADD_OBJECTS_ID','g-r','r-i','i-z','MAG_AUTO_I','MAG_AUTO_R']]
    df_w2 = df2[['coadd_objects_id','coadd_objects_id','coadd_objects_id','g-r','r-i','i-z','mag_auto_i','mag_auto_r']]
    
    df_w1.to_csv('train_spec_weight',index=False,header=False,sep=' ')
    df_w2.to_csv('gold_phot_weight' ,index=False,header=False,sep=' ')
    
    os.system('./weights_code/weights.x train_spec_weight gold_phot_weight 5')
    os.system('mv nnweight.prw  weights_for_spec.prw &')
    
    weights = np.loadtxt('weights_for_spec.prw')
    weights = weights[:,1:3]
    weights[:,0] = weights[:,0].astype(int)
        
    weights = pd.DataFrame(weights,columns=['coadd_objects_id','weight'])
    weights.weight = weights.weight.astype('float32')
        
    assert weights.duplicated('coadd_objects_id').sum() == 0, 'Duplicates in coadd_objects_ids, suggest you run get_weights_c again' 
    
    return  weights
    
    
if __name__ == '__main__':
    
    df_spectra = read_spectra()
    df_WL_catalogue = read_WL_catalogue()
    df_LSS_catalogue = read_LSS_catalogue()
    
    LSS_weights = get_weights_c(df_spectra, df_LSS_catalogue)
    WL_weights = get_weights_c(df_spectra, df_WL_catalogue)
    
    
    