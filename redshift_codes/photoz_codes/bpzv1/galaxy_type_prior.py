import sys
import numpy as np
import copy


"""
Author: Ben Hoyle
-- based heavily on bpz/ by Benitez 2000

This code is a new implementation of the galaxy type prior codes that are availabe from within BPZ

They can be used for any galaxy

Todo:
Include other priors as additional functions.
    call like:

    mag = 22.343
    GALPROIR = GALAXYTYPE_PRIOR(
                    z=np.arange(0.01, 3.5, 0.01),
                    tipo_prior='sed_prior_file.cosmos_Laigle', -- where does the prior sit e.g. sed_prior_file.py def cosmos_Laigle()
                    mag_bins=np.arange(18, 24.1, 0.1), --this is the standard bining. Results *do* depend on this choice
                    template_type_list=['E/S0', 'Spiral', 'Spiral', 'Irr','Irr','Irr'] #We need to know how many templates and of which galaxy type to determine adundances.
                    )

    print 'mag', mag
    g = GALPROIR.calculate_priors(mag)

    return as dictionary of priors for each galaxy type.
    print g['Irr'][0: 50]
    print g['E/S0'][0: 50]
    print g['Spiral'][0: 50]

    g = GALPROIR.evaluate(mag, frac_prior_from_types={'E/S0':1, 'Irr':0, 'Spiral':0})
    frac_prior_from_types tells us how to interpolate from 'E/S0' -> 'Spiral' -> 'Irr'. Also allowed (but perhaps is not physical) any combination of the above. eg. {'E/S0':0.5, 'Irr':0.5, 'Spiral':0}
        of {'E/S0':0.3, 'Irr':0.3, 'Spiral':0.3}
"""

def get_function(function_string):
    import importlib
    module, function = function_string.rsplit('.', 1)
    module = importlib.import_module(module)
    function = getattr(module, function)
    return function

class GALAXYTYPE_PRIOR:
    """
    GALAXYTYPE_PRIOR()
    generates a set of priors for a range of magnitudes, redshifts, and galaxy SED types. Allows fractional interpolation between the SEDs

    priors can be read from anywhere, e.g. sed_proir_file.cosmos_Laigle
    == sed_proir_filecosmos_Laigle()

    """

    def __init__(self, z=np.arange(0, 3.51, 0.01), 
                 tipo_prior='sed_proir_file.cosmos_Laigle', 
                 mag_bins=np.arange(18, 24.01, 0.01),
                 template_type_list=None):

        if template_type_list is None:
            print ("Please list the sed_types we will use: e.g.")
            print ("['E/S0', 'Spiral','Spiral', 'Irr','Irr','Irr']")
            print ("each type will have equal prob with other similar types")
            sys.exit()

        #what are the templates types
        self.template_type_list = template_type_list

        #how many in total (for later P normalisation)
        self.num_tmp_total = len(self.template_type_list)

        #how many per template type [for later P normalisation]
        self.num_tmp_type = {}
        for i in np.unique(template_type_list):
            self.num_tmp_type[i] = np.sum([i == tmp for tmp in template_type_list])

        #get the redshift range
        self.z = z

        #Extract our priors for a file, here.:
        self.a, self.zo, self.km, self.k_t, self.fo_t, self.momin_hdf = get_function(tipo_prior)()

        #mag bins to evaluate proir
        self.mag_bins = mag_bins

        #results will be stored here
        self.prior = None

        #construct the dictionary "self.prior"
        self.prepare_prior_dictionary()

    def prepare_prior_dictionary(self):
        """This function populates a dictionary of prior values. Called once at class creation"""
        self.prior = {}
        for m in self.mag_bins:
            self.prior[np.round(m, decimals=4)] = self.calculate_priors(m)

    def prepare_prior_dictionary_types(self, template_type_dict):
        """pre constructs all magnitude / template fractions priors
        for  a faster look-up. Magbins are stored as int(mag*100)"""
        dic = {}
        for m in self.mag_bins:
            ky = int(m*100)
            dic[ky] = {}
            for i in template_type_dict.keys():
                dic[ky][i] = self.evaluate(m, frac_prior_from_types=template_type_dict[i])
        return dic

    def evaluate(self, imag, frac_prior_from_types=None):
        "produced a prior for i-band magnitude imag, using ratios of frac_types: {'E/S0': 1, 'Spiral': 0, 'Irr': 0}"

        if frac_prior_from_types is None:
            print ("Which galaxy type / fraction do you want")
            print ("e.g for ellipticals {'E/S0':1}")
            print ("e.g for ellipticals {'E/S0':0.5, 'Spiral':0.5}")
            sys.exit()
        frac_types = copy.copy(frac_prior_from_types)

        for gal_type in frac_types.keys():
            if gal_type not in ['E/S0', 'Spiral', 'Irr']:
                print ("galaxy type is not known. Must be 'E/S0', 'Spiral', 'Irr'")
                sys.exit()

        #which Galaxy types do we have
        mag_keys = self.prior.keys()

        #which is the nearest i-mag that we have calucalted a prior
        near_key = np.argmin(np.abs(np.array(mag_keys) - imag))
       
        #store priors in a tempory array
        priors = self.prior[mag_keys[near_key]]

        #interpolate between types, based on frac_types
        norm = np.sum(frac_types.values())

        #the final prior is a linear combination of the different frac_prior_from_types
        final_prior = np.zeros_like(priors[frac_types.keys()[0]])
        for gal_typ in frac_types:
            frac_types[gal_typ] /= norm
            final_prior += priors[gal_typ] * frac_types[gal_typ]

        #renormalise prior
        final_prior /= np.sum(final_prior)

        return final_prior


    def calculate_priors(self, mag):
        """This an updated wrapper to the prior function, as found in the BPZ cosmos_Laigle_prior.py file

         HDFN prior from Benitez 2000
        for Ellipticals, Spirals, and Irregular/Starbursts
        Returns an array pi[z[:],:nt]
        The input magnitude is F814W AB ~= i mag
        """ 

        momin_hdf = copy.copy(self.momin_hdf)
        m = np.clip(np.round(mag, decimals=1), 18, 24)

        # See Table 1 of Benitez00 https://arxiv.org/pdf/astro-ph/9811189v1.pdf
        #and eq 29 /30
        a = copy.copy(self.a)
        zo = copy.copy(self.zo)
        km = copy.copy(self.km)
        k_t = copy.copy(self.k_t)
        z = copy.copy(self.z)

        # Fractions expected at m = 20 of each type
        # 35% E/S0
        # 50% Spiral
        # 15% Irr
        fo_t = copy.copy(self.fo_t)

        #we divide the probs of each galaxy type, by the number of SEDs for that type
        for gal_typ in fo_t:
            fo_t[gal_typ] = fo_t[gal_typ] / self.num_tmp_type[gal_typ]

        #the evolving fractions for each galaxy type is given by
        f_t = {}
        for gal_typ in fo_t:
            f_t[gal_typ] = fo_t[gal_typ] * np.exp(-1.0 * k_t[gal_typ] * (m - momin_hdf))

        #the remaining fraction  are Irr, again / # of Irr SED types
        f_t['Irr'] = np.round(1.0 - np.sum([f_t[gal_typ] * self.num_tmp_type[gal_typ] for gal_typ in f_t.keys()]), decimals=4) / self.num_tmp_type['Irr']

        #calculate probs P(T|mag)
        p_T_m0 = {}
        for gal_typ in ['E/S0', 'Spiral']:
            p_T_m0[gal_typ] = fo_t[gal_typ] * np.exp(-1.0*k_t[gal_typ] * (m - momin_hdf))

        #remaining galaxy-type prob goes to Irr galaxies.
        p_T_m0['Irr'] = 1.0 - np.sum([fo_t[gt] * np.exp(-1.0*k_t[gt] * (m - momin_hdf)) for gt in k_t.keys()])

        p_z_tmo = {}

        for gal_typ in ['E/S0', 'Spiral', 'Irr']:
            z_mt = zo[gal_typ] + km[gal_typ] * (m - momin_hdf)
            z_mt = np.clip(z_mt, 0.01, 15)
            expon = np.clip(np.power(z / z_mt, a[gal_typ]), 0., 700.)
            p_z_tmo[gal_typ] = np.power(z, a[gal_typ]) * np.exp(-1.0 * expon)
            
        #print 'p_T_m0[E/S0]', p_z_tmo['E/S0']
        #mutliply p_T_m0 by the prior on galaxy type f_t
        #and normalise

        for gal_typ in ['E/S0', 'Spiral', 'Irr']:
            p_z_tmo[gal_typ] /= np.sum(p_z_tmo[gal_typ])
            p_z_tmo[gal_typ][p_z_tmo[gal_typ] < 1e-2 / len(self.z)] = 0

            #comparision with BPZ exact! until here.
            if gal_typ == 'Spiral' and np.abs(m-22.343) < 0.05 and False:
                print 'm',m
                print 'f_t', f_t
                print p_z_tmo[gal_typ][0:50]

            p_z_tmo[gal_typ] /= np.sum(p_z_tmo[gal_typ])
            p_z_tmo[gal_typ] *= f_t[gal_typ]
            p_z_tmo[gal_typ] /= np.sum(p_z_tmo[gal_typ])

        #exact with BPZ to here [if we renormalise BPZ output]
        return p_z_tmo

if __name__ == '__main__':
    mag = 22.343
    GALPROIR = GALAXYTYPE_PRIOR(
                    z=np.arange(0.01, 3.5, 0.01),
                    tipo_prior='sed_prior_file.cosmos_Laigle',
                    mag_bins=np.arange(18, 24.1, 0.1),
                    template_type_list=['E/S0', 'Spiral', 'Spiral', 'Irr','Irr','Irr']
                    )

    #g = GALPROIR.evaluate(mag, frac_prior_from_types={'E/S0':1, 'Irr':0, 'Spiral':0})
    print 'mag', mag
    g = GALPROIR.calculate_priors(mag)
    print g['Irr'][0:50]

