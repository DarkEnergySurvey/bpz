
"""Below are functions encoding different priors we may have or want to have """


def des_sva1_prior():
    """   See Carles Sanchez bpz_session.tgz from will  bpz_session/prior_sva1_weights.py
    """
    a = {'E/S0': 2.370, 'Spiral': 1.977, 'Irr': 1.222}
    zo = {'E/S0': 0.358, 'Spiral': 0.262, 'Irr': 0.057}
    km = {'E/S0': 0.091, 'Spiral': 0.099, 'Irr': 0.101}
    k_t = {'E/S0': 0.300, 'Spiral': 0.097}
    fo_t = {'E/S0': 0.34, 'Spiral': 0.27}
    momin_hdf = 19
    minmag = 19

    return a, zo, km, k_t, fo_t, momin_hdf, minmag


def des_y1_prior():
    """     HDFN prior from Benitez 2000
    for Ellipticals, Spirals, and Irregular/Starbursts
    Returns an array pi[z[:],:nt]
    The input magnitude is F814W AB ~= i mag
    """
    # See Table 1 of Benitez00 https://arxiv.org/pdf/astro-ph/9811189v1.pdf
    #and eq 29 /30

    a = {'E/S0': 2.460, 'Spiral': 1.836, 'Irr': 1.180}
    zo = {'E/S0': 0.542, 'Spiral': 0.399, 'Irr': 0.134}
    km = {'E/S0': 0.112, 'Spiral': 0.101, 'Irr': 0.143}
    k_t = {'E/S0': 0.296, 'Spiral': 0.156}
    fo_t = {'E/S0': 0.291, 'Spiral': 0.550}
    momin_hdf = 20
    minmag = 18

    return a, zo, km, k_t, fo_t, momin_hdf, minmag

def des_y1pls_prior():
    """   Will's COSMOS trained prior (from Laigle etal. 2016 data)
    for Ellipticals, Spirals, and Irregular/Starbursts
    Returns an array pi[z[:],:nt]
    The input magnitude is i mag
    Adjusted mix of types to allow prior to go to m_i=17
    """
    # See Table 1 of Benitez00 https://arxiv.org/pdf/astro-ph/9811189v1.pdf
    #and eq 29 /30

    a = {'E/S0': 2.460, 'Spiral': 1.836, 'Irr': 1.180}
    zo = {'E/S0': 0.206, 'Spiral': 0.096, 'Irr': -0.295}
    km = {'E/S0': 0.112, 'Spiral': 0.101, 'Irr': 0.143}
    k_t = {'E/S0': 0.296, 'Spiral': 0.156}
    fo_t = {'E/S0': 0.707, 'Spiral': 0.878}
    momin_hdf = 17
    minmag = 17

    return a, zo, km, k_t, fo_t, momin_hdf, minmag


#Add any new priors here. Or we could BHM over them.
def new_priors():
    """    Example new prior. Must have same params as cosmos_Laigle_proirs and is called with calculate_priors
        the code must
        return a, zo, km, k_t, fo_t, momin_hdf

        With a structure shown below
    """
    # See Table 1 of Benitez00 https://arxiv.org/pdf/astro-ph/9811189v1.pdf
    #and eq 29 /30

    a = {'E/S0': 0, 'Spiral': 0, 'Irr': 0}
    zo = {'E/S0': 0., 'Spiral': 0., 'Irr': 0.}
    km = {'E/S0': 0., 'Spiral': 0., 'Irr': 0.}
    k_t = {'E/S0': 0., 'Spiral': 0.}
    fo_t = {'E/S0': 0., 'Spiral': 0.}
    momin_hdf = 20
    #what is the minimum magnitude that we will allow
    #clip all other values to this.
    minmag = 18
    return a, zo, km, k_t, fo_t, momin_hdf, minmag












"""
#Add any new priors here. Or we can BHM over them.
def example_proirs():
    "    Example new prior. Must have same params as cosmos_Laigle_proirs and is called with calculate_priors
        the code must 
        return a, zo, km, k_t, fo_t, momin_hdf

        With a structure shown below
    ""
    # See Table 1 of Benitez00 https://arxiv.org/pdf/astro-ph/9811189v1.pdf
    #and eq 29 /30 

    a = {'E/S0': 0, 'Spiral': 0, 'Irr': 0}
    zo = {'E/S0': 0., 'Spiral': 0., 'Irr': 0.}
    km = {'E/S0': 0., 'Spiral': 0., 'Irr': 0.}
    k_t = {'E/S0': 0., 'Spiral': 0.}
    fo_t = {'E/S0': 0., 'Spiral': 0.}
    momin_hdf = 20

    return a, zo, km, k_t, fo_t, momin_hdf
"""
