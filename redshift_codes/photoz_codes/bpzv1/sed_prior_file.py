
"""Below are functions encoding different priors we may have or want to have """

def cosmos_Laigle():
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

    return a, zo, km, k_t, fo_t, momin_hdf

#Add any new priors here. Or we can BHM over them.
def new_proirs():
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

    return a, zo, km, k_t, fo_t, momin_hdf












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