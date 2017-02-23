# maybe change constraint to have mean wander over +/- 1 px? 
import cap_mpfit as mpfit
import numpy as np
import math as math

def line_plus_cont(x, p):
    return p[0]*np.exp(-(x-p[1])**2/(2*p[2]**2)) + p[3] + p[4] * x

def fit_line(p, fjac=None, x=None, y=None, err=None):
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    model = line_plus_cont(x, p)
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return([status, (y-model)/err])

def two_cont_break(x, p):
    # model is a linear continuum, plus a step at predefined x (x is actually two separate segments put together - the break occupies a number of pixels inbetween) and a second slope. Choose origin of second slope to be the break point, define intecept of second slope to be zero? Or should I just use two flat segments? Start with this and improve if needed.
    return p[math.ceil((x-p[2])/(p[2]*2))]
    
def fit_break(p, fjac=None, x=None, y=None, err=None):
    # 
    model = [two_cont_break(x[i], p) for i in range(len(x))]
    status = 0
    return([status, (y-model)/err])


def assemble_SN(wave, spec, noise, lines, breaks):
    # spectra should be linear, so get dlambda
    dlambda = wave[1]-wave[0]
    line_SN = []
    break_SN = []

    # 7 ang is roughly 1 pixel - which is what most will come out to be.
    parinf = [{'parname':'amplitude','fixed':0},{'parname':'mean','fixed':1},{'parname':'sigma','fixed':0},{'parname':'cont_level','fixed':0},{'parname':'slope','fixed':0}]
    for line in lines:
        wl = line[0]
        line_type = line[1]
        # amplitude, mean, sigma, continuum level, continuum slope
        # set some of these based on object magnitude?
        # chisq of first two iterations doesn't make any sense (in fact I'm not sure any of them do).
        mask = np.where((wave >= wl-50.) & (wave <= wl+50))[0]
        if len(mask) >= int(100/dlambda)-1:
            #print(spec[mask])
            p0 = [max(spec[mask]), wl, 7., 0., 0.]
            fa = {'x':wave[mask], 'y':spec[mask], 'err':noise[mask]}
            m = mpfit.mpfit(fit_line, p0, functkw=fa, parinfo=parinf, quiet=1)
            if line_type == 'E':
                m.params[0] = max([m.params[0], 0.])
            elif line_type == 'A':
                m.params[0] = min([m.params[0], 0.])
            SN_peak = abs(m.params[0])/m.perror[0]
            SN_flux = SN_peak * (0.67 / 0.6)
            line_SN.append(SN_flux)
        else:
            line_SN.append(0.)
 
    parinf = [{'parname':'cont_one'},{'parname':'cont_two'},{'parname':'break_point','fixed':1}]
    for feature in breaks:
        mask1 = np.where((wave >= feature[0]) & (wave <= feature[1]))[0]
        mask2 = np.where((wave >= feature[2]) & (wave <= feature[3]))[0]
        if (len(mask1) >= int((feature[1]-feature[0])/dlambda)-1) & (len(mask2) >= int((feature[3]-feature[2])/dlambda)-1):
            mask = np.hstack((mask1,mask2))
            p0 = [np.mean(spec[mask1]), np.mean(spec[mask2]), (wave[mask1][-1]+wave[mask2][0])/2.]
            fa = {'x':wave[mask], 'y':spec[mask], 'err':noise[mask]}
            m = mpfit.mpfit(fit_break, p0, functkw=fa, parinfo=parinf, quiet=1)
            break_str = m.params[1] - m.params[0]
            break_uncert = np.sqrt(m.perror[0]**2 + m.perror[1]**2)
            SN_break = break_str / break_uncert
            break_SN.append(max([SN_break,0.]))
        else:
            break_SN.append(0.)

    return line_SN, break_SN


"""

           #print 'status = ', m.status
           #if (m.status <= 0): print 'error message = ', m.errmsg
           #print 'parameters = ', m.params
           #print 'uncert = ', m.perror


###
# Lenz & Ayres 1992
# scaling law for Gaussians in range of noise models.
# S/N = C * (FWHM/d_lambda)^(1/2) * S/N_0
# for integrated line flux, C = 0.67 (some dependence on noise model)
# cannot get f_L error from error propegation of parameter errors. Boo.
# Probably a better way.
# also must take into account the S/N of cont??
# Could instead scale back from known error on the amplitude:
# C = 0.6 for S/N of peak from Lenz & Ayres. and fwhm is 0.7




# nice for simple things... Probably lots of power in it, but think I will use mpfit.
#from lmfit.models import GaussianModel

#mod = GaussianModel()

#pars = mod.guess(y, x=x)
#out  = mod.fit(y, pars, x=x)
#print(out.fit_report(min_correl=0.25))



#       MPFIT: file cap_mpfit.py included in the distribution
#import cap_mpfit as mpfit


#mp = mpfit.mpfit(self._fitfunc, parinfo=parinfo, quiet=1, ftol=1e-4)



#import mpfit


# 892:908
"""
