import numpy as np
import cgs_const as cgs


def psi1(Mbh):
    return 0.8 + 0.26 * (Mbh / cgs.M_sun / 1e6) ** 0.5

def psi2(ms):
    ms_sol = ms / cgs.M_sun
    return (1.47 + np.exp((ms_sol - 0.669) / 0.137)) / (1. + 2.34 * np.exp((ms_sol - 0.669) / 0.137))

def rt_ryu(Mbh, star):
    '''
    Fitting function for the tidal radius from Ryu et al 2020
    '''
    return psi1(Mbh) * psi2(star.ms) * (Mbh / star.ms) ** (1. / 3.) * star.rad

def rt_naive(Mbh, star):
    '''
    Simple estimate of the tidal radius
    '''
    return (Mbh / star.ms) ** (1. / 3.) * star.rad


