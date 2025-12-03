import astropy.constants as const
import numpy as np

import astropy.units as u

class Nucleus(object):
    def __init__(self, Mbh, mstar, rho0, r0, gamma):
        self.Mbh = Mbh
        self.mstar = mstar
        self.rho0 = rho0
        self.r0 = r0
        self.gamma = gamma

    def n(self, r):
        return (self.rho0 / self.mstar) * (r / self.r0)**-self.gamma

    def sigma(self, r):
        return (const.G * self.Mbh / r / (1 + self.gamma))**.5

    def coll_rate(self, r, rtarget, mtarget):
        return (self.n(r) * np.pi * rtarget**2. * self.sigma(r) * (1. + (2. * const.G * (mtarget + self.mstar) / (self.sigma(r)**2. * rtarget)).to(""))).to(u.yr**-1.)

    def trx(self, r):
        return (0.34 * self.sigma(r)**3. / (const.G**2. * self.n(r) * self.mstar**2. * np.log(self.Mbh / self.mstar))).to(u.yr)