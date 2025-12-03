import astropy.constants as const
import numpy as np


class Species(object):
    def __init__(self, rho0, r0, mstar, gamma):
        self.rho0 = rho0
        self.r0 = r0
        self.mstar = mstar
        self.gamma = gamma

    def rho(self, r):
        return self.rho0 * (r / self.r0) ** self.gamma

    def n(self, r):
        return self.rho(r) / self.mstar


class Nucleus(object):
    def __init__(self, Mbh):
        self.species = []
        self.Mbh = Mbh
        self.gamma0 = 2.0

    def add_species(self, rho0, r0, mstar, gamma):
        self.species.append(Species(rho0, r0, mstar, gamma))

    def sigma(self, r):
        return (const.G * self.Mbh / (1.0 + self.gamma0) / r) ** 0.5

    def trx(self, r):
        # ntot = np.sum([ss.n(r) for ss in self.species])
        nm2 = sum([ss.n(r) * ss.mstar**2.0 for ss in self.species])
        return (
            0.34
            * self.sigma(r) ** 3.0
            / (const.G**2 * nm2 * np.log(self.Mbh.to("Msun").value))
        )
