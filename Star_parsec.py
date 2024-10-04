from astropy.io import ascii
import numpy as np
import cgs_const as cgs


def get_closest(val, absc):
    diff = np.abs((val - absc) / val)
    idx = np.argsort(diff)[0]

    return absc[idx]


def remnant_mass(msi):
    if msi < 8 * cgs.M_sun:
        return 0.6 * cgs.M_sun
    elif msi < 20. * cgs.M_sun:
        return 1.4 * cgs.M_sun
    else:
        return 10. * cgs.M_sun


def log_interp(xx, x, y, right=None):
    return 10.**np.interp(np.log10(xx), np.log10(x), np.log10(y), right=np.log10(right))


def semi_log_interp(xx, x, y, right=None):
    return np.interp(np.log10(xx), np.log10(x), y, right=right)


class Star(object):
    def __init__(self, msi, age, ptrack="parsec.dat"):
        """
        Simple container for stellar data. Everything should be in cgs units...

        :param float msi: ZAMS stellar mass
        :param float age: Stellar age

        """
        self.msi = msi
        self.age = age
        self.ms = None
        self.rad = None
        self.teff = None
        self.type = None
        self.lum = None
        self.K = None
        self.times = None
        self.track = None

        self.ptrack = ptrack
        self.get_star_track()
        self.evolve_star(age)

    def get_star_track(self):
        """
        Get evolutionary track for star using sse as a backend for now

        :return: Return table with stellar data as a function of time
        :rtype: Astropy table
        """
        track = ascii.read(self.ptrack)
        #Implicitly assumes time column should be sorted (OK?)
        self.times = np.unique(track['logAge'])
        self.track = track

    def evolve_star(self, t):
        """
        Evolve star to age t...
        """
        self.age = t
        snap1 = get_closest(np.log10(t / cgs.year), self.times)
        filt_age = np.isclose(self.track['logAge'], snap1)
        track = self.track[filt_age]
        mini_cgs = track["Mini"] * cgs.M_sun
        mass_cgs = track["Mass"] * cgs.M_sun

        # self.ms = np.interp(np.log10(snap1), track['Mini'], track['Mass'], right=remnant_mass(self.msi)) * cgs.M_sun
        self.ms = log_interp(self.msi, mini_cgs, mass_cgs, right=remnant_mass(self.msi))
        self.K = semi_log_interp(self.msi, mini_cgs, track['Kmag'], right=-1.0)
        self.type = semi_log_interp(self.msi, mini_cgs, track['label'], right=20.0)

        self.lum = 10.**semi_log_interp(self.msi, mini_cgs, track['logL'], right=-np.inf) * cgs.L_sun
        self.teff = 10.**semi_log_interp(self.msi, mini_cgs, track['logTe'], right=0.)
        if self.type > 10:
            self.teff = -1.0
        self.rad = (self.lum / (4.0 * np.pi * cgs.sigma_sb * self.teff**4.))**.5

