# from sci_analysis import interpolate
from astropy.io import ascii
import numpy as np
from scipy.interpolate import interp1d

from bash_command import bash_command as bc
import cgs_const as cgs
import read_mist_models
from astropy.table import Table
from astropy.io import ascii

def to_mist_file(s1):
    return "{0:05d}M.track".format(int(100 * s1))

def get_mist_prop_iso(isocmd, mass_age_data, field):
    ords = np.ones(len(mass_age_data)) * np.inf

    for ii, row in enumerate(mass_age_data):
        #If age is less than 1e5 year continue-since minimum age of isochrone is usually 1e5 yr...
        if row[1] / cgs.year < 1e5:
            continue
        ##Log Age of star
        log_my_age = np.log10(row[1] / cgs.year)
        age_ind = isocmd.age_index(log_my_age)
        my_isochrone = isocmd.isocmds[age_ind]

        ords[ii] = np.interp(np.log10(row[0] / cgs.M_sun), np.log10(my_isochrone['initial_mass']), my_isochrone[field], right=np.inf)

    return ords

def init_to_final(mi):
    if mi < 8.0 * cgs.M_sun:
        return 0.6 * cgs.M_sun
    elif mi < 20.0 * cgs.M_sun:
        return 1.4 * cgs.M_sun
    else:
        return 10.0 * cgs.M_sun

class Star(object):
    def __init__(self, msi, age, metallicity = 0.02, eep_base = "eep_data/", iso_path = "iso/", iso_file="MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.4_JWST.iso.cmd"):
        """
        Simple container for stellar data. Everything should be in cgs units...

        :param msi float: ZAMS stellar mass
        :param age float: Stellar age
        :param metallicity float: Stellar metallicity (NOT CURRENTLY USED)

        """
        self.msi = msi
        self.metallicity = metallicity
        self.eep_base = eep_base
        self.iso_path = iso_path
        self.track = self.get_star_track(iso_file)
        self.evolve_star(age)


    ##Make class methods??
    def get_star_track(self, my_iso_file):
        """
        Get evolutionary track for star using sse as a backend for now

        :param my_star Star: Star object for which to get track

        :return: Return table with stellar data as a function of time
        :rtype: Astropy table
        """
        isocmd = read_mist_models.ISOCMD(my_iso_file)


        return isocmd

    def evolve_star(self, t):
        """

        """
        track = self.track
        self.age = t
        ##Mist tracks don't include remnant phase. Hack to include such phases in a mc scheme...
        mt = get_mist_prop_iso(track, [[self.msi, t]], "star_mass")[0]
        lum = 10.**get_mist_prop_iso(track, [[self.msi, t]], "log_L")[0] * cgs.L_sun
        teff = 10.**get_mist_prop_iso(track, [[self.msi, t]], "log_Teff")[0]
        type = get_mist_prop_iso(track, [[self.msi, t]], "phase")[0]

        if not np.isinf(mt):
            self.ms = mt * cgs.M_sun
            self.rad = np.sqrt(lum / (4. * np.pi * cgs.sigma_sb * teff**4.))
            self.type = type
        ##Mass is not properly set here so we cannot really track evolution of remnants!!
        ##We can have a crude map here to deal with this...
        else:
            self.ms = init_to_final(self.msi)
            self.type = 20
            self.rad = 0

