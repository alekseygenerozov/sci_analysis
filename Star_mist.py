# from sci_analysis import interpolate
from astropy.io import ascii
import numpy as np
from scipy.interpolate import interp1d

from bash_command import bash_command as bc
import cgs_const as cgs
import read_mist_models
from astropy.table import Table


def to_mist_file(s1):
    return "{0:05d}M.track".format(int(100 * s1))

class Star(object):
    def __init__(self, msi, age, metallicity = 0.02, eep_base = "eep_data/"):
        """
        Simple container for stellar data. Everything should be in cgs units...

        :param msi float: ZAMS stellar mass
        :param age float: Stellar age
        :param metallicity float: Stellar metallicity

        """
        self.msi = msi
        self.metallicity = metallicity
        self.eep_base = eep_base
        self.get_star_track()
        self.track = self.get_star_track()

        self.rad = 10.**np.interp(age, self.track['Tev(Myr)'], self.track['log10(R)']) * cgs.R_sun
        self.evolve_star(age)


    ##Make class methods??
    def get_star_track(self):
        """
        Get evolutionary track for star using sse as a backend for now

        :param my_star Star: Star object for which to get track

        :return: Return table with stellar data as a function of time
        :rtype: Astropy table
        """

        ##Find the closest eep files to use...
        ms_grid = np.genfromtxt(self.eep_base + "/masses") 
        ms_grid_cgs = ms_grid / 100 * cgs.M_sun
        ii = np.where(ms_grid_cgs <= self.msi)[0][-1]
        jj = np.where(ms_grid_cgs > self.msi)[0][0]
        ##Interpolate track using iso repository--using template...
        with open("input.example_template", "r") as ff:
            temp = ff.read()
            temp = temp.replace("TT1", ms_grid[ii] + "M.track")
            temp = temp.replace("TT2", ms_grid[jj] + "M.track")
            temp = temp.replace("EEP_BASE", self.eep_base)    
        with open("input.example", "w") as ff:
            ff.write(temp)

        with open("input.tracks_template", "r") as ff:
            temp = ff.read()
            temp = temp.replace("MM", str(self.msi / cgs.M_sun))
        with open("input.tracks", "w") as ff:
            ff.write(temp)
        bc.bash_command("/home/aleksey/software/iso/make_track ./input.tracks")
        ##Read in data and rename columns
        ##Set initial time to 0.
        track = read_mist_models.EEP("interpTrack")
        ages = track.eeps['star_age']
        ages = ages - ages[0]
        Mt  = track.eeps['star_mass']
        type = track.eeps['phase']
        log_R = track.eeps['log_R']
        track = Table(data=(ages / 1e6, Mt, log_R, type), names=['Tev(Myr)', 'Mt', 'log10(R)', 'type'])

        return track

    def evolve_star(self, t):
        """

        """
        track = self.track
        self.age = t
        ##Mist tracks don't include remnant phase. Hack to include such phases in a mc scheme...
        if t / (1e6 * cgs.year) > np.max(track['Tev(Myr)']):
            self.type = 20
            self.rad = 0
        else:
            self.ms = np.interp(t / (1e6 * cgs.year), track['Tev(Myr)'], track['Mt']) * cgs.M_sun
            self.rad = 10.**np.interp(t / (1e6 * cgs.year), track['Tev(Myr)'], track['log10(R)']) * cgs.R_sun
            self.type = float(interp1d(track['Tev(Myr)'], track['type'], kind='nearest')(t / (1e6 * cgs.year)))




