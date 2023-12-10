# from sci_analysis import interpolate
from astropy.io import ascii
import numpy as np
from scipy.interpolate import interp1d

from bash_command import bash_command as bc
import cgs_const as cgs

class Star(object):
    def __init__(self, msi, age, metallicity = 0.02, template_file = "sse_template"):
        """
        Simple container for stellar data. Everything should be in cgs units...

        :param msi float: ZAMS stellar mass
        :param age float: Stellar age
        :param metallicity float: Stellar metallicity

        """
        self.template_file = template_file
        self.msi = msi
        self.metallicity = metallicity
        self.track = self.get_star_track()

        self.ms = msi
        track = self.track
        self.rad = 10.**np.interp(age, track['Tev(Myr)'], track['log10(R)']) * cgs.R_sun
        self.evolve_star(age)


    ##Make class methods??
    def get_star_track(self):
        """
        Get evolutionary track for star using sse as a backend for now

        :param my_star Star: Star object for which to get track
        :param template_file str (sse_template): Template file which we will use for evolution (Specific format is expected)

        :return: Return table with stellar data as a function of time
        :rtype: Astropy table
        """

        with open(self.template_file, "r") as ff:
            temp = ff.read()
            temp = temp.replace("MM", str(self.msi / cgs.M_sun))
            temp = temp.replace("ZZ", str(self.metallicity))

        with open("evolve.in", "w") as ff:
            ff.write(temp)

        bc.bash_command("./sse evolve.in")
        track = ascii.read("evolve.dat", comment="\*")
        track = track[track['Tev(Myr)'] >= 0.0]

        return track

    def evolve_star(self, t):
        """

        """
        track = self.track

        self.age = t
        self.ms = np.interp(t / (1e6 * cgs.year), track['Tev(Myr)'], track['Mt']) * cgs.M_sun
        self.rad = 10.**np.interp(t / (1e6 * cgs.year), track['Tev(Myr)'], track['log10(R)']) * cgs.R_sun
        self.type = float(interp1d(track['Tev(Myr)'], track['type'], kind='nearest')(t / (1e6 * cgs.year)))


