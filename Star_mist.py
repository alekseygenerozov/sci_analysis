# from sci_analysis import interpolate
from astropy.io import ascii
import numpy as np
from scipy.interpolate import interp1d

from bash_command import bash_command as bc
import cgs_const as cgs
import read_mist_models
from astropy.table import Table

import os

def init_to_final(mi):
    if mi < 8.0 * cgs.M_sun:
        return 0.6 * cgs.M_sun
    elif mi < 20.0 * cgs.M_sun:
        return 1.4 * cgs.M_sun
    else:
        return 10.0 * cgs.M_sun

def to_mist_file(s1):
    return "{0:05d}M.track".format(int(100 * s1))


class Star(object):
    def __init__(self, msi, age, metallicity = 0.02, eep_base = "eep_data/", iso_path = "iso/", phot_sys="", filters=()):
        """
        Simple container for stellar data. Everything should be in cgs units...

        :param msi float: ZAMS stellar mass
        :param age float: Stellar age
        :param metallicity float: Stellar metallicity (NOT CURRENTLY USED)
        :param iso str (iso/): Path to ISO code
        :param eep_base str (eep_data/): Path to eep_data
        :param phot_sys str (): Photometric system to use

        """
        self.msi = msi
        self.metallicity = metallicity
        self.eep_base = eep_base
        self.iso_path = iso_path
        self.phot_sys = phot_sys
        self.filters = np.atleast_1d(filters)
        self.phot = {}
        self.get_star_track()
        ##Set properties to the input age
        self.evolve_star(age)


    def get_star_track(self):
        """
        Get evolutionary track for star

        :param my_star Star: Star object for which to get track

        :return: Return table with stellar data as a function of time
        :rtype: Astropy table
        """

        ##Find the closest eep files to use...
        ms_grid_str = np.genfromtxt(self.eep_base + "/masses", dtype=str)
        ms_grid = np.genfromtxt(self.eep_base + "/masses")
        ms_grid_cgs = ms_grid / 100 * cgs.M_sun
        ii = np.where(ms_grid_cgs <= self.msi)[0][-1]
        jj = np.where(ms_grid_cgs > self.msi)[0][0]
        

        ##Interpolate track using iso repository--using template...
        with open(self.iso_path + "/input.example_template", "r") as ff:
            temp = ff.read()
            temp = temp.replace("TT1", ms_grid_str[ii] + "M.track")
            temp = temp.replace("TT2", ms_grid_str[jj] + "M.track")
            temp = temp.replace("EEP_BASE", self.eep_base)
        print("writing input files for iso", self.iso_path, os.getcwd())
        with open("./input.example", "w") as ff:
            print("test")
            ff.write(temp)

        with open(self.iso_path + "/input.tracks_template", "r") as ff:
            temp = ff.read()
            temp = temp.replace("MM", str(self.msi / cgs.M_sun))
        with open("input.tracks", "w") as ff:
            ff.write(temp)
        # bc.bash_command("cp {0}/input.nml .".format(self.iso_path))
        bc.bash_command("cp {0}/bc_table.list .".format(self.iso_path))
        bc.bash_command(self.iso_path +  f"/make_track input.tracks {self.phot_sys}")
        ##Read in data and rename columns
        ##Set initial time to 0.
        track = read_mist_models.EEP("interpTrack")
        track.eeps["star_age"] = (track.eeps["star_age"] - track.eeps["star_age"][0]) / 1e6
        self.track = Table(track.eeps)['star_age', 'star_mass', 'log_R', 'phase', 'log_Teff', 'log_L']

        # bc.bash_command("rm interpTrack")
        # bc.bash_command("rm input.tracks")
        # bc.bash_command("rm input.example")

        if self.phot_sys:
            cmd_track = read_mist_models.EEPCMD("interpTrack.cmd")
            cmd_track.eepcmds["star_age"] = (cmd_track.eepcmds["star_age"] - cmd_track.eepcmds["star_age"][0]) / 1e6
            cols = ["star_age"]
            cols = list(np.concatenate((cols, self.filters)))
            self.track_cmd = Table(cmd_track.eepcmds)[cols]
            # bc.bash_command("rm interpTrack.cmd")

    def evolve_star(self, t):
        """
        Evolve star to age t

        If star is dead it will have type=20, teff=-1, rad=0, and luminosity 0.
        """
        track = self.track
        self.age = t
        ##Mist tracks don't include remnant phase. Hack to include such phases in a mc scheme...
        if t / (1e6 * cgs.year) > np.max(track['star_age']):
            self.ms = init_to_final(self.msi)
            self.type = 20
            self.rad = 0.0
            self.teff = -1.0
            self.lum = 0.0
            if self.phot_sys:
                for filt in self.filters:
                    self.phot[filt] = np.inf
        else:
            self.ms = np.interp(t / (1e6 * cgs.year), track['star_age'], track['star_mass']) * cgs.M_sun
            self.rad = 10.**np.interp(t / (1e6 * cgs.year), track['star_age'], track['log_R']) * cgs.R_sun
            self.type = float(interp1d(track['star_age'], track['phase'], kind='nearest')(t / (1e6 * cgs.year)))
            self.teff = 10.**np.interp(t / (1e6 * cgs.year), track['star_age'], track['log_Teff'])
            self.lum = 10.**np.interp(t / (1e6 * cgs.year), track['star_age'], track['log_L']) * cgs.L_sun
            if self.phot_sys:
                track = self.track_cmd
                for filt in self.filters:
                    self.phot[filt] = np.interp(t / (1e6 * cgs.year), track["star_age"],  track[filt])



