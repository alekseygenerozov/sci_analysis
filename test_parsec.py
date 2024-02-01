import sys
sys.path.append("/home/aleksey/code/python")
from sci_analysis import Star_parsec
import cgs_const as cgs

s1 = Star_parsec.Star(1.5 * cgs.M_sun, 1e9 * cgs.year,
                      ptrack="/home/aleksey/Dropbox/projects/Hagai_projects/mpe_collab/star_evolve_check/parsec.dat")
breakpoint()