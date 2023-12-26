import cgs_const as cgs
import numpy as np
from numba import njit

@njit('float64(float64)')
def f(e):
	return (1. + 73. / 24. * e ** 2. + 37. / 96. * e ** 4.) / (1 + e) ** (7. / 2.)

@njit('float64(float64)')
def g(e):
	return (1. + 7. / 8. * e ** 2.) / (1 + e) **2.

@njit('float64(float64)')
def h(e):
	return e * (1. + 121 / 304. * e**2.) / (1. + e) ** (5. / 2.)

def rsch(mtot):
	return 2. * cgs.G * mtot / cgs.c**2.

def change_energy_gw_orb(a0, e0, m1, m2):
	"""
	Change in orbital binding energy over a single orbit. Peters 1964
	"""
	rp = a0 * (1. - e0)
	mtot = m1 + m2
	return 8. * np.pi / (5. * 2. ** .5) * f(e0) * (m1 * m2)**2. / (mtot)**3. * cgs.c ** 2. * (rp / rsch(mtot)) ** (-3.5)

def change_a_gw_orb(a0, e0, m1, m2):
	"""
	Change orbital semimajor axis over a single orbit. Peters 1964
	"""
	return -change_energy_gw_orb(a0, e0, m1, m2) * 2. * a0**2. / (cgs.G * m1 * m2)

def change_ecc_gw_orb(a0, e0, m1, m2):
	"""
	Change in eccentricity over a single orbit. Peters 1964
	"""
	mtot = m1 + m2
	rp = a0 * (1. - e0)
	return -304. / 15. * np.pi / 2.**1.5 * h(e0) * m1 * m2 / mtot**2. * (rp / rsch(mtot))**(-2.5)

def change_j_red_gw_orb(a0, e0, m1, m2):
	"""
	Change in reduced angular momentum over a single orbit. Peters 1964
	"""
	return -change_ecc_gw_orb(a0, e0, m1, m2) * (e0 / np.sqrt(1. - e0**2.))

def Tc(a0, m1, m2):
	return 5. * cgs.c**5. * a0**4. / (256. * cgs.G**3. * m1 * m2 * (m1 + m2))

def gw_inspiral_time(a0, e0, m1, m2):
	"""
	Gravitational wave inspiral time from Ilya Mandel + 2021 fit to Peters results...
	"""
	return Tc(a0, m1, m2) * (1.0 + 0.27 * e0**10. + 0.33 * e0**20. + 0.2 * e0**1000.) * (1. - e0**2.)**3.5

