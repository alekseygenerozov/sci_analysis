import numpy as np
import cgs_const as cgs

def t_quad(m1, m2, m3, a1, a2, e2):
	P1 = 2. * np.pi * (a1**3. / (cgs.G * (m1 + m2)))**.5
	P2 = 2. * np.pi * (a2**3. / (cgs.G * (m1 + m2 + m3)))**.5

	return 16. / (30. * np.pi) * (m1 + m2 + m3) / m3 * P2**2. / P1 * (1. - e2**2.)**1.5


def t_oct(m1, m2, m3, a1, a2, e2):
	P1 = 2. * np.pi * (a1**3. / (cgs.G * (m1 + m2)))**.5
	P2 = 2. * np.pi * (a2**3. / (cgs.G * (m1 + m2 + m3)))**.5
	eps = a1 / a2 * e2 / (1 - e2**2.)

	return t_quad(m1, m2, m3, a1, a2, e2) / eps 