import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

def log_interp(x, xs, ys, **kwargs):
	return np.exp(IUS(np.log(xs), np.log(ys), **kwargs)(np.log(x)))

def log_integral(x1, x2, xs, ys, **kwargs):
	'''
	Compute \int_{u1}^{u2} e^u y(u) du, where u is log(x), u1=log(x1), and u2=log(x2). 
	y(u) is conputed from an interpolating spline constructed from xs and ys
	'''
	us=np.log(xs)
	return IUS(us, ys*np.exp(us), **kwargs).integral(np.log(x1), np.log(x2))
