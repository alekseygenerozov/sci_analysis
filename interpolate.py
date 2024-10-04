import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.interpolate import interp1d

def log_interp(x, xs, ys, **kwargs):
	return np.exp(IUS(np.log(xs), np.log(ys), **kwargs)(np.log(x)))

def log_integral(x1, x2, xs, ys, **kwargs):
	'''
	Compute \int_{u1}^{u2} e^u y(u) du, where u is log(x), u1=log(x1), and u2=log(x2). 
	y(u) is conputed from an interpolating spline constructed from xs and ys
	'''
	us=np.log(xs)
	return IUS(us, ys*np.exp(us), **kwargs).integral(np.log(x1), np.log(x2))

def get_crossing(absc, ords, thres):
	'''
	Find location where ords first crosses through thres -- Either from above or from below

	Useful for finding zero-crossings for numerical data
	'''
	x = ords - thres
	test=np.array([x[i]*x[i+1] for i in range(len(x)-1)])
	idx=np.where(test<=0)[0][0]

	return interp1d([x[idx], x[idx+1]], [absc[idx], absc[idx+1]])(0.)