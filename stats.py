import numpy as np
from scipy.stats import ks_2samp, anderson_ksamp
from sklearn.neighbors import KernelDensity
# from scipy.stats import gaussian_kde


def get_all_comp_tests(s1, s2):
    """
    s1, s2 array-like -- Two samples to compare
    """
    comps = np.zeros(4)
    for ii, alt1 in enumerate(("two-sided", "less", "greater")):
        comps[ii] = (ks_2samp(s1, s2, alternative=alt1).pvalue)
    comps[3] = anderson_ksamp((s1, s2)).pvalue
    return comps

def get_likelihood_empirical(pdf_discrete, obs, transform=None, bandwidth=None):
    """

    :param pdf_discrete: PDF represented by discrete points
    :param obs: Observations
    :param transform: Variable transformation (defaults to None)

    :return: Likelihood
    """
    pdf_discrete_copy = np.copy(pdf_discrete)
    obs_copy = np.copy(obs)
    if transform is not None:
        pdf_discrete_copy = transform(pdf_discrete_copy)
        obs_copy = transform(obs_copy)
    if bandwidth is None:
        bandwidth = 'scott'

    pdf_estimate = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(pdf_discrete_copy[:, np.newaxis])
    return np.exp(np.sum(pdf_estimate.score_samples(obs_copy[:, np.newaxis])))