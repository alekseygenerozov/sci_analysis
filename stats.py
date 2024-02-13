from scipy.stats import ks_2samp, anderson_ksamp
import numpy as np
def get_all_comp_tests(s1, s2):
    """
    s1, s2 array-like -- Two samples to compare
    """
    comps = np.zeros(4)
    for ii, alt1 in enumerate(("two-sided", "less", "greater")):
        comps[ii] = (ks_2samp(s1, s2, alternative=alt1).pvalue)
    comps[3] = anderson_ksamp((s1, s2)).pvalue
    return comps