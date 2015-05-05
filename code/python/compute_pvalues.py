import numpy as np
import pandas as pd
from selection.algorithms.forward_step import forward_stepwise
from scipy.stats import norm as ndist

def compute_pvalues(y, X, sigma=1.):
    """
    Parameters
    ----------

    y : np.float(n)
        The target, in the model $y = X\beta$

    X : np.float((n, p))
        The data, in the model $y = X\beta$

    sigma : np.float
        Standard deviation of the gaussian distribution :
        The covariance matrix is
        `sigma**2 * np.identity(X.shape[0])`.
        Defauts to 1.

    Returns
    -------

    results : pd.DataFrame
        DataFrame with variables (variable_id, 
                            selective_pvalue,
                            saturated_pvalue,
                            nominal_pvalue)   

    """

    n, p = X.shape

    FS = forward_stepwise(X, y, covariance=sigma**2 * np.identity(n))
    results = []

    for i in range(min(n, p)):
        FS.next()
        var_select, pval_select = FS.model_pivots(i+1, alternative='twosided',
                                                  which_var=[FS.variables[-1]],
                                                  saturated=False,
                                                  burnin=2000,
                                                  ndraw=8000)[0]
        pval_saturated = FS.model_pivots(i+1, alternative='twosided',
                                         which_var=[FS.variables[-1]],
                                         saturated=True)[0][1]

        # now, nominal ones

        LSfunc = np.linalg.pinv(FS.X[:,FS.variables])
        Z = np.dot(LSfunc[-1], FS.Y) / (np.linalg.norm(LSfunc[-1]) * sigma)
        pval_nominal = 2 * ndist.sf(np.fabs(Z))
        results.append((var_select, pval_select, pval_saturated, pval_nominal))
            
    results = np.array(results).T
    return pd.DataFrame({'variable_selected': results[0].astype(np.int),
                         'selected_pvalue': results[1],
                         'saturated_pvalue': results[2],
                         'nominal_pvalue': results[3]})

def completion_index(selected, active_set):
    """
    Compute completion index from a sequence of 
    selected variables.

    The completion index is the first index
    of selected that contains the full active_set.

    Parameters
    ----------

    selected : []
        Sequence of selected variables.

    active_set : set
        Set of active variables.

    Returns
    -------

    idx : int
        Completion index.

    >>> selected = [1,3,2,4,6,7,8,23,11,5]
    >>> active = [1,4,8]
    >>> completion_index(selected, active)
    6
    """
    active_set = set(active_set)
    for i in range(len(selected)):
        if active_set.issubset(selected[:i]):
            return i-1
    return len(selected)-1

if __name__ == "__main__":
    from selection.algorithms.lasso import instance
    X, y, beta, active, sigma = instance(p=10, snr=4, rho=0.3)
    R = compute_pvalues(y, X, sigma=sigma)
    print completion_index(R['variable_selected'], active)
