from copy import copy

import numpy as np
import pandas as pd
from selection.algorithms.forward_step import forward_step
from selection.distributions.discrete_family import discrete_family
from selection.constraints.affine import gibbs_test
from scipy.stats import norm as ndist

def compute_pvalues(y, X, sigma=1., maxstep=np.inf,
                    compute_maxT_identify=True,
                    burnin=2000,
                    ndraw=8000,
                    accept_reject_params=(100,15,2000)):
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

    compute_maxT_identify : bool
        If True, compute the maxT test having identified
        the variable and sign (i.e. conditioning
        on the variable added and its sign).

    Returns
    -------

    results : pd.DataFrame
        DataFrame with variables (variable_id, 
                            selective_pvalue,
                            saturated_pvalue,
                            nominal_pvalue)   

    """

    n, p = X.shape

    FS_identity = forward_step(X, y, covariance=sigma**2 * np.identity(n))
    FS_maxT = forward_step(X, y, covariance=sigma**2 * np.identity(n))
    FS_identity_U = forward_step(X, y, covariance=np.identity(n))
    FS_maxT_U = forward_step(X, y, covariance=np.identity(n))
    
    iter(FS_identity); iter(FS_maxT)
    iter(FS_identity_U); iter(FS_maxT_U)
    results = []

    for i in range(min([n, p, maxstep])):

        # take a step of FS

        pval_maxT = FS_maxT.next(compute_pval=True,
                                 use_identity=False,
                                 ndraw=ndraw,
                                 burnin=burnin)
        
        pval_maxT_U = FS_maxT_U.next(compute_pval=True,
                                     use_identity=False,
                                     ndraw=ndraw,
                                     burnin=burnin,
                                     sigma_known=False)
        
        if compute_maxT_identify:

            pval_maxT_identify = FS_identity.next(compute_pval=True,
                                                  use_identity=True,
                                                  ndraw=ndraw,
                                                  burnin=burnin)

            pval_maxT_identify_U = FS_identity_U.next(compute_pval=True,
                                                      use_identity=True,
                                                      ndraw=ndraw,
                                                      burnin=burnin,
                                                      sigma_known=False)

        else:
            FS_identity.next(compute_pval=False)
            pval_maxT_identify = np.random.sample()
        var_select, pval_saturated = FS_maxT.model_pivots(i+1, alternative='twosided',
                                                          which_var=[FS_maxT.variables[-1]],
                                                          saturated=True)[0]

        # now, nominal ones

        LSfunc = np.linalg.pinv(FS_maxT.X[:,FS_maxT.variables])
        Z = np.dot(LSfunc[-1], FS_maxT.Y) / (np.linalg.norm(LSfunc[-1]) * sigma)
        pval_nominal = 2 * ndist.sf(np.fabs(Z))
        results.append((var_select, pval_maxT_identify, pval_saturated, pval_nominal, pval_maxT, pval_maxT_U, pval_maxT_identify_U))
            
    results = np.array(results).T
    return pd.DataFrame({'variable_selected': results[0].astype(np.int),
                         'maxT_identify_pvalue': results[1],
                         'saturated_pvalue': results[2],
                         'nominal_pvalue': results[3],
                         'maxT_pvalue': results[4],
                         'maxT_unknown_pvalue': results[5],
                         'maxT_identify_unknown_pvalue': results[6]}), FS_maxT

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
    X, y, beta, active, sigma = instance(n=100, p=40, snr=5, rho=0.3)
    R, FS = compute_pvalues(y, X, sigma=sigma, maxstep=20)
    print completion_index(R['variable_selected'], active)
