from copy import copy

import numpy as np
import pandas as pd
from selection.algorithms.forward_step import forward_stepwise
from selection.distributions.discrete_family import discrete_family
from selection.constraints.affine import gibbs_test
from scipy.stats import norm as ndist

def maxT(FS, sigma=1, burnin=2000, ndraw=8000):
    """
    Sample from constraints, reporting largest
    inner product with T.

    Parameters
    ----------

    FS : selection.algorithms.forward_step.forward_stepwise

    """
    n, p = FS.X.shape

    projection = P = FS.P[-1]

    # recomputing this is a little wasteful

    if P is not None:
        RX = FS.X-P(FS.X)
    else:
        RX = FS.X
    scale = np.sqrt((FS.X**2).sum(0))

    if len(FS.variables) >= 1:
        con = copy(FS.constraints())
        linear_part = FS.X[:,FS.variables]
        observed = np.dot(linear_part.T, FS.Y)
        LSfunc = np.linalg.pinv(linear_part)
        conditional_con = con.conditional(linear_part.T,
                                          observed)

        _, Z, W, _ = gibbs_test(conditional_con,
                                FS.Y,
                                LSfunc[-1],
                                alternative='twosided',
                                sigma_known=True,
                                burnin=burnin,
                                ndraw=ndraw,
                                UMPU=False,
                                how_often=-1,
                                use_random_directions=False,
                                use_constraint_directions=False)

    else: # first step

        Z = np.random.standard_normal((ndraw, n)) * sigma
        W = np.ones(ndraw)

    null_statistics = np.fabs(np.dot(Z, RX) / scale[None,:]).max(1)
    dfam = discrete_family(null_statistics, W)
    observed = np.fabs(np.dot(RX.T, FS.Y) / scale).max()
    pvalue = dfam.cdf(0, observed)
    pvalue = max(2 * min(pvalue, 1 - pvalue), 0)
    return pvalue

def compute_pvalues(y, X, sigma=1., maxstep=np.inf,
                    compute_maxT_identify=True):
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

    FS = forward_stepwise(X, y, covariance=sigma**2 * np.identity(n))
    results = []

    for i in range(min([n, p, maxstep])):

        # maxT computed before next constraints
        # are added

        pval_maxT = maxT(FS, sigma)

        # take a step of FS

        FS.next()
        
        if compute_maxT_identify:
            pval_maxT_identify = FS.model_pivots(i+1, alternative='twosided',
                                          which_var=[FS.variables[-1]],
                                          saturated=False,
                                          burnin=2000,
                                          ndraw=8000)[0][1]
        else:
            pval_maxT_identify = np.random.sample()
        var_select, pval_saturated = FS.model_pivots(i+1, alternative='twosided',
                                         which_var=[FS.variables[-1]],
                                         saturated=True)[0]

        # now, nominal ones

        LSfunc = np.linalg.pinv(FS.X[:,FS.variables])
        Z = np.dot(LSfunc[-1], FS.Y) / (np.linalg.norm(LSfunc[-1]) * sigma)
        pval_nominal = 2 * ndist.sf(np.fabs(Z))
        results.append((var_select, pval_maxT_identify, pval_saturated, pval_nominal, pval_maxT))
            
    results = np.array(results).T
    return pd.DataFrame({'variable_selected': results[0].astype(np.int),
                         'maxT_identify_pvalue': results[1],
                         'saturated_pvalue': results[2],
                         'nominal_pvalue': results[3],
                         'maxT_pvalue': results[4]}), FS

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
