from __future__ import division

import os
import numpy as np, pandas as pd
from itertools import product
from selection.tests.instance import gaussian_instance, AR_instance

from rules import simple_stop, forward_stop, strong_stop
from compute_pvalues import compute_pvalues, completion_index

from selection.sampling.truncnorm import BoundViolation

from matplotlib import mlab

def summary(variables, pvalues, active, rule, alpha):
    """
    Compute R, V_var, V_model, screening
    for a sequence of pvalues.
    """

    completion_idx = completion_index(variables, active)
    R = rule(pvalues, alpha)
    V_var = R - len(set(active).intersection(variables[:R]))
    V_model = max(R - 1 - completion_idx, 0)
    screen = R > completion_idx
    FWER_model = (V_model > 0)
    FDP_model = V_model / max(R, 1)
    FDP_var = V_var / max(R, 1)
    S_var = R - V_var
    return R, V_var, V_model, screen, FWER_model, FDP_model, FDP_var, S_var

def simulate(n=100, p=40, rho=0.3, 
             signal=5,
             do_knockoff=False,
             full_results={},
             alpha=0.05,
             s=7,
             correlation='equicorrelated',
             maxstep=np.inf,
             compute_maxT_identify=True):

    if correlation == 'equicorrelated':
        X, y, _, active, sigma = gaussian_instance(n=n,
                                                   p=p,
                                                   rho=rho,
                                                   signal=signal,
                                                   s=s,
                                                   random_signs=False)
    elif correlation == 'AR':
        X, y, _, active, sigma = AR_instance(n=n,
                                             p=p,
                                             rho=rho,
                                             signal=signal,
                                             s=s)
    else:
        raise ValueError('correlation must be one of ["equicorrelated", "AR"]')

    full_results.setdefault('n', []).append(n)
    full_results.setdefault('p', []).append(p)
    full_results.setdefault('rho', []).append(rho)
    full_results.setdefault('s', []).append(len(active))
    full_results.setdefault('signal', []).append(signal)

    return run(y, X, sigma, active, 
               do_knockoff=do_knockoff,
               full_results=full_results,
               maxstep=maxstep,
               compute_maxT_identify=compute_maxT_identify,
               alpha=alpha)

def run(y, X, sigma, active,
        full_results={},
        do_knockoff=False,
        alpha=0.05,
        maxstep=np.inf,
        compute_maxT_identify=True,
        burnin=2000,
        ndraw=8000):

    n, p = X.shape
    results, FS = compute_pvalues(y, X, active, sigma, maxstep=maxstep,
                                  compute_maxT_identify=compute_maxT_identify,
                                  burnin=burnin,
                                  ndraw=ndraw)
    completion_idx = completion_index(results['variable_selected'], active)
    full_results.setdefault('completion_idx', []).append(completion_idx)

    for column in results.columns:
        for i in range(results.shape[0]):
            full_results.setdefault('%s_%d' % (str(column), i+1), []).append(results[column][i])

    for i in range(len(active)):
        full_results.setdefault('active_%d' % (i+1,), []).append(active[i])

    full_results.setdefault('alpha', []).append(alpha)

    if do_knockoff:

        # this will probably not work on miller
        import rpy2.robjects as rpy
        from rpy2.robjects import numpy2ri
        rpy.conversion.py2ri = numpy2ri.numpy2ri

        numpy2ri.activate()
        rpy.r.assign('X', X)
        rpy.r.assign('y', y)

        # knockoff

        rpy.r.assign('alpha', alpha)

        knockoff = np.array(rpy.r("""
        library(knockoff)
        knockoff.filter(X = X, y = y, fdr=alpha, threshold="knockoff")$selected
    """)) - 1
        numpy2ri.deactivate()
        knockoff_R = knockoff.shape[0]
        knockoff_V = knockoff_R - len(set(active).intersection(knockoff))
        knockoff_screen = knockoff_R > completion_idx

        knockoff_plus = np.array(rpy.r("""
        knockoff.filter(X = X, y = y, fdr=alpha, threshold="knockoff+")$selected
    """)) - 1
        numpy2ri.deactivate()
        knockoff_plus_R = knockoff_plus.shape[0]
        knockoff_plus_V = knockoff_plus_R - len(set(active).intersection(knockoff_plus))
        knockoff_plus_screen = knockoff_plus_R > completion_idx

        full_results.setdefault('knockoff_R', []).append(knockoff_R)
        full_results.setdefault('knockoff_V', []).append(knockoff_V)
        full_results.setdefault('knockoff_screen', []).append(knockoff_screen)

        full_results.setdefault('knockoff_plus_R', []).append(knockoff_plus_R)
        full_results.setdefault('knockoff_plus_V', []).append(knockoff_plus_V)
        full_results.setdefault('knockoff_plus_screen', []).append(knockoff_plus_screen)

    for pval, rule_ in product(['maxT_identify_pvalue',
                                'maxT_identify_unknown_pvalue',
                                'maxT_unknown_pvalue',
                                'saturated_pvalue',
                                'nominal_pvalue',
                                'nominalT_pvalue',
                                'maxT_pvalue'],
                               zip([simple_stop, 
                                    strong_stop,
                                    forward_stop],
                                   ['simple',
                                    'strong',
                                    'forward'])):
        rule, rule_name = rule_
        (R, 
         V_var, 
         V_model, 
         screen,
         FWER_model,
         FDP_model,
         FDP_var,
         S_var) = summary(np.asarray(results['variable_selected']),
                          results[pval],
                          active,
                          rule, 
                          alpha)

        pval_name = '_'.join(pval.split('_')[:-1])
        for (n, value) in zip(['R', 'V_var', 'V_model', 'FDP_model', 'FDP_var', 'S_var', 'FWER_model', 'screen'],
                              [R, V_var, V_model, FDP_model, FDP_var, S_var, FWER_model, screen]):
            full_results.setdefault('%s_%s_%s' % (pval_name, rule_name, n), []).append(value)

    return full_results, FS

def batch(outfile, nsim, **simulate_args):

    full_results = {}
    simulate_args['full_results'] = full_results
    for i in range(nsim):
        try:
            print 'run %d' % (i+1)
            simulate(**simulate_args)
            df = pd.DataFrame(full_results)
            df.to_csv(outfile, index=False)
        except KeyboardInterrupt: 
            raise KeyboardInterrupt('halting due to keyboard interruption')
        except BoundViolation:
            print("bound violation in sampling -- new instance drawn")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description= '''
Run a batch of simulations.
        ''')
    parser.add_argument('--nsim',
                        help='How many simulations to run.', type=int)
    parser.add_argument('--seed',
                        help='Random number seed (optional).', type=int)
    parser.add_argument('--alpha',
                        help='Level for FDR and FWER control.', type=float,
                        default=0.05)
    parser.add_argument('--signal',
                        help='Signal size.', type=float,
                        default=7)
    parser.add_argument('--nsample',
                        help='Sample size.', type=int,
                        default=100)
    parser.add_argument('--nfeature',
                        help='Number of features.', type=int,
                        default=40)
    parser.add_argument('--sparsity',
                        help='Sparsity level.', type=int,
                        default=7)
    parser.add_argument('--outfile',
                        help='Where to store results.')
    parser.add_argument('--maxstep',
                        help='How many steps should we take?',
                        type=int,
                        default=-1)
    parser.add_argument('--rho',
                        help='Correlation parameter (for AR or equi-correlated)',
                        type=float,
                        default=0.1)
    parser.add_argument('--correlation',
                        help='One of ["equicorrelated", "AR"]',
                        default='equicorrelated')


    args = parser.parse_args()
    if args.maxstep < 0:
        args.maxstep = np.inf

    if args.seed is not None:
        np.random.seed(args.seed)

    batch(args.outfile, 
          args.nsim, 
          maxstep=args.maxstep,
          alpha=args.alpha,
          n=args.nsample,
          p=args.nfeature,
          s=args.sparsity,
          rho=args.rho,
          correlation=args.correlation,
          do_knockoff=(args.nsample > 2 * args.nfeature))
