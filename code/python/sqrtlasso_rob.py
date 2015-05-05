#!/software/free/statdept/bin/python

import sys
sys.path.append('/home/jtaylo/.local/lib/python2.7/site-packages')

from argparse import ArgumentParser

import numpy as np, pandas as pd
from selection.algorithms.sqrt_lasso import sqrt_lasso, choose_lambda

def main():
    parser = ArgumentParser(
        description= '''
        Run LASSO at fixed lambda and save p-values and intervals.

        Objective is

        \|Y-X\beta\|_2 + \lambda \|\beta\|_1
        
        ''')
    parser.add_argument('--X',
                        help='Design matrix (CSV).')
    parser.add_argument('--Y',
                        help='Response.')
    parser.add_argument('--lam',
                        help='Multiple of theoretical lambda for sqrt-LASSO.', default=0.9, type=float)
    parser.add_argument('--outfile',
                        help='Where to results.')

    args = parser.parse_args()
    X = np.loadtxt(args.X, delimiter=',')
    Y = np.loadtxt(args.Y, delimiter=',')
    
    L = sqrt_lasso(Y, X, args.lam * choose_lambda(X))
    L.fit(min_its=150, tol=1.e-10)
    
    intervals = []
    pvalues = []
    if L.active is not None and args.outfile:
        intervals = np.array(L.active_gaussian_intervals,
                             np.dtype([('variable', np.int),
                                       ('lower', np.float),
                                       ('upper', np.float)]))
        pvals = np.array(L.active_gaussian_pval,
                             np.dtype([('variable', np.int),
                                       ('pvalue', np.float)]))
        if np.any(intervals['variable'] != pvals['variable']):
            raise ValueError('variable lists are different!')

        df = pd.DataFrame({'variable':pvals['variable'],
                           'pvalue':pvals['pvalue'],
                           'lower':intervals['lower'],
                           'upper':intervals['upper']})
        df.reindex(columns=['variable', 'pvalue', 'lower', 'upper'])
        df.to_csv(args.outfile, index=False)
    else:
        if not args.outfile:
            sys.stderr.write('Active set was empty!\n')
        
if __name__ == "__main__":
    main()

                
