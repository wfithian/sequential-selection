#!/software/free/statdept/bin/python

# These lines are specific to miller so that the python used is my local one

# Usage example:
# ./fwdstep_rob.py --X=X.csv --Y=Y.csv --active=active.csv --outfile=results.csv --sigma=5. --outfile=results.csv

import sys, os

# add the current directory to the path so we can find the 
# sequential-selection functions

sys.path.append(os.path.abspath(os.path.dirname("__file__")))

# The rest of the script should run as long as selection is installed
# and the selective-selection python directory is on the path

from argparse import ArgumentParser

import numpy as np, pandas as pd
from selection.algorithms.sqrt_lasso import sqrt_lasso, choose_lambda

def main():
    parser = ArgumentParser(
        description= '''
        Run forward step a full min(n, p) steps, computing
        selected, saturated and nominal pvalues.

        ''')
    parser.add_argument('--X',
                        help='Design matrix (CSV filename). REQUIRED')
    parser.add_argument('--Y',
                        help='Response (CSV filename). REQUIRED')
    parser.add_argument('--lam',
                        help='What proportion of population lambda should we use?', type=float, default=0.5)
    parser.add_argument('--outfile',
                        help='Where to store output.')

    try: 			
        args = parser.parse_args()
        X = np.loadtxt(args.X, delimiter=',')
        Y = np.loadtxt(args.Y, delimiter=',')
    except:
        parser.print_help()
        return 

    results = {}

    L = args.lam * choose_lambda(X)
    SQL = sqrt_lasso(Y, X, args.lam)
    SQL.fit(min_its=100)

    V = SQL.active_gaussian_pval
    V['variable'] += 1
    if args.outfile:
        pd.DataFrame(V).to_csv(args.outfile, index=False)

        
if __name__ == "__main__":
    main()

                
