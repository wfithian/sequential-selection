#!/software/free/statdept/bin/python

# These lines are specific to miller so that the python used is my local one

# Usage example:
# ./fwdstep_rob.py --X=X.csv --Y=Y.csv --active=active.csv --outfile=results.csv --sigma=5. --outfile=results.csv

import sys, os
sys.path.append('/home/jtaylo/.local/lib/python2.7/site-packages')

# add the current directory to the path so we can find the 
# sequential-selection functions

sys.path.append(os.path.abspath(os.path.dirname("__file__")))

# The rest of the script should run as long as selection is installed
# and the selective-selection python directory is on the path

from argparse import ArgumentParser

import numpy as np, pandas as pd
from compare_methods import run

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
    parser.add_argument('--sigma',
                        help='Noise standard deviation. REQUIRED', type=float)
    parser.add_argument('--maxstep',
                        help='How many steps. Defaults to min(n,p)', type=int)
    parser.add_argument('--active',
                        help='Active set (CSV filename). Necessary ' + 
                        'to compute completion index. REQUIRED. ASSUMES 1-BASED INDEXING!')    
    parser.add_argument('--outfile',
                        help='Where to store output.')
    parser.add_argument('--burnin',
                        help='How many burnin steps for sampling.', default=2000, type=int)
    parser.add_argument('--ndraw',
                        help='How many sampling steps.', default=8000, type=int)

    try: 			
        args = parser.parse_args()
        X = np.loadtxt(args.X, delimiter=',')
        Y = np.loadtxt(args.Y, delimiter=',')
        # ASSUMING 1-BASED INDEXING!!!!
        active = np.loadtxt(args.active, delimiter=',').astype(np.int) - 1
    except:
        parser.print_help()
        return 

    if args.sigma is None:
        raise ValueError('sigma is needed to compute saturated p-values!')

    full_results = {}

    maxstep = args.maxstep or np.inf

    run(Y, X, args.sigma, active,
        full_results=full_results,
        maxstep=args.maxstep,
        burnin=args.burnin,
        ndraw=args.ndraw)

    # change output variables to 1-based indices

    for k in full_results.keys():
        if 'variable_selected' in k:
            full_results[k][0] = full_results[k][0] + 1

    if args.outfile:
        pd.DataFrame(full_results).to_csv(args.outfile, index=False)
        
if __name__ == "__main__":
    main()

                
