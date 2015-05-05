import numpy as np, csv, hashlib, os.path
import uuid
from selection.algorithms.tests.test_forward_step import test_data_carving_IC
from data_carving import split_vals, vals, df, dname

# how many points do we want for each fraction

min_sample_size = 100000

vals = vals + [('split_frac', 0.9),
               ('ndraw', 6000),
               ('burnin', 5000)]

dtype = np.dtype([('n', np.int), 
                  ('p', np.int), 
                  ('s', np.int), 
                  ('rho', np.float), 
                  ('snr', np.float), 
                  ('split_frac', np.float), 
                  ('ndraw', np.int), 
                  ('burnin', np.int), 
                  ('method', 'S5'), 
                  ('null', np.bool), 
                  ('pval', np.float),
                  ('uuid', 'S40')])

for i in range(5000):
    for split_frac in split_vals[::-1]:
        opts = dict(vals)
        opts['split_frac'] = split_frac
        identifier = str(uuid.uuid1())
        fname = '%s/results_split_%0.2f.npy' % (dname, split_frac)
        opts['df'] = np.inf # degrees of freedom for noise

        test = not os.path.exists(fname)
        if not test:
            prev_results = np.load(fname)
            if split_frac not in [0.3, 0.4, 0.5]:
                test = prev_results.shape[0] < min_sample_size
            elif split_frac in [0.3, 0.4]:
                test = prev_results.shape[0] < 5000
            else:
                test = prev_results.shape[0] < 10000

        if test:
            (null_carve, 
             null_split, 
             alt_carve,
             alt_split,
             counter) = test_data_carving_IC(**opts)[:-2]

            params = [v for _, v in vals]
            results = []
            if os.path.exists('%s/screening.npy' % dname):
                prev_results = np.load('%s/screening.npy' % dname)
                screening = np.empty(prev_results.shape[0]+1, 
                                      prev_results.dtype)
                screening[:-1] = prev_results
                screening[-1] = (split_frac, counter, identifier)
                np.save('%s/screening.npy' % dname, screening)
            else:
                dtype_screen = np.dtype([('split', np.float),
                                         ('counter', np.float),
                                         ('uuid', 'S40')])
                screening = np.array([(split_frac, counter, identifier)],
                                     dtype_screen)
                np.save('%s/screening.npy' % dname, screening)

            results.extend([tuple(params) + ('carve', True, p, identifier) 
                            for p in null_carve])
            results.extend([tuple(params) + ('split', True, p, identifier)
                            for p in null_split])
            results.extend([tuple(params) + ('carve', False, p, identifier) 
                            for p in alt_carve])
            results.extend([tuple(params) + ('split', False, p, identifier) 
                            for p in alt_split])

            rec_results = np.array(results, dtype)
            if os.path.exists(fname):
                prev_results = np.load(fname)
                rec_results = np.hstack([prev_results, rec_results])
            np.save(fname, rec_results)
            print rec_results.shape, 1. / screening[screening['split'] == split_frac]['counter'].mean(), fname

