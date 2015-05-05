import numpy as np, os
import matplotlib.pyplot as plt
from matplotlib.mlab import rec2csv

split_vals = ([0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.875, 0.9,
              0.925, 0.95, 0.96, 0.97, 0.98, 0.99, 1.])
vals = [('n', 100),
        ('p', 50),
        ('s', 7),
        ('rho', 0.3),
        ('snr', 7.)
        ]
opts = dict(vals)

df = np.inf
if df == np.inf:
    dname = 'gaussian'
else:
    dname = 'df_%d' % df

def summary(save=True):

    results = []

    try:
        screen = np.load('%s/screening.npy' % dname)

        for split_frac in sorted(np.unique(split_vals)):
            fname = '%s/results_split_%0.2f.npy' % (dname, split_frac)
            data = np.load(fname)
            null_carve = np.array([d['pval'] for d in data if d['method'] == 'carve' 
                                   and d['null'] == True])
            null_split = np.array([d['pval'] for d in data if d['method'] == 'split' 
                                   and d['null'] == True])
            alt_carve = np.array([d['pval'] for d in data if d['method'] == 'carve' 
                                   and d['null'] == False])
            alt_split = np.array([d['pval'] for d in data if d['method'] == 'split' 
                                   and d['null'] == False])

            power_carve = np.mean(alt_carve < 0.05)
            power_split = np.mean(alt_split < 0.05)
            level_carve = np.mean(null_carve < 0.05)
            level_split = np.mean(null_split < 0.05)

            p_screen = 1. / np.mean(screen[screen['split'] == split_frac]['counter'])
            result = (split_frac, level_carve, level_split, power_carve, power_split, p_screen)
            results.append(result)

            print split_frac

        R = np.array(results, np.dtype([ \
                    ('split', np.float),
                    ('level_carve', np.float),
                    ('level_split', np.float),
                    ('power_carve', np.float),
                    ('power_split', np.float),
                    ('p_screen', np.float)]))

        if save:
            np.save('%s/summary.npy' % dname, R)
            rec2csv(R, '%s/summary.csv' % dname)
            os.system('cd %s; R CMD BATCH ../makeRplots.r' % dname)

        plt.clf()
        plt.plot(R['split'], R['p_screen'])
        plt.plot(R['split'], R['power_split'], label='split alt')
        plt.plot(R['split'], R['power_carve'], label='carve alt')
        plt.plot(R['split'], R['level_split'], label='split null')
        plt.plot(R['split'], R['level_carve'], label='carve null')
        plt.plot([0.5,1],[0.05,0.05], 'k--')
        plt.gca().set_xlim([0.5,1])
        plt.legend(loc='upper left')
        plt.savefig('%s/summary.pdf' % dname)

    except:
        print 'no results yet'
        pass
