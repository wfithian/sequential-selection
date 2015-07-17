import os
from itertools import product
import pandas as pd, numpy as np

os.system('R CMD BATCH maxT.R')



for alpha in ['05', '10', '20']:

    sim_results = pd.read_csv('../snr_5_alpha_%s.csv' % alpha)

    getvar = lambda var: sim_results['_'.join([name, var])]

    names = []
    results = []
    #error_types = ['pscreen', 'FWER.mod', 'FWER.mod.cond', 'FDR.var', 'FDR.mod']

    name_map = dict(zip(['nominal', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                        ['Nominal', 'MaxT', 'MaxT identify', 'MaxT unknown', 'Saturated']))
    for test, rule in product(['nominal', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                              ['simple', 'forward', 'strong']):
        name = '_'.join([test, rule])
        screen = np.mean(getvar("screen"))
        FDR_var = np.mean(getvar("FDP_var"))
        FDR_model = np.mean(getvar("FDP_model"))
        FDR_model_cond = np.mean(getvar("FDP_model") * getvar("screen")) / screen
        FWER_model = np.mean(getvar("FWER_model"))
        FWER_model_cond = FWER_model / screen
        S_var = np.mean(getvar("S_var"))

        results.append((screen, FWER_model, FWER_model_cond, FDR_var, FDR_model, FDR_model_cond, S_var))
        names.append('%s %s' % (name_map[test], rule))

    # now knockoffs

    R = sim_results['knockoff_R']
    V = sim_results['knockoff_V']
    results.append((np.mean(sim_results['knockoff_screen']),
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.mean(V * 1. / np.maximum(R, 1)),
                    np.mean(R - V)))
    names.append('Knockoff')

    R = sim_results['knockoff_plus_R']
    V = sim_results['knockoff_plus_V']
    results.append((np.mean(sim_results['knockoff_plus_screen']),
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.mean(V * 1. / np.maximum(R, 1)),
                    np.mean(R - V)))
    names.append('Knockoff+')

    results = np.array(results)
    error_df = pd.DataFrame({'screen':results[:,0],
                             'fwer.mod':results[:,1],
                             'fwer.mod.cond':results[:,2],
                             'fdr.var':results[:,3],
                             'fdr.mod':results[:,4],
                             'fdr.mod.cond':results[:,5],
                             's.var':results[:,6]},
                            index=names)

    error_df = error_df.reindex_axis(['screen', 'fwer.mod', 'fwer.mod.cond', 'fdr.mod', 'fdr.mod.cond', 'fdr.var', 's.var'], axis=1)
    file('../../error_rates_%s.html' % alpha, 'w').write(error_df.to_html(float_format = lambda v: '%0.3f' % v))

    error_df.columns = pd.Index([r'$p_{\text{screen}}$',
                                 r'$\text{FWER}_{\text{mod}}$',
                                 r'$\text{FWER}_{\text{mod}} \vert \text{screen}$',
                                 r'$\text{FDR}_{\text{var}}$',
                                 r'$\text{FDR}_{\text{model}}',
                                 r'$\text{FDR}_{\text{mod}} \vert \text{screen}$',
                                 r'$\text{S}_{\text{var}}'])

    file('../../error_rates_%s.tex' % alpha, 'w').write(error_df.to_latex(float_format = lambda v: '%0.3f' % v).replace('\\_', '_'))

    print error_df.to_latex(float_format = lambda v: '%0.3f' % v)
