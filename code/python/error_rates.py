import os
from itertools import product
import pandas as pd, numpy as np

for alpha in ['05', '10', '20']:

    sim_results = pd.read_csv('../snr_5_alpha_%s.csv' % alpha)

    getvar = lambda var: sim_results['_'.join([name, var])]

    names = []
    results = []
    guarantees = []

    name_map = dict(zip(['nominal', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                        ['Nominal', 'MaxT', 'MaxT identify', 'MaxT unknown', 'Saturated']))
    for test, rule in product(['nominal', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                              ['simple', 'forward']):
        name = '_'.join([test, rule])
        screen = np.mean(getvar("screen"))
        FDR_var = np.mean(getvar("FDP_var"))
        FDR_model = np.mean(getvar("FDP_model"))
        FDR_model_cond = np.mean(getvar("FDP_model") * getvar("screen")) / screen
        FWER_model = np.mean(getvar("FWER_model"))
        FWER_model_cond = FWER_model / screen
        S_var = np.mean(getvar("S_var"))

        if test in ['maxT', 'maxT_identify', 'maxT_unknown']:
            if rule == 'simple':
                guarantees.append((False, True, True, True, False, False))
            else:
                guarantees.append((False, False, False, True, False, False))
        else:
            guarantees.append((False, False, False, False, False, False))
        results.append((screen, FWER_model, FWER_model_cond, FDR_model, FDR_var, S_var))
        names.append('%s %s' % (name_map[test], rule))

    # now knockoffs

    R = sim_results['knockoff_R']
    V = sim_results['knockoff_V']
    results.append((np.mean(sim_results['knockoff_screen']),
                    np.nan,
                    np.nan,
                    np.nan,
                    np.mean(V * 1. / np.maximum(R, 1)),
                    np.mean(R - V)))
    names.append('Knockoff')
    guarantees.append((False, False, False, False, False, False))

    R = sim_results['knockoff_plus_R']
    V = sim_results['knockoff_plus_V']
    results.append((np.mean(sim_results['knockoff_plus_screen']),
                    np.nan,
                    np.nan,
                    np.nan,
                    np.mean(V * 1. / np.maximum(R, 1)),
                    np.mean(R - V)))
    names.append('Knockoff+')
    guarantees.append((False, False, False, False, True, False))

    results = np.array(results)
    guarantees = np.array(guarantees)
    error_df = pd.DataFrame({'screen':results[:,0],
                             'fwer.mod':results[:,1],
                             'fwer.mod.cond':results[:,2],
                             'fdr.mod':results[:,3],
                             'fdr.var':results[:,4],
                             's.var':results[:,5]},
                            index=names)

    error_df = error_df.reindex_axis(['screen', 'fwer.mod', 'fwer.mod.cond', 'fdr.mod', 'fdr.var', 's.var'], axis=1)
    file('../../tables/error_rates_%s.html' % alpha, 'w').write(error_df.to_html(float_format = lambda v: '%0.3f' % v))

    def table_generator():
        for result, guar, name in zip(results, guarantees, error_df.index):
            val = [name]
            for _r, g in zip(result, guar):
                if np.isnan(_r):
                    val.append('NA')
                elif not g:
                    val.append('%0.3f' % _r)
                else:
                    val.append(r'\guarantee{%0.3f}' % _r)
            yield ' & '.join(val) + r' \\ '

    table = r'''
\renewcommand{\guarantee}[1]{{\color{blue} #1}}
\begin{tabular}{lrrrrrr}
\toprule
{} &  $p_{\text{screen}}$ &  $\text{FWER}_{\text{mod}}$ &  $\text{FWER}_{\text{mod}} \vert \text{screen}$ &  $\text{FDR}_{\text{model}}$ &  $\text{FDR}_{\text{var}}$ &  $\text{S}_{\text{var}}$ \\ \midrule
%s 
\bottomrule
\end{tabular}''' % '\n'.join(table_generator())

    file('../../tables/error_rates_%s.tex' % alpha, 'w').write(table)

    print table
