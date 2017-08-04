import os
from itertools import product
import pandas as pd, numpy as np

def produce_tables(results, outbase):
    sim_results = pd.read_csv(results)

    getvar = lambda var: sim_results['_'.join([name, var])]

    names = []
    results = []
    stderrs = []
    guarantees = []

    name_map = dict(zip(['nominalT', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                        ['Nominal', 'Max-z', 'Max-z-identify', 'Max-t', 'Saturated']))
    rule_map = dict(zip(['simple', 'forward'],
                        ['BasicStop', 'ForwardStop']))
    for test, rule in product(['nominalT', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                              ['simple', 'forward']):
        name = '_'.join([test, rule])
        n_sim = len(getvar("screen"))
        screen = np.mean(getvar("screen"))
        screen_se = np.sqrt(np.var(getvar("screen"))/n_sim)
        FDR_var = np.mean(getvar("FDP_var"))
        FDR_var_se = np.sqrt(np.var(getvar("FDP_var"))/n_sim)
        FDR_model = np.mean(getvar("FDP_model"))
        FDR_model_se = np.sqrt(np.var(getvar("FDP_model"))/n_sim)
        FDR_model_cond = np.mean(getvar("FDP_model") * getvar("screen")) / screen
        FDR_model_cond_se = np.nan #np.sqrt(np.mean((getvar("FDP_model") * getvar("screen") - FDR_model_cond)^2) / screen / (screen*n_sim))
        FWER_model = np.mean(getvar("FWER_model"))
        FWER_model_se = np.sqrt(np.var(getvar("FWER_model"))/n_sim)
        FWER_model_cond = FWER_model / screen
        FWER_model_cond_se = np.sqrt(FWER_model_cond * (1 - FWER_model_cond) / (n_sim * screen))
        S_var = np.mean(getvar("S_var"))
        S_var_se = np.sqrt(np.var(getvar("S_var"))/n_sim)

        if test in ['maxT', 'maxT_identify', 'maxT_unknown']:
            if rule == 'simple':
                guarantees.append((False, True, True, True, False, False))
            elif rule == 'strong':
                guarantees.append((False, True, False, True, False, False))
            else:
                guarantees.append((False, False, False, True, False, False))
        else:
            guarantees.append((False, False, False, False, False, False))
        results.append((screen, FWER_model, FWER_model_cond, FDR_model, FDR_var, S_var))
        stderrs.append((screen_se, FWER_model_se, FWER_model_cond_se, FDR_model_se, FDR_var_se, S_var_se))
        names.append('%s & %s' % (name_map[test], rule_map[rule]))

    # now knockoffs

    if 'knockoff_R' in sim_results:

        R = sim_results['knockoff_R']
        V = sim_results['knockoff_V']
        knock_FDP = V * 1. / np.maximum(R, 1)
        knock_S_var = R - V
        results.append((np.mean(sim_results['knockoff_screen']),
                        np.nan,
                        np.nan,
                        np.nan,
                        np.mean(knock_FDP),
                        np.mean(knock_S_var)))
        stderrs.append((np.sqrt(np.var(sim_results['knockoff_screen'])/n_sim),
                        np.nan,
                        np.nan,
                        np.nan,
                        np.sqrt(np.var(knock_FDP)/n_sim),
                        np.sqrt(np.var(knock_S_var)/n_sim)))
        names.append('Knockoff &')
        guarantees.append((False, False, False, False, False, False))

        R = sim_results['knockoff_plus_R']
        V = sim_results['knockoff_plus_V']
        knock_FDP = V * 1. / np.maximum(R, 1)
        knock_S_var = R - V
        results.append((np.mean(sim_results['knockoff_plus_screen']),
                        np.nan,
                        np.nan,
                        np.nan,
                        np.mean(knock_FDP),
                        np.mean(knock_S_var)))
        stderrs.append((np.sqrt(np.var(sim_results['knockoff_screen'])/n_sim),
                        np.nan,
                        np.nan,
                        np.nan,
                        np.sqrt(np.var(knock_FDP)/n_sim),
                        np.sqrt(np.var(knock_S_var)/n_sim)))
        names.append('Knockoff+ &')
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
    file('%s.html' % outbase, 'w').write(error_df.to_html(float_format = lambda v: '%0.3f' % v))

    def table_generator():
        for result, stderr, guar, name in zip(results, stderrs, guarantees, error_df.index):
            if name in ['Max-z-identify & BasicStop', 'Max-z-identify & ForwardStop']:
                continue
            val = [name]
            for _r, se, lse, g in zip(result, stderr, largest_stderr, guar):
                if np.isnan(_r):
                    val.append('---')
                elif not g:
                    if lse < 0.01:
                        res_part = r'%0.3f' % _r
                    else:
                        res_part = r'%0.2f' % _r
                    val.append(res_part)
                    #se_part = r'(%0.3f)' % se
                    #val.append(' '.join([res_part, se_part]))
                else:
                    if lse < 0.01:
                        res_part = r'\guarantee{%0.3f}' % _r
                    else:
                        res_part = r'\guarantee{%0.2f}' % _r
                    val.append(res_part)
                    #se_part = r'(%0.3f)' % se
                    #val.append(' '.join([res_part, se_part]))
            yield ' & '.join(val) + r' \\ '
            
    largest_stderr = stderrs[1]
    for stderr in stderrs:
        largest_stderr = [max(l, s) for l, s in zip(largest_stderr, stderr)]

    table = r'''
\newcommand{\guarantee}[1]{{\bf #1}}
{\small 
\begin{tabular}{|l l|cccccc|}
\hline
{} & {} &  $\P(\hk \geq k_0)$ &  $\text{FWER}$ &  $\text{cFWER}$ &  $\text{FDR}$ &  $\text{FDR}^{\text{full}}$ &  $\E[S^{\text{full}}]$ \\ \hline
%s \hline
\end{tabular}}''' % '\n'.join(table_generator())

    file('%s.tex' % outbase, 'w').write(table)

    print table
    print largest_stderr

if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser(
        description= '''
Produce LaTeX and HTML tables from simulation results.
        ''')
    parser.add_argument('--results',
                        help='CSV file containing the results.')
    parser.add_argument('--outbase',
                        help='Base for output filenames.')

    args = parser.parse_args()
    produce_tables(args.results, args.outbase)

        #  zip(['../snr_5_alpha_%s.csv' % alpha for alpha in ['05', '10', '20']] + ['../snr_7_alpha_20_sparsity7_p200.csv'],
        #                             ['../../tables/error_rates_%s' % alpha for alpha in ['05', '10', '20']] + ['../../error_rates_p200'])[:-1]:

