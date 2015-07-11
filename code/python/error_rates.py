from itertools import product
import pandas as pd, numpy as np

sim_results = pd.read_csv('../snr_5_alpha_05.csv')

names = []
results = []
error_types = ['pscreen', 'FWER.mod', 'FDR.var', 'FDR.mod']

name_map = dict(zip(['nominal', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                    ['Nominal', 'MaxT', 'MaxT identify', 'MaxT unknown', 'Saturated']))
for test, rule in product(['nominal', 'maxT', 'maxT_identify', 'maxT_unknown', 'saturated'],
                          ['simple', 'forward', 'strong']):
    name = '_'.join([test, rule])
    screen = np.mean(sim_results['_'.join([name, 'screen'])])

    R = sim_results['_'.join([name, 'R'])]
    V_var = sim_results['_'.join([name, 'V_var'])]
    FDR_var = np.mean(V_var * 1. / np.maximum(R, 1))

    V_model = sim_results['_'.join([name, 'V_model'])]
    S_var = np.mean(R - V_var)
    FDR_model = np.mean(V_model * 1. / np.maximum(R, 1))
    FWER_mod = np.mean(V_model > 0)

    results.append((screen, FWER_mod, FDR_var, FDR_model, S_var))
    names.append('%s %s' % (name_map[test], rule))

# now knockoffs

R = sim_results['knockoff_R']
V = sim_results['knockoff_V']
results.append((np.mean(sim_results['knockoff_screen']),
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
                np.mean(V * 1. / np.maximum(R, 1)),
                np.mean(R - V)))
names.append('Knockoff+')

results = np.array(results)
error_df = pd.DataFrame({'screen':results[:,0],
                         'fwer.mod':results[:,1],
                         'fdr.var':results[:,2],
                         'fdr.mod':results[:,3],
                         's.var':results[:,4]},
                        index=names)

error_df = error_df.reindex_axis(['screen', 'fwer.mod', 'fdr.mod', 'fdr.var', 's.var'], axis=1)
error_df.columns = pd.Index(['$p_{\text{screen}}$',
                             '$\text{FWER}_{\text{mod}}$',
                             '$\text{FDR}_{\text{var}}$',
                             '$\text{FDR}_{\text{model}}',
                             '$\text{S}_{\text{var}}'])

file('../../error_rates.tex', 'w').write(error_df.to_latex(float_format = lambda v: '%0.3f' % v))

print error_df.to_latex(float_format = lambda v: '%0.3f' % v)
