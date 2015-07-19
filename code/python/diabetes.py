import numpy as np
import pandas as pd
from rules import forward_stop
from selection.algorithms.forward_step import forward_step
from compute_pvalues import compute_pvalues

import rpy2.robjects as rpy
from rpy2.robjects import numpy2ri
rpy.conversion.py2ri = numpy2ri.numpy2ri

numpy2ri.activate()
rpy.r('''library(lars)
library(knockoff)
data(diabetes)
X = diabetes$x2
Y = diabetes$y
D.lm = lm(Y ~ X)
sigma = sqrt(sum(resid(D.lm)^2 / D.lm$df.resid))
knock_plus = knockoff.filter(X = X, y = Y, fdr=0.1, threshold="knockoff+")$selected
print('knockoff+')
print(knock_plus)

knock = knockoff.filter(X = X, y = Y, fdr=0.1, threshold="knockoff")$selected
print('knockoff')
print(knock)

''')

X = np.array(rpy.r('X'))
n, p = X.shape

Y = np.array(rpy.r('Y'))
varnames = list(np.array(rpy.r('colnames(X)')))
sigma = np.array(rpy.r('sigma'))[0]

pvals = compute_pvalues(Y, X, sigma=sigma, 
                        maxstep=55,
                        compute_maxT_identify=False,
                        burnin=10000,
                        ndraw=40000,
                        accept_reject_params=(100,15,10000),
                        shortcut=True)[0]

pvals['Variable'] = [varnames[i] for i in pvals['variable_selected']]
pvals['Step'] = np.arange(pvals.shape[0]) + 1
pvals['variable_selected'] += 1
pvals = pvals.reindex_axis(['Step', 'Variable', 'variable_selected',
                            'nominalT_pvalue',
                            'saturated_pvalue',
                            'maxT_unknown_pvalue',
                            'maxT_pvalue'], axis=1)
pvals.columns = pd.Index(['Step', 'Variable', 'Column number', 'Nominal pvalue', 'Saturated pvalue', 'MaxT pvalue', 'MaxT pvalue plugin sigma'])

forward_stop_U = forward_stop(pvals['MaxT pvalue'], 0.1)
forward_stop_N = forward_stop(pvals['Nominal pvalue'], 0.1)
forward_stop_S = forward_stop(pvals['Saturated pvalue'], 0.1)

print 'maxT unknown:', forward_stop_U
print 'nominal:',  forward_stop_N
print 'saturated:', forward_stop_S

pvals = pvals[:20]

# R pvalues

Rpval = []
model_str = ''
for i in range(pvals.shape[0]):
    model_str = '+'.join([' X[,%d] ' % v for v in pvals['Column number'][:(i+1)]])
    Rstr = 'summary(lm(Y ~ %s))$coef[,4]' % model_str
    Rpval.append(np.array(rpy.r(Rstr))[-1])

print 'checking whether nominal agrees with R:', np.linalg.norm(np.array(Rpval) - pvals['Nominal pvalue']) / np.linalg.norm(pvals['Nominal pvalue']) 

# save the HTML table

file('../../tables/diabetes.html', 'w').write(pvals.to_html(float_format = lambda v : '%0.2f' % v, index=False))

pvals = pvals.reindex_axis(['Step', 'Variable', 'Nominal pvalue', 'Saturated pvalue', 'MaxT pvalue'], axis=1)
print pvals

# save the LaTeX table

file('../../tables/diabetes.tex', 'w').write(pvals.to_latex(float_format = lambda v : '%0.2f' % v, index=False))

numpy2ri.deactivate()
