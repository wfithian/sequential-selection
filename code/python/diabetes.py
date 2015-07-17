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
knock = knockoff.filter(X = X, y = Y, fdr=0.1, threshold="knockoff+")$selected
print(knock)
''')
X = np.array(rpy.r('X'))
Y = np.array(rpy.r('Y'))
sigma = np.array(rpy.r('sigma'))[0]

FS_maxT_U = forward_step(X, Y, covariance=sigma**2 * np.identity(Y.shape[0])) # not really using covariance here...
iter(FS_maxT_U)

FS_maxT = forward_step(X, Y, covariance=sigma**2 * np.identity(Y.shape[0])) 
iter(FS_maxT)

results = []

for i in range(64):
    pval_maxT = FS_maxT.next(compute_pval=True,
                             use_identity=False,
                             ndraw=80000,
                             burnin=20000,
                             sigma_known=True,
                             alternative='greater')
    pval_maxT_U = FS_maxT_U.next(compute_pval=True,
                                 use_identity=False,
                                 ndraw=80000,
                                 burnin=20000,
                                 sigma_known=False,
                                 alternative='twosided')
    results.append((pval_maxT, pval_maxT_U))
    print results[-1]

results = np.array(results)
pval_table = pd.DataFrame({'known':results[:,0],
                           'unknown':results[:,1]})

print pval_table.to_html(float_format = lambda v : '%0.3f' % v)

forward_stop_U = forward_stop(results[:,1], 0.1)
forward_stop_U = forward_stop(results[:,0], 0.1)


