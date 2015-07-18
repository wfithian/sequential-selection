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
#X = scale(X, T, T)
Y = diabetes$y
#Y = Y - mean(Y)
D.lm = lm(Y ~ X)
sigma = sqrt(sum(resid(D.lm)^2 / D.lm$df.resid))
knock = knockoff.filter(X = X, y = Y, fdr=0.1, threshold="knockoff+")$selected
print(knock)
''')
X = np.array(rpy.r('X'))
n, p = X.shape

Y = np.array(rpy.r('Y'))
sigma = np.array(rpy.r('sigma'))[0]

FS_maxT_U = forward_step(X, Y, covariance=sigma**2 * np.identity(Y.shape[0])) # not really using covariance here...
iter(FS_maxT_U)

FS_maxT = forward_step(X, Y, covariance=sigma**2 * np.identity(Y.shape[0])) 
iter(FS_maxT)

FS_maxT2 = forward_step(X, Y, covariance=sigma**2 * np.identity(Y.shape[0])) 
iter(FS_maxT2)

results = []

for i in range(64):
    pval_maxT = FS_maxT.next(compute_pval=True,
                             use_identity=False,
                             ndraw=40000,
                             burnin=10000,
                             sigma_known=True,
                             accept_reject_params=(100, 15, 10000))
    pval_maxT_U = FS_maxT_U.next(compute_pval=True,
                                 use_identity=False,
                                 ndraw=40000,
                                 burnin=10000,
                                 sigma_known=False)

    X_cur = X.copy()[:,FS_maxT.variables]
    X_cur -= X_cur.mean(0)[None,:]
    P_cur = np.dot(X_cur, np.linalg.pinv(X_cur))
    sigma_cur = np.linalg.norm(Y - np.dot(P_cur, Y) - Y.mean()) / np.sqrt(n - np.diag(P_cur).sum())

    FS_maxT2.covariance[:] = sigma_cur**2 * np.identity(n)
    pval_maxT2 = FS_maxT2.next(compute_pval=True,
                               use_identity=False,
                               ndraw=40000,
                               burnin=10000,
                               sigma_known=True,
                               accept_reject_params=())


    results.append((pval_maxT, pval_maxT2, pval_maxT_U, sigma_cur))
    print results[-1]

results = np.array(results)
pval_table = pd.DataFrame({'known':results[:,0],
                           'unknown':results[:,1],
                           'using.sigma.cur':results[:,2],
                           'sigma.cur':results[:,3]})

file('diabetes.html', 'w').write(pval_table.to_html(float_format = lambda v : '%0.3f' % v))

forward_stop_U = forward_stop(results[:,1], 0.1)
forward_stop_K = forward_stop(results[:,0], 0.1)
forward_stop_K2 = forward_stop(results[:,2], 0.1)

print forward_stop_U, forward_stop_K, forward_stop_K2


