import numpy as np, os
import rpy2.robjects as rpy
from rpy2.robjects import numpy2ri
rpy.conversion.py2ri = numpy2ri.numpy2ri
from matplotlib import mlab

from selection.algorithms.tests.test_forward_step import test_full_pvals

numpy2ri.activate()
rpy.r("""
strongstop <- function(p.values,alpha) {
   d <- length(p.values)
   lhs <- exp(rev(cumsum(rev(log(p.values)/(1:d))))) # LHS from G'Sell et al.
   rhs <- alpha * (1:d) / d # RHS from G'Sell et al.
   return(max(c(0,which(lhs <= rhs))))
}

forwardstop <- function(p, alpha) {
   m <- length(p)
   sums <- -(1/(1:m))*cumsum(log(1-p))
   return(max(c(0, which(sums < alpha))))
}
"""
)
numpy2ri.deactivate()

def simple_stop(pvalues, alpha):
    if not np.all(pvalues <= alpha):
        return np.min(np.nonzero(pvalues > alpha)[0])
    else:
        return pvalues.shape[0]

def strong_stop(pvalues, alpha):
    numpy2ri.activate()
    rpy.r.assign('P', pvalues)
    rpy.r.assign('alpha', alpha)
    rpy.r('IDX = strongstop(P, alpha)')
    value = int(np.array(rpy.r('IDX'))[0])
    numpy2ri.deactivate()
    return value

def forward_stop(pvalues, alpha):
    numpy2ri.activate()
    rpy.r.assign('P', pvalues)
    rpy.r.assign('alpha', alpha)
    rpy.r('IDX = forwardstop(P, alpha)')
    value = int(np.array(rpy.r('IDX'))[0])
    numpy2ri.deactivate()
    return value

def run_one(n=100, p=40, rho=0.3, snr=4, alpha=0.05):
    (X, 
     y, 
     beta, 
     active, 
     sigma, 
     results, 
     completion_index) = test_full_pvals(n=n,
                                         p=p,
                                         snr=snr,
                                         rho=rho)
    numpy2ri.activate()
    rpy.r.assign('X', X)
    rpy.r.assign('y', y)

    # knockoff

    rpy.r.assign('alpha', alpha)
    knockoff = np.array(rpy.r("""
library(knockoff)
knockoff.filter(X = X, y = y, fdr=alpha)$selected
"""))
    knockoff_R = knockoff.shape[0]
    knockoff_V = knockoff_R - len(set(active).intersection(knockoff))

    # simple_stop

    P_selected = results[:,1]
    simple_selected = simple_stop(P_selected, alpha)
    simple_selected_R = simple_selected 
    simple_selected_V_var = simple_selected_R - len(set(active).intersection(results[:simple_selected_R,0]))
    simple_selected_V_model = max(simple_selected_R - completion_index, 0)
    simple_selected_screen = simple_selected_R >= completion_index

    P_saturated = results[:,2]
    simple_saturated = simple_stop(P_saturated, alpha)
    simple_saturated_R = simple_saturated 
    simple_saturated_V_var = simple_saturated_R - len(set(active).intersection(results[:simple_saturated_R,0]))
    simple_saturated_V_model = max(simple_selected_R - completion_index, 0)
    simple_saturated_screen = simple_saturated_R >= completion_index

    P_nominal = results[:,3]
    simple_nominal = simple_stop(P_nominal, alpha)
    simple_nominal_R = simple_nominal 
    simple_nominal_V_var = simple_nominal_R - len(set(active).intersection(results[:simple_nominal_R,0]))
    simple_nominal_V_model = max(simple_selected_R - completion_index, 0)
    simple_nominal_screen = simple_nominal_R >= completion_index

    # strong_stop

    P_selected = results[:,1]
    strong_selected = strong_stop(P_selected, alpha)
    strong_selected_R = strong_selected 
    strong_selected_V_var = strong_selected_R - len(set(active).intersection(results[:strong_selected_R,0]))
    strong_selected_V_model = max(strong_selected_R - completion_index, 0)
    strong_selected_screen = strong_selected_R >= completion_index

    P_saturated = results[:,2]
    strong_saturated = strong_stop(P_saturated, alpha)
    strong_saturated_R = strong_saturated
    strong_saturated_V_var = strong_saturated_R - len(set(active).intersection(results[:strong_saturated_R,0]))
    strong_saturated_V_model = max(strong_saturated_R - completion_index, 0)
    strong_saturated_screen = strong_saturated_R >= completion_index

    P_nominal = results[:,3]
    strong_nominal = strong_stop(P_nominal, alpha)
    strong_nominal_R = strong_nominal 
    strong_nominal_V_var = strong_nominal_R - len(set(active).intersection(results[:strong_nominal_R,0]))
    strong_nominal_V_model = max(strong_nominal_R - completion_index, 0)
    strong_nominal_screen = strong_nominal_R >= completion_index

    # forward_stop

    P_selected = results[:,1]
    forward_selected = forward_stop(P_selected, alpha)
    forward_selected_R = forward_selected 
    forward_selected_V_var = forward_selected_R - len(set(active).intersection(results[:forward_selected_R,0]))
    forward_selected_V_model = max(forward_selected_R - completion_index, 0)
    forward_selected_screen = forward_selected_R >= completion_index

    P_saturated = results[:,2]
    forward_saturated = forward_stop(P_saturated, alpha)
    forward_saturated_R = forward_saturated
    forward_saturated_V_var = forward_saturated_R - len(set(active).intersection(results[:forward_saturated_R,0]))
    forward_saturated_V_model = max(forward_saturated_R - completion_index, 0)
    forward_saturated_screen = forward_saturated_R >= completion_index

    P_nominal = results[:,3]
    forward_nominal = forward_stop(P_nominal, alpha)
    forward_nominal_R = forward_nominal 
    forward_nominal_V_var = forward_nominal_R - len(set(active).intersection(results[:forward_nominal_R,0]))
    forward_nominal_V_model = max(forward_nominal_R - completion_index, 0)
    forward_nominal_screen = forward_nominal_R >= completion_index

    results = (list(results.T.reshape(-1)) + 
               list(active) + 
               [completion_index] + 
               [simple_saturated_V_var,
                simple_saturated_V_model,
                simple_saturated_R,
                simple_saturated_screen,
                simple_selected_V_var,
                simple_selected_V_model,
                simple_selected_R,
                simple_selected_screen,
                simple_nominal_V_var,
                simple_nominal_V_model,
                simple_nominal_R,
                simple_nominal_screen] +
               [strong_saturated_V_var,
                strong_saturated_V_model,
                strong_saturated_R,
                strong_saturated_screen,
                strong_selected_V_var,
                strong_selected_V_model,
                strong_selected_R,
                strong_selected_screen,
                strong_nominal_V_var,
                strong_nominal_V_model,
                strong_nominal_R,
                strong_nominal_screen] +
               [forward_saturated_V_var,
                forward_saturated_V_model,
                forward_saturated_R,
                forward_saturated_screen,
                forward_selected_V_var,
                forward_selected_V_model,
                forward_selected_R,
                forward_selected_screen,
                forward_nominal_V_var,
                forward_nominal_V_model,
                forward_nominal_R,
                forward_nominal_screen] +
               [knockoff_V, knockoff_R])

    results = tuple([int(v) for v in results[:p]] + results[p:])
    dtype = ([('var_%d' % i, np.int) for i in range(1, p+1)] +
             [('select_%d' % i, np.float) for i in range(1, p+1)] +
             [('saturated_%d' % i, np.float) for i in range(1, p+1)] +
             [('nominal_%d' % i, np.float) for i in range(1, p+1)] + 
             [('active_%d' % i, np.int) for i in range(1, active.shape[0]+1)]
             + [('completion_index', np.int),
                ('simple_saturated_V_var', np.int), 
                ('simple_saturated_V_model', np.int), 
                ('simple_saturated_R', np.int),
                ('simple_saturated_screen', np.int),
                ('simple_selected_V_var', np.int), 
                ('simple_selected_V_model', np.int), 
                ('simple_selected_R', np.int),
                ('simple_selected_screen', np.int),
                ('simple_nominal_V_var', np.int), 
                ('simple_nominal_V_model', np.int), 
                ('simple_nominal_R', np.int),
                ('simple_nominal_screen', np.int)]
             + [('strong_saturated_V_var', np.int), 
                ('strong_saturated_V_model', np.int), 
                ('strong_saturated_R', np.int),
                ('strong_saturated_screen', np.int),
                ('strong_selected_V_var', np.int), 
                ('strong_selected_V_model', np.int), 
                ('strong_selected_R', np.int),
                ('strong_selected_screen', np.int),
                ('strong_nominal_V_var', np.int), 
                ('strong_nominal_V_model', np.int), 
                ('strong_nominal_R', np.int),
                ('strong_nominal_screen', np.int)]
             + [('forward_saturated_V_var', np.int), 
                ('forward_saturated_V_model', np.int), 
                ('forward_saturated_R', np.int),
                ('forward_saturated_screen', np.int),
                ('forward_selected_V_var', np.int), 
                ('forward_selected_V_model', np.int), 
                ('forward_selected_R', np.int),
                ('forward_selected_screen', np.int),
                ('forward_nominal_V_var', np.int), 
                ('forward_nominal_V_model', np.int), 
                ('forward_nominal_R', np.int),
                ('forward_nominal_screen', np.int)]
             + [('knockoff_V', np.int), 
                ('knockoff_R', np.int)])
    value = np.array([results], dtype)
    numpy2ri.deactivate()

    return value


def batch(fbase, nsim=100, n=100, p=40, rho=0.3, snr=4):
    value = []
    for _ in range(nsim):
        try:
            value.append(run_one(n=n, p=p, rho=rho, snr=snr, alpha=0.05))
        except:
            pass
        if not os.path.exists(fbase + '.npy'):
            V = np.hstack(value)
            np.save(fbase + '.npy', V)
            mlab.rec2csv(V, fbase + '.csv')
        else:
            V = np.load(fbase + '.npy')
            V = np.hstack([V] + value[-1:])
            np.save(fbase + '.npy', V)
            mlab.rec2csv(V, fbase + '.csv')
            
        print np.mean(V['forward_selected_screen']), np.mean(V['forward_saturated_screen']), V.shape
if __name__ == "__main__":
    batch('test', nsim=1000, p=40, snr=5)
        
