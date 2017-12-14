import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import random, numpy, scipy
from scipy.stats import * # pareto, expon, norm

from patch import *
from rvs import *

def plot_pareto():
  a = 100
  dist = pareto(a)
  # dist = expon()
  # dist = norm(1)
  # dist = uniform()
  """
  # x = numpy.linspace(dist.ppf(0.01), dist.ppf(0.99), 100)
  x = numpy.linspace(1, 5, 100)
  plot.plot(x, dist.cdf(x), 'r-', lw=5, alpha=0.6, label='pareto pdf')
  """
  def D_pareto(d, u_l):
    # return [(1-d)*(1-u)**(-1/a) + d for u in u_l]
    def F(u):
      return 1 - u**(-a)
    def F_inv(u):
      return (1 - u)**(-1/a)
    return [F_inv(u) - F_inv(u+F(d)*(1-u) ) + d for u in u_l]
  
  color = iter(cm.rainbow(numpy.linspace(0, 1, 10) ) )
  # for d in numpy.arange(0.1, 1, 0.1):
  for d in numpy.arange(0.1, 10, 1):
    u_l = numpy.linspace(0.01, 0.99, 100)
    D_l = [dist.ppf(u) - dist.ppf(u+dist.cdf(d)*(1-u) ) + d for u in u_l]
    D_l_ = D_pareto(d, u_l)
    
    c = next(color)
    plot.plot(u_l, D_l, 'o', lw=5, color=c, alpha=0.2, label=r'd= {0:.2f}'.format(d) )
    # plot.plot(u_l, D_l_, '.', lw=5, color=c, alpha=0.5, label=r'd= {0:.2f}'.format(d) )
  plot.legend()
  plot.ylabel(r'$x$')
  plot.ylabel(r'$D_r(\Delta)$')
  plot.savefig("plot_pareto.png")
  plot.gcf().clear()

# #####################  (k, d)  ################### #
def sim_arepeat_k_c(task_t_rv, d, k, c, num_run, w_relaunch=False):
  stat_id__trial_stat_l_m = {'T': [], 'C': [], 'C_wc': [], 'T_2': [], 'C_2': [] }
  for i in range(num_run):
    i__t_l_m = {i:[[0, task_t_rv.gen_sample() ] ] for i in range(k) }
    job_compl_t = max([t_l[0][1] for i,t_l in i__t_l_m.items() ] )
    if job_compl_t > d:
      for i,t_l in i__t_l_m.items():
        if t_l[0][1] > d:
          if w_relaunch:
            t_l[0][1] = d + task_t_rv.gen_sample()
          for i in range(c):
            t_l.append([d, d + task_t_rv.gen_sample() ] )
    
    E_T = max([min([t[1] for t in t_l] ) for i,t_l in i__t_l_m.items() ] )
    stat_id__trial_stat_l_m['T'].append(E_T)
    stat_id__trial_stat_l_m['T_2'].append(E_T**2)
    stat_id__trial_stat_l_m['C'].append(
      sum([sum([t[1]-t[0] for t in t_l] ) for i,t_l in i__t_l_m.items() ] ) )
    E_C_wc = 0
    for i,t_l in i__t_l_m.items():
      compl_t = min([t[1] for t in t_l] )
      for t in t_l:
        E_C_wc += min(compl_t, t[1])-t[0]
    stat_id__trial_stat_l_m['C_wc'].append(E_C_wc)
    stat_id__trial_stat_l_m['C_2'].append(E_C_wc**2)
  
  return stat_id__trial_stat_l_m

# ##################  (l, k, n, d)  ################ #
def sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run, w_relaunch=False):
  if l < k:
    l = k
  log(DEBUG, "k= {}, l= {}, n= {}".format(k, l, n) )
  # if l < k or n < l:
  #   return 1
  stat_id__trial_stat_l_m = {'T': [], 'C': [], 'C_wc': [], 'T_2': [], 'C_2': [] }
  for i in range(num_run):
    compl_t_l = [[0, task_t_rv.gen_sample() ] for i in range(l) ]
    compl_t_l.sort(key=lambda tup: tup[1] )
    if compl_t_l[k-1][1] > d:
      if w_relaunch:
        for t in compl_t_l:
          if t[1] > d: t[1] = d + task_t_rv.gen_sample()
      compl_t_l += [[d, d + task_t_rv.gen_sample() ] for i in range(n-l) ]
      compl_t_l.sort(key=lambda tup: tup[1] )
    
    E_T = compl_t_l[k-1][1]
    stat_id__trial_stat_l_m['T'].append(E_T)
    stat_id__trial_stat_l_m['T_2'].append(E_T**2)
    stat_id__trial_stat_l_m['C'].append(sum([t[1]-t[0] for t in compl_t_l] ) )
    E_C = sum([min(t[1], E_T)-t[0] for t in compl_t_l] )
    stat_id__trial_stat_l_m['C_wc'].append(E_C)
    stat_id__trial_stat_l_m['C_2'].append(E_C**2)
  
  return stat_id__trial_stat_l_m

# ##################  (l, k, n, d)  ################ #
def sim_arepeat_k_nd0_wrelaunch(task_t_rv, k, n, d, num_run):
  log(DEBUG, "k= {}, n= {}, d= {}".format(k, n, d) )
  stat_id__trial_stat_l_m = {'T': [], 'C': [], 'C_wc': [], 'T_2': [], 'C_2': [] }
  for i in range(num_run):
    compl_t_l = [task_t_rv.gen_sample() for i in range(n) ]
    compl_t_l.sort()
    if compl_t_l[k-1] > d:
      for i, t in enumerate(compl_t_l):
        if t > d: compl_t_l[i] = d + task_t_rv.gen_sample()
      compl_t_l.sort()
    
    E_T = compl_t_l[k-1]
    stat_id__trial_stat_l_m['T'].append(E_T)
    stat_id__trial_stat_l_m['T_2'].append(E_T**2)
    stat_id__trial_stat_l_m['C'].append(sum([t for t in compl_t_l] ) )
    E_C = sum([min(t, E_T) for t in compl_t_l] )
    stat_id__trial_stat_l_m['C_wc'].append(E_C)
    stat_id__trial_stat_l_m['C_2'].append(E_C**2)
  
  return stat_id__trial_stat_l_m

if __name__ == "__main__":
  plot_pareto()
  