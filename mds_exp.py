import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt, itertools

from mds_sim import *
from mds_models import *

def sim_mds_nk(num_frun, ar, n, k, dist_m, r=None, preempt=False, fi_l=[] ):
  E_T_sum = 0
  for f in range(num_frun):
    log(WARNING, "ar= {}, n= {}, k= {}, dist_m= {}".format(ar, n, k, dist_m) )
    env = simpy.Environment()
    pg = MDS_PG(env, "pg", ar)
    mdsq = MDSQ("mdsq", env, k, range(n), dist_m)
    # monitor = MDSQMonitor(env, mdsq, lambda: 1)
    pg.out = mdsq
    env.run(until=50000)
    
    st_l = mdsq.jsink.st_l
    E_T_sum += float(sum(st_l) )/len(st_l)
    
    total_num_wins = sum([n for i, n in mdsq.jsink.qid__num_win_map.items() ] )
    qid__win_freq_map = {i:float(n)/total_num_wins for i, n in mdsq.jsink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
    
    total_n_types = sum(mdsq.servtype__num_m)
    fi_l[:] = [n/total_n_types for t,n in enumerate(mdsq.servtype__num_m) ]
    print("fi_l= {}".format(fi_l) )
    
    # print("\n")
    # # print("monitor.state__counter_map= {}".format(pprint.pformat(monitor.state__counter_map) ) )
    # total_counter = sum([c for rs, c in monitor.state__counter_map.items() ] )
    # polled_state__freq_map = {rs:float(c)/total_counter for rs, c in monitor.state__counter_map.items() }
    # print("polled_state__freq_map= {}".format(pprint.pformat(polled_state__freq_map) ) )
    # print("----------------------------------------")
  E_T = E_T_sum/num_frun
  if E_T > 100: return None
  return E_T

def plot_mds_n_2():
  n = 5 # 10 # 5
  k = 2 # 6 # 4
  # dist_m = {'dist': 'Exp', 'mu': 1}
  # dist_m = {'dist': 'SExp', 'D': 1, 'mu': 1}
  dist_m = {'dist': 'Pareto', 'l': 1, 'a': 2}
  ar_ub = mds_exact_bound_on_arr_rate(n, k, dist_m) # mds_innerbound_on_ar(n, k, dist_m)
  log(WARNING, "n= {}, k= {}, dist_m={}".format(n, k, dist_m) )
  
  ar_l = []
  ET_sm_l, ET_sim_l, ET_approx_l = [], [], []
  ET_approx2_l = []
  
  num_frun = 1
  for ar in numpy.linspace(0.05, ar_ub, 10):
    pi_l = []
    ET_sim = sim_mds_nk(num_frun, ar, n, k, dist_m, fi_l=pi_l)
    print("ET_sim= {}".format(ET_sim) )
    if ET_sim is None: break
    ar_l.append(ar)
    ET_sim_l.append(ET_sim)
    
    ET_sm = ET_mds_nk_sm(ar, n, k, dist_m)
    print("ET_sm= {}".format(ET_sm) )
    ET_sm_l.append(ET_sm)
    
    ET_approx = ET_mds_n2_approx(ar, n, dist_m)
    print("ET_approx= {}".format(ET_approx) )
    ET_approx_l.append(ET_approx)
    
    ET_approx_ = ET_mds_nk_approx(ar, n, k, dist_m, pi_l)
    print("ET_approx_= {}".format(ET_approx_) )
    ET_approx2_l.append(ET_approx_)
  plot.plot(ar_l, ET_sm_l, color=next(dark_color), label='Split-merge', marker=next(marker), linestyle=':', mew=2)
  # print("ET_sim_l= {}".format(pprint.pformat(ET_sim_l) ) )
  plot.plot(ar_l, ET_sim_l, color=next(dark_color), label='Simulation', marker=next(marker), linestyle=':', mew=2)
  plot.plot(ar_l, ET_approx_l, color=next(dark_color), label='Approx', marker=next(marker), linestyle=':', mew=2)
  # plot.plot(ar_l, ET_approx2_l, color=next(dark_color), label='Approx_', marker=next(marker), linestyle=':', mew=2)
  
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel(r'$E[T]$')
  plot.title(r'$n= {}$, $k= {}$, $V \sim {}$'.format(n, k, dist_m) )
  plot.savefig("plot_mds_{}_{}.png".format(n, k) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_mds_n_2()
  
