import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt, itertools

from mds_sim import *
from mds_models import *

def sim_mds_nk(num_frun, ar, n, k, dist_m, r=None, preempt=False, fi_l=[] ):
  ET_sum = 0
  for f in range(num_frun):
    log(WARNING, "ar= {}, n= {}, k= {}, dist_m= {}".format(ar, n, k, dist_m) )
    env = simpy.Environment()
    pg = MDS_PG(env, "pg", ar)
    mdsq = MDSQ("mdsq", env, k, range(n), dist_m)
    # monitor = MDSQMonitor(env, mdsq, lambda: 1)
    pg.out = mdsq
    env.run(until=50000*10)
    
    st_l = mdsq.jsink.st_l
    ET_sum += float(sum(st_l) )/len(st_l)
    
    total_num_wins = sum([n for i, n in mdsq.jsink.qid__num_win_map.items() ] )
    qid__win_freq_map = {i:float(n)/total_num_wins for i, n in mdsq.jsink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
    
    total_numtypes = sum(mdsq.servtype__num_m)
    fi_l[:] = [n/total_numtypes for t,n in enumerate(mdsq.servtype__num_m) ]
    print("fi_l= {}".format(fi_l) )
    
    # print("\n")
    # # print("monitor.state__counter_map= {}".format(pprint.pformat(monitor.state__counter_map) ) )
    # total_counter = sum([c for rs, c in monitor.state__counter_map.items() ] )
    # polled_state__freq_map = {rs:float(c)/total_counter for rs, c in monitor.state__counter_map.items() }
    # print("polled_state__freq_map= {}".format(pprint.pformat(polled_state__freq_map) ) )
    # print("----------------------------------------")
  ET = ET_sum/num_frun
  if ET > 100: return None
  return ET

def plot_mds_n_2():
  n = 4 # 10 # 5
  k = 2 # 6 # 4
  # dist_m = {'dist': 'Exp', 'mu': 1}
  dist_m = {'dist': 'SExp', 'D': 1, 'mu': 1}
  # dist_m = {'dist': 'Pareto', 'l': 1, 'a': 2}
  
  if dist_m['dist'] == 'Exp':
    ar_ub = 0.9*mds_exact_bound_on_arr_rate(n, k, dist_m) # mds_innerbound_on_ar(n, k, dist_m)
  else:
    ar_ub = 0.85*mds_exact_bound_on_arr_rate(n, k, dist_m)
  log(WARNING, "n= {}, k= {}, dist_m={}".format(n, k, dist_m) )
  
  ar_l = []
  ET_sm_l, ET_sim_l, ET_approx_l = [], [], []
  ET_approx2_l = []
  ET_varkigauri_lb_l = []
  
  num_frun = 1
  for ar in [*numpy.linspace(0.05, 0.6*ar_ub, 4, endpoint=False), *numpy.linspace(0.6*ar_ub, ar_ub, 7) ]:
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
    
    # ET_approx_ = ET_mds_nk_approx(ar, n, k, dist_m, pi_l)
    # print("ET_approx_= {}".format(ET_approx_) )
    # ET_approx2_l.append(ET_approx_)
    
    ET_varkigauri_lb = ET_mds_nk_varkigauri_lb(ar, n, k, dist_m)
    print("ET_varkigauri_lb= {}".format(ET_varkigauri_lb) )
    ET_varkigauri_lb_l.append(ET_varkigauri_lb)
  plot.plot(ar_l, ET_sm_l, color=next(dark_color), label='Split-merge upper bound', marker=next(marker), linestyle=':', mew=2)
  # print("ET_sim_l= {}".format(pprint.pformat(ET_sim_l) ) )
  plot.plot(ar_l, ET_sim_l, color=next(dark_color), label='Simulation', marker=next(marker), linestyle=':', mew=2)
  plot.plot(ar_l, ET_approx_l, color=next(dark_color), label=r'$M/G/1$ approximation', marker=next(marker), linestyle=':', mew=2)
  # plot.plot(ar_l, ET_approx2_l, color=next(dark_color), label='Approx_', marker=next(marker), linestyle=':', mew=2)
  if dist_m['dist'] == 'Exp':
    plot.plot(ar_l, ET_varkigauri_lb_l, color=next(dark_color), label='Varki-Gauri lower bound', marker=next(marker), linestyle=':', mew=2)
  plot.legend()
  plot.xlabel(r'$\lambda$', fontsize=14)
  plot.ylabel(r'$E[T]$', fontsize=14)
  plot.title(r'$n= {}$, $k= {}$, $V \sim {}$'.format(n, k, dist_to_latex(dist_m) ) )
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_mds_{}_{}_{}.pdf".format(n, k, dist_m['dist'] ) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_mds_n_2()
  
