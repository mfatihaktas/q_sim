import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt, itertools

from multiq_sim import *

from mds_sim import *
from mds_models import *

def scale_distm(dist_m, k):
  m = dist_m.copy()
  m['l'] /= k
  return m

# Simulating a system of N servers with all arrivals completing as mds[n, k]
def plot_Nserver_mdsnkjobs():
  N = 10
  tsize_dist = DUniform(1, 1)
  dist_m = {'dist': 'Pareto', 'l': 1, 'a': 2}
  
  def sim(ar, k, n, dist_m):
    k_dist = DUniform(k, k)
    sl_dist = rv_from_m(dist_m) # Pareto(l=1, a=2)
    sching_m = {'t': 'coded', 'r': n/k}
    
    env = simpy.Environment()
    jg = JG(env, ar, k_dist, tsize_dist)
    mq = MultiQ(env, N, sching_m, sl_dist)
    jg.out = mq
    jg.init()
    env.run(until=50000*1)
    return sum(mq.jtime_l)/len(mq.jtime_l)
  
  def plot_varyingk(ar, n):
    print("> ar= {}, n= {}".format(ar, n) )
    k_l, ET_l = [], []
    # for k in range(1, n+1):
    for k in range(n, n+1):
      print("k= {}".format(k) )
      k_l.append(k)
      
      ET = sim(ar, k, n, scale_distm(dist_m, k) )
      print("ET= {}".format(ET) )
      ET_l.append(ET)
    plot.plot(k_l, ET_l, color=next(dark_color), label=r'$\lambda= {}$'.format(ar), marker=next(marker), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$k$', fontsize=14)
    plot.title(r'$N= {}$, $n= {}$, $V \sim {}$'.format(N, n, dist_to_latex(dist_m) ) )
  
  def plot_varyingn(ar, k):
    print("> ar= {}, k= {}".format(ar, k) )
    dist_m = scale_distm(dist_m, k)
    n_l, ET_l = [], []
    for n in range(k, N+1):
      print("n= {}".format(n) )
      n_l.append(n)
      
      ET = sim(ar, k, n, dist_m)
      print("ET= {}".format(ET) )
      ET_l.append(ET)
      plot.plot(n_l, ET_l, color=next(dark_color), label=r'$\lambda= {}$'.format(ar), marker=next(marker), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$n$', fontsize=14)
    plot.title(r'$N= {}$, $k= {}$, $V \sim {}$'.format(N, k, dist_to_latex(dist_m) ) )
  
  print("N= {}, dist_m= {}".format(N, dist_m) )
  # plot_varyingk(ar=1.6, n=4)
  # plot_varyingk(ar=3, n=4)
  plot_varyingk(ar=4, n=4)
  
  plot.legend()
  plot.ylabel('Average download time', fontsize=14)
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_N{}_mdsnkjobs.pdf".format(N) )
  log(WARNING, "done; N= {}".format(N) )

# ###########################################  ISIT'18  ########################################## #
# Simulating a system of N/n decoupled mds[n, k] storage
def sim_mdsnk(num_frun, ar, n, k, dist_m, r=None, preempt=False, fi_l=[] ):
  ET_sum = 0
  for f in range(num_frun):
    log(WARNING, "ar= {}, n= {}, k= {}, dist_m= {}".format(ar, n, k, dist_m) )
    env = simpy.Environment()
    pg = MDS_PG(env, "pg", ar)
    mdsq = MDSQ("mdsq", env, k, range(n), dist_m)
    # monitor = MDSQMonitor(env, mdsq, lambda: 1)
    pg.out = mdsq
    env.run(until=50000*1)
    
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
  # if ET > 100: return None
  return ET

def plot_Nnk():
  N = 20
  dist_m = {'dist': 'Pareto', 'l': 10, 'a': 2}
  num_frun = 10 # 10
  
  def plot_varyingk(ar, n):
    print("> ar= {}, n= {}".format(ar, n) )
    ar_ = ar * n/N
    k_l, ET_l = [], []
    for k in range(1, n+1):
      print("k= {}".format(k) )
      k_l.append(k)
      
      ET = sim_mdsnk(num_frun, ar_, n, k, scale_distm(dist_m, k) )
      print("ET= {}".format(ET) )
      ET_l.append(ET)
    plot.plot(k_l, ET_l, color=next(dark_color), label=r'$\lambda= {}$'.format(ar), marker=next(marker), linestyle=':', mew=2)
    plot.xlabel(r'$k$', fontsize=14)
    # plot.title(r'$N= {}$, $n= {}$, $V \sim {}$'.format(N, n, dist_to_latex(dist_m) ) )
    # plot.title(r'$n= {}$, $V \sim {}$'.format(n, dist_to_latex(dist_m) ) )
    plot.title(r'$n= {}$, $S= {}$, Slowdown$\sim Pareto(1, \alpha=2)$'.format(n, dist_m['l'] ) )
  
  def plot_varyingn(ar, k):
    print("> ar= {}, k= {}".format(ar, k) )
    dist_m = scale_distm(dist_m, k)
    n_l, ET_l = [], []
    # for n in range(k, N+1):
    for n in [2, 4, 5, 10, 20]:
      print("n= {}".format(n) )
      n_l.append(n)
      
      ar_ = ar * n/N
      ET = sim_mdsnk(num_frun, ar_, n, k, dist_m)
      print("ET= {}".format(ET) )
      ET_l.append(ET)
    plot.plot(n_l, ET_l, color=next(dark_color), label=r'$\lambda= {}$'.format(ar), marker=next(marker), linestyle=':', mew=2)
    plot.xlabel(r'$n$', fontsize=14)
    plot.title(r'$N= {}$, $k= {}$, $V \sim {}$'.format(N, k, dist_to_latex(dist_m) ) )
  
  print("N= {}, dist_m= {}".format(N, dist_m) )
  # plot_varyingk(ar=2, n=5)
  # plot_varyingk(ar=2.5, n=5)
  # plot_varyingk(ar=3, n=5)
  # plot_varyingk(ar=3.25, n=5)
  # plot_varyingk(ar=3.5, n=5)
  
  # plot_varyingk(ar=0.2, n=5)
  # plot_varyingk(ar=0.8, n=5)
  # plot_varyingk(ar=1.4, n=5)
  
  plot_varyingk(ar=2/10, n=5)
  plot_varyingk(ar=2.5/10, n=5)
  plot_varyingk(ar=2.75/10, n=5)
  plot_varyingk(ar=2.95/10, n=5)
  
  # plot_varyingn(ar=5, k=2)
  
  plot.legend()
  plot.ylabel('Average download time', fontsize=14)
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_N{}_nk.pdf".format(N) )
  log(WARNING, "done; N= {}.".format(N) )

def plot_mds_n2_wrtn():
  N = 20
  k = 2
  dist_m = {'dist': 'Pareto', 'l': 1, 'a': 2}
  
  ar_ub = mds_exactbound_on_ar(N, k, dist_m)
  log(WARNING, "N= {}, k= {}, ar_ub= {}, dist_m={}".format(N, k, ar_ub, dist_m) )
  
  num_frun = 1
  def plot_(ar):
    n_l = []
    ET_sm_l, ET_sim_l = [], []
    
    print("> ar= {}".format(ar) )
    for n in [2, 4, 5, 10, 20]:
    # for n in [1, 2, 5, 10]:
      n_l.append(n)
      
      ar_ = ar*n/N
      ET_sim = sim_mdsnk(num_frun, ar_, n, k, dist_m)
      print("ET_sim= {}".format(ET_sim) )
      if ET_sim is None: break
      ET_sim_l.append(ET_sim)
      
      # ET_sm = ET_mds_nk_sm(ar, n, k, dist_m)
      # print("ET_sm= {}".format(ET_sm) )
      # ET_sm_l.append(ET_sm)
    plot.plot(n_l, ET_sim_l, color=next(dark_color), label=r'$\lambda= {}$'.format(ar), marker=next(marker), linestyle=':', mew=2)
    # plot.plot(n_l, ET_sm_l, color=next(dark_color), label='Split-merge upper bound', marker=next(marker), linestyle=':', mew=2)
    # print("ET_sim_l= {}".format(pprint.pformat(ET_sim_l) ) )
  ar_l = [0.05, 0.25, 0.45, 0.65, 0.85] # [0.05]
  ar_l = [2*ar for ar in ar_l]
  for ar in ar_l:
    plot_(ar)
  plot.legend()
  plot.xlabel(r'$n$', fontsize=14)
  # plot.ylabel(r'$E[T]$', fontsize=14)
  plot.ylabel(r'Average download time', fontsize=14)
  plot.title(r'$N= {}$, $k= {}$, $V \sim {}$'.format(N, k, dist_to_latex(dist_m) ) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.6, def_size[1]/1.3)
  fig.tight_layout()
  plot.savefig("plot_mds_N{}_k{}_{}.pdf".format(N, k, dist_m['dist'] ) )
  log(WARNING, "done; N= {}, k= {}".format(N, k) )

def plot_mds_n2():
  n = 2 # 3 # 10 # 5
  k = 2 # 6 # 4
  dist_m = {'dist': 'Exp', 'mu': 1}
  # dist_m = {'dist': 'SExp', 'D': 1, 'mu': 0.5}
  # dist_m = {'dist': 'Pareto', 'l': 1, 'a': 2}
  
  if dist_m['dist'] == 'Exp':
    ar_ub = 0.9*mds_exactbound_on_ar(n, k, dist_m) # mds_innerbound_on_ar(n, k, dist_m)
  else:
    ar_ub = 0.75*mds_exactbound_on_ar(n, k, dist_m) # 0.7*
  log(WARNING, "n= {}, k= {}, dist_m={}".format(n, k, dist_m) )
  
  ar_l = []
  ET_sm_l, ET_sim_l, ET_approx_l = [], [], []
  ET_approx2_l = []
  ET_varkigauri_lb_l = []
  
  num_frun = 5
  for ar in [*numpy.linspace(0.05, 0.6*ar_ub, 4, endpoint=False), *numpy.linspace(0.6*ar_ub, ar_ub, 7) ]:
    pi_l = []
    ET_sim = sim_mdsnk(num_frun, ar, n, k, dist_m, fi_l=pi_l)
    print("ET_sim= {}".format(ET_sim) )
    if ET_sim is None: break
    ar_l.append(ar)
    ET_sim_l.append(ET_sim)
    
    ET_sm = ET_mds_nk_sm(ar, n, k, dist_m)
    ET_sm = ET_sm if ET_sm is not None and ET_sm < 15 else None
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
  plot.plot(ar_l, ET_sm_l, color=next(dark_color), label='Upper bound of [1]', marker=next(marker), linestyle=':', mew=2)
  # print("ET_sim_l= {}".format(pprint.pformat(ET_sim_l) ) )
  plot.plot(ar_l, ET_sim_l, color=next(dark_color), label='Simulation', marker=next(marker), linestyle=':', mew=2)
  plot.plot(ar_l, ET_approx_l, color=next(dark_color), label=r'$M/G/1$ approximation', marker=next(marker), linestyle=':', mew=2)
  # plot.plot(ar_l, ET_approx2_l, color=next(dark_color), label='Approx_', marker=next(marker), linestyle=':', mew=2)
  if dist_m['dist'] == 'Exp':
    plot.plot(ar_l, ET_varkigauri_lb_l, color=next(dark_color), label='Lower bound of [1]', marker=next(marker), linestyle=':', mew=2)
  plot.legend()
  plot.xlabel(r'$\lambda$', fontsize=14)
  plot.ylabel(r'Average download time', fontsize=14)
  plot.title(r'$n= {}$, $k= {}$, $V \sim {}$'.format(n, k, dist_to_latex(dist_m) ) )
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_mds_{}_{}_{}.pdf".format(n, k, dist_m['dist'] ) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  # plot_mds_n2_wrtn()
  plot_mds_n2()
  # plot_Nnk()
  
  # plot_Nserver_mdsnkjobs()
