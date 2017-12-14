import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import pprint, math, numpy

from mds_sim import *
from commonly_used import *

num_frun = 1
def sim_mds_nrk(ar, n, r, k, dist_m):
  ar = ar/(n/r)
  
  ET_sum = 0
  for f in range(num_frun):
    log(WARNING, "ar= {}, n= {}, r= {}, k= {}, dist_m= {}".format(ar, n, r, k, dist_m) )
    env = simpy.Environment()
    pg = MDS_PG(env, "pg", ar)
    mdsq = MDSQ("mdsq", env, k, range(r), dist_m)
    # monitor = MDSQMonitor(env, mdsq, lambda: 1)
    pg.out = mdsq
    env.run(until=50000 * 2)
    
    st_l = mdsq.jsink.st_l
    ET_sum += float(sum(st_l) )/len(st_l)
    
    total_num_wins = sum([r for i, r in mdsq.jsink.qid__num_win_map.items() ] )
    qid__win_freq_map = {i:float(r)/total_num_wins for i, r in mdsq.jsink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  ET = ET_sum/num_frun
  if ET > 100: return None
  return ET

# ############################################  Model  ########################################### #
def ar_ub(n, r, k, dist_m):
  EV = EXm_n_k(1, r, k, dist_m)
  return float(1/EV)*n/r

def ET_mds_nrk(ar, n, r, k, dist_m):
  ar = ar/(n/r)
  EV = EXm_n_k(1, r, k, dist_m)
  if ar*EV >= 1: return None
  EV2 = EXm_n_k(2, r, k, dist_m)
  ET = EV + ar/2 * EV2/(1 - ar*EV)
  
  if ET < 0 or ET >= 100: return None
  return ET

def plot_n_r_k():
  n = 10
  # dist_m = {'dist': 'Exp', 'mu': 1}
  # dist_m = {'dist': 'SExp', 'D': 0.1, 'mu': 1}
  # dist_m = {'dist': 'Pareto', 'l': 2, 'a': 2}
  dist_m = {'dist': 'Pareto', 'l': 1, 'a': 2}
  log(WARNING, "n= {}, dist_m={}".format(n, dist_m) )
  
  def plot_varying_ar(r, k, sim):
    dist_m_ = scale_dist(dist_m, k)
    ub = ar_ub(n, r, k, dist_m_)
    print("r= {}, k= {}, ar_ub= {}".format(r, k, ub) )
    
    x_l, y_l = [], []
    for ar in numpy.linspace(0.05, ub, 20):
      if sim:
        ET = sim_mds_nrk(ar, n, r, k, dist_m_)
        label = 'Sim, $r= {}$, $k= {}$'.format(r, k)
      else:
        ET = ET_mds_nrk(ar, n, r, k, dist_m_)
        label = 'Split-merge, $r= {}$, $k= {}$'.format(r, k)
      if ET is None: break
      x_l.append(ar)
      y_l.append(ET)
    plot.plot(x_l, y_l, color=next(dark_color), label=label, marker=next(marker), linestyle=':', mew=2)
    plot.xlabel(r'$\lambda$', fontsize=14)
  
  def plot_varying_r(ar, k, sim):
    dist_m_ = scale_dist(dist_m, k)
    print("ar= {}, k= {}".format(ar, k) )
    
    x_l, y_l = [], []
    for r in range(k, n+1):
      if sim:
        ET = sim_mds_nrk(ar, n, r, k, dist_m_)
      else:
        ET = ET_mds_nrk(ar, n, r, k, dist_m_)
      if ET is None: break
      x_l.append(r)
      y_l.append(ET)
    sim_str = 'Sim' if sim else 'Split-merge'
    plot.plot(x_l, y_l, color=next(dark_color), label=r'{}, $\lambda= {0:.2f}$'.format(sim_str, ar), marker=next(marker), linestyle=':', mew=2)
    plot.xlabel(r'$r$', fontsize=14)
  
  def plot_varying_k(ar, r, sim):
    print("ar= {}, r= {}".format(ar, r) )
    
    x_l, y_l = [], []
    for k in range(1, r+1):
      dist_m_ = scale_dist(dist_m, k)
      if sim:
        ET = sim_mds_nrk(ar, n, r, k, dist_m_)
      else:
        ET = ET_mds_nrk(ar, n, r, k, dist_m_)
      if ET is None: break
      x_l.append(k)
      y_l.append(ET)
    sim_str = 'Sim' if sim else 'Split-merge'
    label = r'{}, $\lambda= {}$'.format(sim_str, '{0:.2f}'.format(ar) )
    plot.plot(x_l, y_l, color=next(dark_color), label=label, marker=next(marker), linestyle=':', mew=2)
    plot.xlabel(r'$k$', fontsize=14)
  
  def plot_varying_kr(ar, sim):
    print("ar= {}".format(ar) )
    k_r_grid = numpy.zeros((n, n) )
    k_l, r_l = [], []
    for r in range(1, n+1):
      for k in range(1, r+1):
        dist_m_ = scale_dist(dist_m, k)
        if sim:
          ET = sim_mds_nrk(ar, n, r, k, dist_m_)
        else:
          ET = ET_mds_nrk(ar, n, r, k, dist_m_)
        if ET is None: break
        k_r_grid[k-1, r-1] = ET
    print("k_r_grid= \n{}".format(k_r_grid) )
    extent = [0.5, n+0.5, 0.5, n+0.5]
    img = plot.imshow(k_r_grid, cmap='gray_r', extent=extent, origin='lower')
    plot.colorbar(img, cmap='gray_r')
    
    sim_str = 'Sim' if sim else 'Split-merge'
    plot.title(r'{}, $\lambda= {}$, $n= {}$, $V \sim {}$'.format(sim_str, ar, n, dist_to_latex(dist_m) ) )
    plot.xticks(range(1, n+1), fontsize=12)
    plot.xlabel(r'$r$', fontsize=14)
    plot.yticks(range(1, n+1), fontsize=12)
    plot.ylabel(r'$k$', fontsize=14)
  
  # ub = ar_ub(n, n, 1, dist_m)
  # for ar in numpy.linspace(0.05, ub*0.9, 5):
  #   plot_varying_r(ar, k=1, sim=False)
  # plot.title(r'$n= {}$, $k= 1$, $V \sim {}$'.format(n, dist_to_latex(dist_m) ) )
  
  r = 5
  ub = ar_ub(n, r, 1, dist_m)
  for ar in numpy.linspace(0.2, ub*0.8, 3):
    plot_varying_k(ar, r, sim=True)
    plot_varying_k(ar, r, sim=False)
  plot.title(r'$n= {}$, $r= {}$, $V \sim {}$'.format(n, r, dist_to_latex(dist_m) ) )
  # plot_varying_kr(ar=0.075, sim=False)
  
  # plot_varying_ar(r=1, k=1)
  # plot_varying_ar(r=2, k=1)
  # plot_varying_ar(r=4, k=1)
  # plot_varying_ar(r=5, k=1)
  # plot_varying_ar(r=10, k=1)
  # plot_varying_ar(r=20, k=1)
  
  # plot_varying_ar(r=1, k=1)
  # plot_varying_ar(r=2, k=2)
  # plot_varying_ar(r=4, k=4)
  # plot_varying_ar(r=5, k=5)
  # plot_varying_ar(r=10, k=10)
  # plot_varying_ar(r=20, k=20)
  
  # plot_varying_ar(r=10, k=10)
  # plot_varying_ar(r=10, k=5)
  # plot_varying_ar(r=10, k=2)
  # plot_varying_ar(r=10, k=1)
  
  plot.legend()
  plot.ylabel(r'$E[T]$', fontsize=14)
  # plot.title(r'$n= {}$, $V \sim {}$'.format(n, dist_to_latex(dist_m) ) )
  plot.savefig("plot_n_r_k.pdf")
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_n_r_k()
