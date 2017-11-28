import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt, itertools

from mixed_sim import *
from mixed_models import *

def sim_mixednet(num_f_run, n, k, qarrdist_m_l):
  ET_sum, ET2_sum = 0, 0
  EL_sum, EL2_sum = 0, 0
  for f in range(num_f_run):
    log(WARNING, "n= {}, k= {}, qarrdist_m_l=\n {}".format(n, k, pprint.pformat(qarrdist_m_l) ) )
    env = simpy.Environment()
    
    pg = MixedPG(env, "pg", qarrdist_m_l)
    mn = MixedNet(env, n, k)
    # monitor = MNMonitor(env, mn, 0.05)
    pg.out = mn
    env.run(until=50000)
    
    ET_ET2 = mn.ET_ET2()
    ET_sum += ET_ET2[0]
    ET2_sum += ET_ET2[1]
    
    EL_EL2 = (0, 0) # monitor.EL_EL2()
    EL_sum += EL_EL2[0]
    EL2_sum += EL_EL2[1]
    
  return ET_sum/num_f_run, ET2_sum/num_f_run, \
         EL_sum/num_f_run, EL2_sum/num_f_run

def plot_mixednet():
  n = 10
  num_f_run = 1
  
  def plot_ET_vs_ar(k):
    pET, pET2, pEL, pEL2 = True, False, False, False
    
    x_l = []
    ET_sim_l, ET_approx_l, ET_ub_l, ET_lb_l = [], [], [], []
    ET2_sim_l, ET2_approx_l = [], []
    
    EL_sim_l, EL_approx_l = [], []
    EL2_sim_l, EL2_approx_l = [], []
    for ar in numpy.linspace(0.2, 2, 10):
      x_l.append(ar)
      print(">> ar= {}".format(ar) )
      
      dist_m = {'dist': 'Exp', 'mu': ar}
      ET_sim, ET2_sim, EL_sim, EL2_sim = \
        sim_mixednet(num_f_run, n, k, qarrdist_m_l=[dist_m for i in range(n) ] )
      if pET:
        print("ET_sim= {}".format(ET_sim) )
        ET_sim_l.append(ET_sim)
        # ET_approx = ET_n_2(n, ar)
        # ET2_approx = ET2_n_2(n, ar)
        pe = pempty_iteratively(n, k, dist_m)
        ET_approx = ET_mg1_approx(n, k, pe, dist_m)
        print("ET_approx= {}".format(ET_approx) )
        ET_approx_l.append(ET_approx)
      elif pET2:
        print("ET2_sim= {}".format(ET2_sim) )
        ET2_sim_l.append(ET2_sim)
        print("ET2_approx= {}".format(ET2_approx) )
        ET2_approx_l.append(ET2_approx)
      elif pEL:
        print("EL_sim= {}".format(EL_sim) )
        EL_sim_l.append(EL_sim)
        EL_approx = EL_n_2(n)
        print("EL_approx= {}".format(EL_approx) )
        EL_approx_l.append(EL_approx)
      elif pEL2:
        print("EL2_sim= {}".format(EL2_sim) )
        EL2_sim_l.append(EL2_sim)
        EL2_approx = EL2_n_2(n)
        print("EL2_approx= {}".format(EL2_approx) )
        EL2_approx_l.append(EL2_approx)
      # ET_ub_l.append(ET_mixednet_ub(n, k, ar) )
      # ET_lb_l.append(ET_mixednet_lb(n, k, ar) )
    # plot.plot(x_l, ET_lb_l, label=r'Lower-bound', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.plot(x_l, ET_ub_l, label=r'Upper-bound', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    if pET:
      plot.plot(x_l, ET_sim_l, label=r'$k= {}$, Simulation'.format(k), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(x_l, ET_approx_l, label=r'$k= {}$, Approx'.format(k), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'$E[D]$', fontsize=13)
    elif pET2:
      plot.plot(x_l, ET2_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(x_l, ET2_approx_l, label=r'Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'$E[D^2]$', fontsize=13)
    elif pEL:
      plot.plot(x_l, EL_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(x_l, EL_approx_l, label=r'Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'$E[L]$', fontsize=13)
    elif pEL2:
      plot.plot(x_l, EL2_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(x_l, EL2_approx_l, label=r'Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'$E[L^2]$', fontsize=13)
    
    # plot.title(r'$n= {}$, $k= {}$'.format(n, k) )
    plot.title(r'$n= {}$'.format(n) )
    plot.xlabel(r'$\lambda$', fontsize=13)
    
  def plot_ET_vs_k(dist_m):
    pET, pET2, pEL, pEL2 = True, False, False, False
    
    k_l = []
    ET_sim_l, ET_approx_l = [], []
    ET2_sim_l, ET2_approx_l = [], []
    
    EL_sim_l, EL_approx_l = [], []
    EL2_sim_l, EL2_approx_l = [], []
    
    for k in range(2, n):
      k_l.append(k)
      print(">> k= {}".format(k) )
      
      ET_sim, ET2_sim, EL_sim, EL2_sim = \
        sim_mixednet(num_f_run, n, k, qarrdist_m_l=[dist_m for i in range(n) ] )
      
      if pET:
        print("ET_sim= {}".format(ET_sim) )
        ET_sim_l.append(ET_sim)
        pe = pempty_iteratively(n, k, dist_m)
        # pe = approx_pempty_iteratively(n, k, ar)
        ET_approx = ET_mg1_approx(n, k, pe, dist_m)
        print("pe= {}, ET_approx= {}".format(pe, ET_approx) )
        ET_approx_l.append(ET_approx)
        plot.xlabel(r'$E[D]$', fontsize=13)
      elif pET2:
        print("ET2_sim= {}".format(ET2_sim) )
        ET2_sim_l.append(ET2_sim)
        pe = 0 # pempty_iteratively(n, k, dist_m)
        # pe = approx_pempty_iteratively(n, k, dist_m)
        ET2_approx = 0 # ET_mg1_approx(n, k, pe, dist_m)
        print("pe= {}, ET2_approx= {}".format(pe, ET2_approx) )
        ET2_approx_l.append(ET2_approx)
    if pET:
      plot.plot(k_l, ET_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(k_l, ET_approx_l, label=r'M/G/1 Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'$E[D]$', fontsize=13)
    elif pET2:
      plot.plot(k_l, ET2_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(k_l, ET2_approx_l, label=r'M/G/1 Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'$E[D^2]$', fontsize=13)
    dist = dist_m['dist']
    if dist == 'Exp':
      X_latex = r'Exp(\mu= {})'.format(dist_m['mu'] )
    elif dist == 'Pareto':
      X_latex = r'Pareto(l= {}, \alpha= {})'.format(dist_m['loc'], dist_m['a'] )
    plot.title(r'$n= {}$, $X \sim {}$'.format(n, X_latex) )
    plot.xlabel(r'$k$', fontsize=12)
  # 
  # plot_ET_vs_ar(k=9)
  # plot_ET_vs_ar(k=7)
  # plot_ET_vs_ar(k=3)
  
  # dist_m = {'dist': 'Exp', 'mu': ar}
  dist_m = {'dist': 'Pareto', 'loc': 1, 'a': 10}
  plot_ET_vs_k(dist_m)
  #
  plot.legend()
  # plot.ylabel(r'Throughput', fontsize=13)
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mixednet_n_{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

# ****************************** Steady state prob  ******************************* #
def sim_steadystate_prob_mixednet(n, k, qlambda_l):
  log(WARNING, "n= {}, k= {}, qlambda_l= {}".format(n, k, qlambda_l) )
  env = simpy.Environment()
  
  pg = MixedPG(env, _id="pg", qlambda_l=qlambda_l)
  mn = MixedNet(env, n, k)
  monitor = MNMonitor(env, mn, 0.05)
  pg.out = mn
  env.run(until=50000)
  
  # for qid, state_counter_map in monitor.qid__state_counter_map_map.items():
  #   total_c = sum([c for s,c in state_counter_map.items() ] )
  #   monitor.qid__state_counter_map_map[qid] = {s:c/total_c for s,c in state_counter_map.items() }
  # print("qid__state_counter_map_map=\n {}".format(pprint.pformat(monitor.qid__state_counter_map_map) ) )
  
  sstate_prob_map = monitor.steadystate_prob_map()
  # log(WARNING, "sstate_prob_map= {}".format(pprint.pformat(sstate_prob_map) ) )
  return sstate_prob_map

def plot_mixednet_steadystate():
  n, k = 100, 3
  
  k_l, Pr_busy_sim_l, Pr_busy_l = [], [], []
  
  # for k in range(1, n+1, 2):
  for k in range(1, n+1, 20):
    k_l.append(k)
    for l in numpy.linspace(1, 1, 1):
      qlambda_l=n*[l]
      steadystate_prob_map = sim_steadystate_prob_mixednet(n, k, qlambda_l)
      print("k= {}, l= {}, steadystate_prob_map= \n{}".format(k, l, pprint.pformat(steadystate_prob_map) ) )
      Pr_busy_sim = 1 - steadystate_prob_map[0]
      print("Pr_busy_sim= {}, Pr_busy= {}".format(Pr_busy_sim, Pr_busy(n, k) ) )
      Pr_busy_sim_l.append(Pr_busy_sim)
      Pr_busy_l.append(Pr_busy(n, k) )
  mew, ms = 3, 5
  plot.plot(k_l, Pr_busy_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  plot.plot(k_l, Pr_busy_l, label=r'Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.legend()
  plot.xlabel(r'$k$', fontsize=13)
  plot.ylabel(r'Steady-state probability of busy', fontsize=13)
  fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mixednet_steadystate_n_{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_mixednet()
  # plot_mixednet_steadystate()
  # Pr_busy_mixednet_approx()
  
  