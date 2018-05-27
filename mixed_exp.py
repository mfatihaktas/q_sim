import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt, itertools

from mixed_sim import *
from mixed_models import *
from mixed_newmodels import *

def plot_delay_dist():
  n = 10
  
  def sim_mu__approx_mu(ar, k):
    dist_m = {'dist': 'Exp', 'mu': ar}
    
    env = simpy.Environment()
    pg = MixedPG(env, "pg", [dist_m for i in range(n) ] )
    mn = MixedNet(env, n, k)
    pg.out = mn
    env.run(until=50000*2)
    
    ET_ET2 = mn.ET_ET2()
    print("ET= {}, ET2= {}".format(ET_ET2[0], ET_ET2[1] ) )
    
    qt_l = mn.qt_l()
    print("len(mn.qt_l)= {}".format(len(qt_l) ) )
    s_l = numpy.sort(qt_l)
    x_l = s_l[::-1]
    i_ = None
    for i in range(len(x_l)-1, 0, -1):
      if x_l[i] > 0.0001: i_ = i; break
      # if x_l[i] > 1: i_ = i; break
    x_l = x_l[:i_]
    
    y_l = numpy.arange(x_l.size)/x_l.size
    # plot.plot(x_l, y_l, label="k= {}".format(k), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    
    x__l, logy_l = [], []
    for i, y in enumerate(y_l):
      if y != 0:
        x__l.append(x_l[i] )
        logy_l.append(math.log(y) )
    [a, b] = numpy.polyfit(x__l, logy_l, 1)
    print("a= {}, b= {}".format(a, b) )
    # y_l = []
    # for x in x_l:
    #   y_l.append(math.exp(a*x + b) )
    # plot.plot(x_l, y_l, label="fitted, k= {}".format(k), marker=next(marker), color=next(dark_color), linestyle='-')
    # plot.xlabel(r'$d$')
    # plot.ylabel(r'$Pr\{D > d\}$')
    # plot.yscale('log')
    mu = tail_exponent(n, k, dist_m)
    print("mu= {}".format(mu) )
    return a, mu
  
  def plot_mu_vs_k(ar):
    print(">> ar= {}".format(ar) )
    k_l, sim_mu_l, mu_l = [], [], []
    for k in range(2, n):
      k_l.append(k)
      sim_mu, mu = sim_mu__approx_mu(ar, k)
      sim_mu_l.append(sim_mu)
      mu_l.append(mu)
    plot.plot(k_l, sim_mu_l, label="Simulation", marker=next(marker), color=next(dark_color), linestyle=':')
    plot.plot(k_l, mu_l, label="Approximation", marker=next(marker), color=next(dark_color), linestyle=':')
    plot.xlabel(r'$k$', fontsize=13)
    plot.ylabel(r'$\mu$', fontsize=14)
    plot.title(r'$n= {}$, $X \sim Exp(\lambda= {})$'.format(n, ar) )
  
  def plot_mu_vs_ar():
    def plot_(k):
      print(">> k= {}".format(k) )
      ar_l, sim_mu_l, mu_l = [], [], []
      for ar in numpy.linspace(0.2, 2, 10):
        ar_l.append(ar)
        sim_mu, mu = sim_mu__approx_mu(ar, k)
        sim_mu_l.append(sim_mu)
        mu_l.append(mu)
      plot.plot(ar_l, sim_mu_l, label=r'$k= {}$, Simulation'.format(k) , marker=next(marker), color=next(dark_color), linestyle=':')
      plot.plot(ar_l, mu_l, label=r'$k= {}$, Approx'.format(k), marker=next(marker), color=next(dark_color), linestyle=':')
    plot_(k=3)
    plot_(k=7)
    plot_(k=9)
    
    plot.xlabel(r'$\lambda$', fontsize=13)
    plot.ylabel(r'$\mu$', fontsize=14)
    plot.title(r'$n= {}$, $X \sim Exp(\lambda)$'.format(n) )
  
  plot_mu_vs_k(ar=1)
  # plot_mu_vs_ar()
  
  fig = plot.gcf()
  fig.tight_layout()
  plot.legend()
  plot.savefig("plot_delay_dist_n_{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

def sim_mixednet(num_frun, n, k, qarrdist_m_l):
  ET_sum, ET2_sum = 0, 0
  EL_sum, EL2_sum = 0, 0
  E_Deanont_sum = 0
  for f in range(num_frun):
    # log(WARNING, "n= {}, k= {}, qarrdist_m_l=\n {}".format(n, k, pprint.pformat(qarrdist_m_l) ) )
    log(WARNING, "n= {}, k= {}".format(n, k) )
    env = simpy.Environment()
    
    pg = MixedPG(env, "pg", qarrdist_m_l)
    attacker = None # AttackOne(env, n, k) # StateSniffer(env, n, k)
    mn = MixedNet(env, n, k, attacker)
    # monitor = MNMonitor(env, mn, 0.05)
    pg.out = mn
    env.run(until=50000*1)
    
    ET_ET2 = mn.ET_ET2()
    ET_sum += ET_ET2[0]
    ET2_sum += ET_ET2[1]
    
    EL_EL2 = (0, 0) # monitor.EL_EL2()
    EL_sum += EL_EL2[0]
    EL2_sum += EL_EL2[1]
    
    E_Deanont = sum(mn.attackt_l)/len(mn.attackt_l) if len(mn.attackt_l) else 0
    print("num_deanon= {}, E_Deanont= {}".format(len(mn.attackt_l), E_Deanont) )
    # E_EstedDeanont = sum(mn.ested_attackt_l)/len(mn.ested_attackt_l)
    # print("E_EstedDeanont= {}".format(E_EstedDeanont) )
    
    E_Deanont_sum += E_Deanont
  
  return ET_sum/num_frun, ET2_sum/num_frun, \
         EL_sum/num_frun, EL2_sum/num_frun, \
         E_Deanont_sum/num_frun

def plot_mixednet():
  num_frun = 1
  
  def plot_ET_vs_ar(n, k):
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
      # ET_sim, ET2_sim, EL_sim, EL2_sim = \
      #  sim_mixednet(num_frun, n, k, qarrdist_m_l=[dist_m for i in range(n) ] )
      ET_sim, ET2_sim, EL_sim, EL2_sim, EDeanont_sim = \
        sim_mixednet(num_frun, n, k, qarrdist_m_l=[dist_m for i in range(n) ] )
      if pET:
        print("ET_sim= {}".format(ET_sim) )
        ET_sim_l.append(ET_sim)
        # ET_approx = ET_n_2(n, ar)
        # ET2_approx = ET2_n_2(n, ar)
        
        ET_mg1_ = ET_mg1approx_(n, k, dist_m)
        print("ET_mg1_= {}".format(ET_mg1_) )
        
        ET_mg1efs = ET_mg1efsapprox(n, k, dist_m)
        print("ET_mg1efs= {}".format(ET_mg1efs) )
        
        ET_mg1efs_ = ET_mg1efsapprox_(n, k, dist_m)
        print("ET_mg1efs_= {}".format(ET_mg1efs_) )
        ET_approx_l.append(ET_mg1efs_)
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
      plot.plot(x_l, ET_approx_l, label=r'$k= {}$, M/G/1/efs approximation'.format(k), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'Average delay', fontsize=14)
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
    plot.title(r'Batch mix, $n= {}$'.format(n) )
    plot.xlabel(r'$\lambda$', fontsize=14)
    
  def plot_ET_vs_k(n, dist_m):
    pET, pET2, pEL, pEL2 = True, False, False, False
    
    k_l = []
    ET_sim_l, ET_approx_l = [], []
    ET2_sim_l, ET2_approx_l = [], []
    
    EL_sim_l, EL_approx_l = [], []
    EL2_sim_l, EL2_approx_l = [], []
    
    for k in range(2, n, 1):
    # for k in numpy.linspace(40, n, 2, endpoint=False):
      k = int(k)
      k_l.append(k)
      print(">> k= {}".format(k) )
      
      ET_sim, ET2_sim, EL_sim, EL2_sim, EDeanont_sim = \
        sim_mixednet(num_frun, n, k, qarrdist_m_l=[dist_m for i in range(n) ] )
      # ET_sim, ET2_sim, EL_sim, EL2_sim, EDeanont_sim = None, None, None, None, None
      
      if pET:
        print("ET_sim= {}".format(ET_sim) )
        ET_sim_l.append(ET_sim)
        
        ET_mg1_ = ET_mg1approx_(n, k, dist_m)
        print("ET_mg1_= {}".format(ET_mg1_) )
        
        ET_mg1efs = ET_mg1efsapprox(n, k, dist_m)
        print("ET_mg1efs= {}".format(ET_mg1efs) )
        
        ET_mg1efs_ = ET_mg1efsapprox_(n, k, dist_m)
        print("ET_mg1efs_= {}".format(ET_mg1efs_) )
        ET_approx_l.append(ET_mg1efs_)
        
        # ETtighterlb = ET_tighterlb(n, k, dist_m)
        # print("ETtighterlb= {}".format(ETtighterlb) )
      elif pET2:
        print("ET2_sim= {}".format(ET2_sim) )
        ET2_sim_l.append(ET2_sim)
        ET2_approx = 0 # ET_mg1approx(n, k, dist_m)
        print("ET2_approx= {}".format(ET2_approx) )
        ET2_approx_l.append(ET2_approx)
    if pET:
      plot.plot(k_l, ET_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(k_l, ET_approx_l, label=r'M/G/1/efs approximation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'Average delay', fontsize=13)
    elif pET2:
      plot.plot(k_l, ET2_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.plot(k_l, ET2_approx_l, label=r'M/G/1 approximation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
      plot.ylabel(r'$E[D^2]$', fontsize=13)
    
    # plot.title(r'$n= {}$, $X \sim {}$'.format(n, dist_to_latex(dist_m) ) )
    plot.title(r'Batch mix, $n= {}$, $\lambda= {}$'.format(n, dist_m['mu'] ), fontsize=13)
    plot.xlabel(r'$k$', fontsize=14)
  
  def plot_ET_vs_EDeanont(dist_m):
    k_l = []
    ET_sim_l, EDeanont_sim_l = [], []
    ET_over_EDeanont_sim_l = []
    
    for k in range(2, n):
      k_l.append(k)
      print(">> k= {}".format(k) )
      
      ET_sim, ET2_sim, EL_sim, EL2_sim, EDeanont_sim = \
        sim_mixednet(num_frun, n, k, qarrdist_m_l=[dist_m for i in range(n) ] )
      
      print("ET_sim= {}".format(ET_sim) )
      ET_sim_l.append(ET_sim)
      EDeanont_sim_l.append(EDeanont_sim)
      
      ET_over_EDeanont_sim_l.append(ET_sim/EDeanont_sim)
    plot.plot(EDeanont_sim_l, ET_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.xlabel(r'$E[Deanont]$', fontsize=13)
    plot.ylabel(r'$E[D]$', fontsize=13)
    
    # plot.plot(k_l, ET_over_EDeanont_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.xlabel(r'$k$', fontsize=13)
    # plot.ylabel(r'$E[D]/E[Deanont]$', fontsize=13)
    # plot.title(r'$n= {}$, $X \sim {}$'.format(n, dist_to_latex(dist_m) ) )
  
  def plot_EDeanont_vs_k(dist_m):
    k_l = []
    EDeanont_sim_l = []
    
    for k in range(2, n):
      k_l.append(k)
      print(">> k= {}".format(k) )
      
      ET_sim, ET2_sim, EL_sim, EL2_sim, EDeanont_sim = \
        sim_mixednet(num_frun, n, k, qarrdist_m_l=[dist_m for i in range(n) ] )
      
      print("EDeanont_sim= {}".format(EDeanont_sim) )
      EDeanont_sim_l.append(EDeanont_sim)
    plot.plot(k_l, EDeanont_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.xlabel(r'$k$', fontsize=13)
    plot.ylabel(r'$E[Deanont]$', fontsize=13)
    plot.title(r'$n= {}$, $X \sim {}$'.format(n, dist_to_latex(dist_m) ) )
  #
  # plot_ET_vs_ar(k=9)
  # plot_ET_vs_ar(k=7)
  # plot_ET_vs_ar(k=3)
  
  n = 40
  dist_m = {'dist': 'Exp', 'mu': 1}
  print("n= {}, X ~ {}".format(n, dist_m) )
  # dist_m = {'dist': 'Pareto', 'loc': 1, 'a': 50}
  plot_ET_vs_k(n, dist_m)
  # plot_ET_vs_ar(n, k=10)
  # plot_ET_vs_ar(n, k=20)
  
  # plot_ET_vs_EDeanont(dist_m)
  # plot_EDeanont_vs_k(dist_m)
  #
  plot.legend()
  # plot.ylabel(r'Throughput', fontsize=13)
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mixednet_n{}.pdf".format(n) )
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
  # plot_delay_dist()
  
  # plot_mixednet_steadystate()
  # Pr_busy_mixednet_approx()
  
  
