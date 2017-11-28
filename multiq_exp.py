import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

from multiq_sim import *

# *********************************************  Sim  ******************************************** #
def sim_multiq(k, num_srun, ar, N, sching_m, k_dist, tsize_dist):
  n = min(N, math.floor(k*sching_m['r'] ) )
  
  E_T_sim, E_C_sim, E_T_tpar, E_C_tpar, E_T_par, E_C_par = 0, 0, 0, 0, 0, 0
  E_jtime_sum = 0
  for s in range(num_srun):
    log(WARNING, "s= {}, ar= {}, N= {}, sching_m= {}".format(s, ar, N, sching_m) )
    
    env = simpy.Environment()
    jg = JG(env, ar, k_dist, tsize_dist)
    if sching_m['t'] == 'red-to-idle':
      mq = MultiQ_RedToIdle(env, N, sching_m)
    else:
      mq = MultiQ(env, N, sching_m)
    jg.out = mq
    jg.init()
    env.run(until=50000*10)
    
    E_jtime = sum(mq.jtime_l)/len(mq.jtime_l)
    # jtime_l = mq.k__jtime_m[k]
    # E_jtime = sum(jtime_l)/len(jtime_l)
    nj_sent_over_departed = jg.nsent/len(mq.jq.deped_jid_l)
    print("k= {}, E_jtime= {}, nj_sent_over_departed= {}".format(k, E_jtime, nj_sent_over_departed) )
    if nj_sent_over_departed > 2:
      return None
    E_jtime_sum += E_jtime
    
    # Plotting the tail of task lifetimes
    # s_l = numpy.sort(mq.tlt_l() )[::-1]
    s_l = numpy.sort(mq.tsl_l() )[::-1]
    i_ = None
    for i in range(len(s_l)-1, 0, -1):
      if s_l[i] > 1.01:
        i_ = i
        break
    s_l = s_l[:i_]
    # y_sim_l = numpy.arange(s_l.size)/s_l.size
    # plot.plot(s_l, y_sim_l, label=r'n-k= {}'.format(sching_m['n-k'] ), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    
    taskt_rv = SimRV(s_l)
    stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(taskt_rv, 0, k, k, n, num_run=10000*10)
    E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
    E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
    print("Sim: E_T= {}, E_C= {}".format(E_T, E_C) )
    E_T_sim += E_T
    E_C_sim += E_C
    
    # l, u, a = fit_tpareto(s_l)
    # # rv = TPareto(l, u, a)
    # # y_l = []
    # # for x in s_l:
    # #   y_l.append(rv.tail(x) )
    # # plot.plot(s_l, y_l, label="TPareto, n-k= {}, a= {}".format(sching_m['n-k'], a), color=next(dark_color), linestyle='-')
    # task_t = "TPareto"
    # task_dist_m = {'l': l, 'u': u, 'a': a}
    # E_T = E_T_k_l_n(task_t, task_dist_m, 0, k, k, n)
    # E_C = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n, w_cancel=True)
    # print("Fitted TPareto: E_T= {}, E_C= {}".format(E_T, E_C) )
    E_T_tpar += 0 # E_T
    E_C_tpar += 0 # E_C
    
    l, a = fit_pareto(s_l)
    # rv = Pareto(l, a)
    # y_l = []
    # for x in s_l:
    #   y_l.append(rv.tail(x) )
    # plot.plot(s_l, y_l, label="Pareto, n-k= {}, a= {}".format(sching_m['n-k'], a), color=next(dark_color), linestyle='-')
    task_t = "Pareto"
    task_dist_m = {'loc': l, 'a': a}
    E_T = E_T_k_l_n(task_t, task_dist_m, 0, k, k, n)
    E_C = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n, w_cancel=True)
    print("Fitted Pareto(l= {}, a= {}):\n E_T= {}, E_C= {}".format(l, a, E_T, E_C) )
    E_T_par += E_T
    E_C_par += E_C
    
    # plot.legend()
    # plot.xscale('log')
    # plot.yscale('log')
    # plot.xlabel(r'Lifetime', fontsize=13)
    # plot.ylabel(r'Tail distribution', fontsize=13)
    # plot.title(r'$T \sim {}$, $\lambda= {}$'.format(tsize_dist, ar) )
    # plot.savefig("plot_psq_n_k_{}.png".format(sching_m['n-k'] ) )
    # plot.gcf().clear()
  E_jtime = E_jtime_sum/num_srun
  E_T_sim, E_C_sim = E_T_sim/num_srun, E_C_sim/num_srun
  E_T_tpar, E_C_tpar = E_T_tpar/num_srun, E_C_tpar/num_srun
  E_T_par, E_C_par = E_T_par/num_srun, E_C_par/num_srun
  
  print(">> E_jtime= {}".format(E_jtime) )
  return {'E_jtime': E_jtime,
          'E_T_sim': E_T_sim, 'E_C_sim': E_C_sim,
          'E_T_tpar': E_T_tpar, 'E_C_tpar': E_C_tpar,
          'E_T_par': E_T_par, 'E_C_par': E_C_par}

def plot_multiq():
  N = 20
  k = 10
  k_dist = DUniform(1, k)
  l, u, a = 1, 10**10, 1.1 # 1.1
  tsize_dist = TPareto(l, u, a)
  tsize_in_latex = r'TPareto(l= {}, u= {}, \alpha= {})'.format(l, u, a)
  # tsize_dist = DUniform(1, 1)
  # tsize_in_latex = r'{}'.format(tsize_dist)
  # ar_ub = 0.7 * N/tsize_dist.mean()/k_dist.mean()
  log(WARNING, "k= {}, N= {}, k_dist= {}, tsize_dist= {}".format(k, N, k_dist, tsize_dist) )
  
  num_srun = 1
  def plot_ESl_vs_ar(sching_t, n_k):
    x_l, y_l = [], []
    sching_m = {'t': sching_t, 'n-k': n_k}
    ar = 0.05
    while True:
      sim_m = sim_multiq(k, num_srun, ar, N, sching_m, k_dist, tsize_dist)
      if sim_m is None: break
      E_jtime = sim_m['E_jtime']
      x_l.append(ar)
      y_l.append(E_jtime)
      ar += 0.1
    
    print("y_l=\n{}".format(pprint.pformat(y_l) ) )
    plot.plot(x_l, y_l, label=r'$n-k= {}$'.format(n_k), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$\lambda$', fontsize=12)
    plot.ylabel(r'$E[Slowdown]$', fontsize=12)
  
  def plot_ESl_vs_r(sching_t, ar):
    x_l, y_l, y_par_l = [], [], []
    # for n_k in range(N - k_dist.u_l):
    # for r in numpy.linspace(1, N/k_dist.u_l, 5):
    # for r in numpy.linspace(1, N/k_dist.l_l, 10):
    E_T_sim_r1 = None
    for r in numpy.arange(1, N/k_dist.l_l, 0.2):
      sching_m = {'t': sching_t, 'r': r}
      sim_m = sim_multiq(k, num_srun, ar, N, sching_m, k_dist, tsize_dist)
      if E_T_sim_r1 is None:
        E_T_sim_r1 = sim_m['E_T_sim']
      if sim_m is None or sim_m['E_T_sim'] > E_T_sim_r1:
        break
      x_l.append(r)
      # y_l.append(sim_m['E_jtime'] )
      y_l.append(sim_m['E_T_sim'] )
      y_par_l.append(sim_m['E_T_par'] )
    plot.plot(x_l, y_l, label='Sim', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(x_l, y_par_l, label='Pareto', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$r$', fontsize=12)
    plot.ylabel(r'$E[Slowdown]$', fontsize=12)
  def plot_EC_vs_ET(ar):
    r_l = []
    x_sim_l, y_sim_l, x_tpar_l, y_tpar_l, x_par_l, y_par_l = [], [], [], [], [], []
    for r in numpy.arange(1, 10, 0.5):
      sching_m = {'t': 'coded', 'r': r}
      sim_m = sim_multiq(k, num_srun, ar, N, sching_m, k_dist, tsize_dist)
      if sim_m is None: break
      r_l.append(r)
      x_sim_l.append(sim_m['E_T_sim'] )
      y_sim_l.append(sim_m['E_C_sim'] )
      x_tpar_l.append(sim_m['E_T_tpar'] )
      y_tpar_l.append(sim_m['E_C_tpar'] )
      x_par_l.append(sim_m['E_T_par'] )
      y_par_l.append(sim_m['E_C_par'] )
    plot.plot(x_sim_l[0], y_sim_l[0], label=r'No redundancy', zorder=2, marker='x', color='blue', mew=3, ms=9)
    plot.plot(x_sim_l, y_sim_l, label=r'Simulation', zorder=1, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(x_tpar_l, y_tpar_l, zorder=0, label=r'Using fitted Truncated-Pareto', color=next(dark_color), linestyle='-.', lw=2)
    plot.plot(x_par_l, y_par_l, zorder=0, label=r'Using fitted Pareto', color=next(dark_color), linestyle='-', lw=2)
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel(r'$E[T]$', fontsize=13)
    plot.ylabel(r'$E[C]$', fontsize=13)
  # plot_ESl_vs_ar('coded', n_k=0)
  # plot_ESl_vs_ar('coded', n_k=2)
  
  # plot_ESl_vs_ar('red-to-idle', n_k=2)
  # plot_ESl_vs_ar('coded', n_k=2)
  
  # plot_ESl_vs_r('red-to-idle', ar=0.2)
  plot_ESl_vs_r('coded', ar=0.4)
  
  # plot_EC_vs_ET(ar=0.4)
  # 
  plot.legend(prop={'size':11} )
  plot.title(r'$N= {}$, $k \sim {}$, $T \sim {}$'.format(N, k_dist, tsize_in_latex) )
  fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_multiq_N_{}.png".format(N) )
  fig.clear()
  log(WARNING, "done; N= {}".format(N) )

if __name__ == "__main__":
  plot_multiq()
  