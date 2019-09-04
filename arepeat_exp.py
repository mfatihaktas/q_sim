import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams['ps.useafm'] = True
# matplotlib.rcParams['pdf.use14corefonts'] = True
# matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import matplotlib.pyplot as plot

from arepeat_models import *
from arepeat_sim import *

def plot_k_nc():
  K = 10
  N = 14
  D, mu = 30, 1
  loc, a = 3, 2
  task_t = "SExp" # "Pareto" # "Exp"
  
  if task_t == "Exp": task_t_in_latex = r'{}(\mu={})'.format(r'\mathrm{Exp}', mu)
  elif task_t == "SExp": task_t_in_latex = r'{}(s={}, \mu={})'.format(r'\mathrm{SExp}', D/K, mu)
  elif task_t == "Pareto": task_t_in_latex = r'{}(s={}, \alpha={})'.format(r'\mathrm{Pareto}', loc, a)
  
  ET = True
  num_srun = 10000
  def plot_(k, n=None, c=None, sim=False):
    l = k
    if task_t == "Exp":
      task_t_rv = Exp(mu)
      task_dist_m = {'mu': mu}
    elif task_t == "SExp":
      task_t_rv = Exp(mu, D/k)
      task_dist_m = {'D': D, 'mu': mu}
    elif task_t == "Pareto":
      task_t_rv = Pareto(loc, a)
      task_dist_m = {'loc': loc, 'a': a}
    elif task_t == "Google": task_t_rv = Google(k)
    
    d_l = []
    y_l, y_approx_l, y_sim_l = [], [], []
    u_l = 3*H(K)/mu + 1
    for d in numpy.arange(0, u_l, 0.5):
      d_l.append(d)
      if task_t != "Google":
        if n is not None:
          y_l.append(
            E_T_k_l_n(task_t, task_dist_m, d, k, l, n) if ET else \
            E_C_k_l_n(task_t, task_dist_m, d, k, l, n, w_cancel=True) )
        elif c is not None:
          y_l.append(
            E_T_k_c(task_t, task_dist_m, d, k, c) if ET else \
            E_C_k_l_n(task_t, task_dist_m, d, k, c, w_cancel=True) )
      if sim:
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=num_srun)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=num_srun)
        y_sim_l.append(
          sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] ) if ET else \
          sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] ) )
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    plot.plot(d_l, y_l, label=label, color=next(dark_color), marker=next(marker), ms=10, mew=0.1, ls=':')
    if sim:
      plot.plot(d_l, y_sim_l, label='Simulation, {}'.format(label), color=next(light_color), marker=next(marker), ms=10, mew=0.1, ls='')
  
  sim = True # False
  # k_step = 2
  # for n in [K, K+1, *range(K+k_step, N, k_step)]:
  #   plot_(K, n=n, sim=sim)
  # plot_(K, n=20, sim=sim)
  
  K = 100
  # plot_(K, n=int(K*1.1), sim=sim)
  # plot_(K, n=int(K*1.4), sim=sim)
  # plot_(K, n=2*K, sim=sim)
  
  plot_(K, c=1, sim=sim)
  plot_(K, c=2, sim=sim)
  
  plot.legend(loc='lower right', fontsize=14, framealpha=0.25, numpoints=1)
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.title(r'${} \sim {}, k= {}$'.format(r'\mathrm{Tasks}', task_t_in_latex, K), fontsize=18)
  plot.xlabel(r'$\Delta$', fontsize=18)
  plot.ylabel('Latency' if ET else 'Cost', fontsize=18)
  prettify(plot.gca() )
  plot.savefig("plot_k_nc_{}_k_{}.png".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )
  
def plot_k_nc_E_C_vs_E_T():
  K = 400 # 1050 # 1000 # 400 # 10
  D, mu = 30, 1
  loc, a = 3, 2
  w_cancel = True
  task_t = "Google" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim SExp(D/k, \mu), D={}, \mu={}'.format(D, mu)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  elif task_t == "Google": task_t_in_latex = r'X \sim Google'
  
  def plot_(k, n=None, c=None, sim=False):
    l = k
    E_T_l, E_T_sim_l, E_T_approx_l = [], [], []
    E_C_l, E_C_sim_l = [], []
    
    u_l = 10*task_t_rv.mean()
    for d in numpy.linspace(0, u_l, 20):
      if task_t != "Google":
        E_T_l.append(E_T_k_l_n(task_t, task_dist_m, d, k, l, n) )
        E_C_l.append(E_C_k_l_n(task_t, task_dist_m, d, k, l, n, w_cancel=w_cancel) )
      if sim:
        if task_t == "Exp":
          task_t_rv = Exp(mu)
          task_dist_m = {'mu': mu}
        elif task_t == "SExp":
          task_t_rv = Exp(mu, D/k)
          task_dist_m = {'D': D, 'mu': mu}
        elif task_t == "Pareto":
          task_t_rv = Pareto(loc, a)
          task_dist_m = {'loc': loc, 'a': a}
        elif task_t == "Google": task_t_rv = Google(k)
        
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=num_srun)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=num_srun)
        E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
        E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )/60/60
        if task_t == "Google":
          E_T = E_T/60/60
          E_C = E_C/60/60
        E_T_sim_l.append(E_T)
        E_C_sim_l.append(E_C)
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    # plot.plot(E_T_l, E_C_l, label=label, color=next(dark_color), ls=':', mew=0.1)
    plot.plot(E_T_sim_l, E_C_sim_l, label='Sim, {}'.format(label), marker=next(marker), color=next(dark_color), ls=':', mew=0.1)
  sim = True
  plot_(K, n=math.floor(K*1.1), sim=sim)
  plot_(K, c=1, sim=sim)
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.8, 0.6) )
  unit = '(s)' if task_t != "Google" else '(hours)'
  plot.xlabel(r'Expected latency $E[T]$ {}'.format(unit) )
  plot.ylabel(r'Expected cost $E[C^c]$ {}'.format(unit) )
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  fig = plot.gcf()
  fig.tight_layout()
  plot.savefig("plot_k_nc_E_C_vs_E_T_{}_k_{}.png".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

# ***************************************  REP vs. CODING  *************************************** #
def plot_EC_vs_ET_wdelay(w_cancel=True):
  num_run = 100000 # 100
  K = 10
  D = 30
  mu = 0.5
  loc, a = 3, 2
  task_t = "Pareto" # "SExp" # "Exp"
  task_t_rv, task_t_in_latex = None, None
  if task_t == "Exp":
    task_t_in_latex = r'{} \sim {}(\mu={})'.format(r'\mathrm{Tasks}', r'\mathrm{Exp}', mu)
    task_t_rv = Exp(mu)
    task_dist_m = {'mu': mu}
  elif task_t == "SExp":
    task_t_in_latex = r'k= {}, {} \sim {}(s={}, \mu={})'.format(K, r'\mathrm{Tasks}', r'\mathrm{SExp}', int(D/K), mu)
    task_t_rv = Exp(mu, D/K)
    task_dist_m = {'D': D, 'mu': mu}
  elif task_t == "Pareto":
    task_t_in_latex = r'k= {}, {} \sim {}(s={}, \alpha={})'.format(K, r'\mathrm{Tasks}', r'\mathrm{Pareto}', loc, a)
    task_t_rv = Pareto(loc, a)
    task_dist_m = {'loc': loc, 'a': a}
  elif task_t == "Google":
    task_t_in_latex = r'k= {}, {} \sim {}'.format(K, r'\mathrm{Tasks}', r'\mathrm{Google}')
    task_t_rv = Google(K)
  
  ann_color = 'black'
  def plot_d_extremes(c=0, n=0):
    l, k, d = K, K, 0
    E_T_l, E_C_l = [], []
    if c:
      for c_ in range(1, c+1):
        if task_t == "SExp":
          E_T = E_T_k_c(task_t, task_dist_m, d, k, c_)
          E_C = E_C_k_c(task_t, task_dist_m, d, k, c_, w_cancel=w_cancel)
        else:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run)
          E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
          E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        if c_ == c:
          plot.annotate(r'$\Delta=0$', ha='center', va='center', xy=(E_T, E_C), xytext=(E_T, E_C-5), color=ann_color, fontsize=20)
          break
      plot.plot(E_T_l, E_C_l, color=ann_color, alpha=0.8, ls='--')
    if n:
      E_T_l.clear()
      E_C_l.clear()
      for n_ in range(k+1, n+1):
        if task_t == "SExp":
          E_T = E_T_k_l_n(task_t, task_dist_m, d, k, l, n_)
          E_C = E_C_k_l_n(task_t, task_dist_m, d, k, l, n_, w_cancel=w_cancel)
        else:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n, num_run)
          E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
          E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        if n_ == n:
          plot.annotate(r'$\Delta=0$', ha='center', va='center', xy=(E_T, E_C), xytext=(E_T+0.1, E_C-5), color=ann_color, fontsize=20)
          # plot.annotate(r'$\Delta=0$', ha='center', va='center', xy=(E_T, E_C), xytext=(E_T+0.75, E_C-8), color=ann_color, fontsize=18)
      plot.plot(E_T_l, E_C_l, color=ann_color, alpha=0.8, ls='--')
  
  def plot_(k=K, n=0, c=0, sim=False):
    l = k
    E_T_l, E_T_sim_l, E_C_l, E_C_sim_l = [], [], [], []
    u_l = 100 # 40 # 25 + 1
    color = next(dark_color)
    if c:
      for d in numpy.arange(0, u_l, 0.5):
        if sim:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run)
          E_T_sim_l.append(sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] ) )
          key = 'C_wc' if w_cancel else 'C_wc'
          E_C_sim_l.append(sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] ) )
        else:
          E_T = E_T_k_c(task_t, task_dist_m, d, k, c)
          E_C = E_C_k_c(task_t, task_dist_m, d, k, c, w_cancel=w_cancel)
          E_T_l.append(E_T)
          E_C_l.append(E_C)
      if sim:
        plot.plot(E_T_sim_l, E_C_sim_l, label=r'Rep,$c={}$'.format(c), color=next(dark_color), marker=next(marker), zorder=2, ms=14, mew=0.1, lw=2, ls='-')
      else:
        plot.plot(E_T_l, E_C_l, label=r'Rep,$c={}$'.format(c), color=next(dark_color), marker=next(marker), zorder=2, ms=14, mew=0.1, lw=2, ls='-')
    elif n:
      for d in numpy.arange(0, u_l, 0.5):
        if sim:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n, num_run)
          E_T_sim_l.append(sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] ) )
          key = 'C_wc' if w_cancel else 'C_wc'
          E_C_sim_l.append(sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] ) )
        else:
          E_T = E_T_k_l_n(task_t, task_dist_m, d, k, l, n)
          E_C = E_C_k_l_n(task_t, task_dist_m, d, k, l, n, w_cancel=w_cancel)
          E_T_l.append(E_T)
          E_C_l.append(E_C)
          # print("d= {}, E_T= {}, E_C= {}".format(d, E_T, E_C) )
      if sim:
        plot.plot(E_T_sim_l, E_C_sim_l, label=r'Code,$n={}$'.format(n), color=next(dark_color), marker=next(marker), zorder=1, ms=14, mew=0.1, ls=':')
      else:
        plot.plot(E_T_l, E_C_l, label=r'Code,$n={}$'.format(n), color=next(dark_color), marker=next(marker), zorder=1, ms=14, mew=0.1, ls=':')
  sim = True # False
  plot_(c=1, sim=sim)
  plot_(c=2, sim=sim)
  # plot_d_extremes(c=2)
  
  plot_(n=K+1, sim=sim)
  plot_(n=K+2, sim=sim)
  for n in range(K+5, 3*K+1, 10):
    plot_(n=n, sim=sim)
  plot_(n=3*K, sim=sim)
  plot_d_extremes(n=3*K)
  
  #
  d, c = 0, 0
  if task_t == "SExp":
    x_nored = E_T_k_c(task_t, task_dist_m, d, K, c)
    y_nored = E_C_k_c(task_t, task_dist_m, d, K, c, w_cancel)
  else:
    stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, K, c, num_run)
    x_nored = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
    y_nored = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
  plot.plot([x_nored], [y_nored], 'o', zorder=3, mew=4, color=ann_color)
  # plot.annotate(r'$\Delta \to \infty$', xy=(x_nored, y_nored), ha='center', va='center', xytext=(x_nored+0.45, y_nored), color=ann_color, fontsize=20)
  if task_t == 'SExp':
    plot.annotate(r'$\Delta \to \infty$', xy=(x_nored, y_nored), ha='center', va='center', xytext=(x_nored+0.5, y_nored), color=ann_color, fontsize=20)
    plot.xlim(([3.5, 9.5] ) )
    plot.ylim(([45, 115] ) )
  elif task_t == 'Pareto':
    plot.annotate(r'$\Delta \to \infty$', xy=(x_nored, y_nored), ha='left', va='center', xytext=(x_nored+0.5, y_nored), color=ann_color, fontsize=20)
  
  # legend = plot.legend()
  legend = plot.legend(loc='center left', bbox_to_anchor=(0.8, 0.7), fontsize=18, framealpha=0.25)
  # frame = legend.get_frame()
  # frame.frameon=True
  # frame.set_edgecolor('black')
  plot.title(r'${}$'.format(task_t_in_latex), fontsize=24)
  plot.xlabel('Latency', fontsize=24)
  plot.ylabel('Cost', fontsize=24)
  prettify(plot.gca() )
  # fig = plot.gcf()
  # plot.savefig("plot_EC_vs_ET_wdelay_{}_k_{}.png".format(task_t, K, w_cancel), bbox_extra_artists=(legend,), bbox_inches='tight', dpi=fig.dpi)
  plot.savefig("plot_EC_vs_ET_wdelay_{}_k_{}.png".format(task_t, K, w_cancel), bbox_inches='tight', dpi=fig.dpi)
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_cost_reduction_vs_redundant_occupation():
  K = 15
  D, mu = 30, 0.5
  loc, a = 3, 1.5 # 3
  task_t = "Pareto" # "SExp"
  if task_t == "SExp":
    task_dist_m = {'D': D, 'mu': mu}
    task_t_in_latex = r'X \sim SExp(D/k, \mu={}), D={}'.format(mu, D)
  elif task_t == "Pareto":
    task_dist_m = {'loc': loc, 'a': a}
    task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  
  def E_redundant_server_occupation(n, k):
    E_O = 0
    # E_O = E_X_n_k_pareto(loc, a, n, n-k)
    if task_t == "Pareto":
      for i in range(1, n-k+1):
        E_O += E_X_n_k_pareto(loc, a, n, i)
    elif task_t == "SExp":
      E_O = D + (H(n) - H(n-k) )/mu
    return E_O
  
  x_l, y_l = [], []
  def plot_(k, n):
    E_C_base = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n=k, w_cancel=True)
    for n_ in numpy.arange(k, n+1, 1):
      E_C = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n_, w_cancel=True)
      E_O = E_redundant_server_occupation(n_, k)
      
      x_l.append(E_O)
      y_l.append(E_C_base-E_C)
    plot.plot(x_l, y_l, label='k= {}'.format(k), color=next(dark_color), marker=next(marker), zorder=0, mew=0.1, ls=':')
    plot.xlabel(r'Redundant server occupation', fontsize=12)
    plot.ylabel(r'Cost reduction', fontsize=12)
  
  def plot_E_X_n_k_for_incing_k(k, n):
    for k_ in numpy.arange(k, n+1, 1):
      x_l.append(k_)
      y_l.append(E_X_n_k_pareto(loc, a, n, k_) )
    plot.plot(x_l, y_l, label='n= {}'.format(n), color=next(dark_color), marker=next(marker), zorder=0, mew=0.1, ls=':')
    plot.xlabel(r'$k$', fontsize=12)
    plot.ylabel(r'$E[X_{n:k}]$', fontsize=12)
  # plot_(k=K, n=2*K)
  plot_E_X_n_k_for_incing_k(k=K, n=2*K)
  
  plot.legend()
  plot.title(r'${}$'.format(task_t_in_latex) )
  
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.1, def_size[1]/1.1)
  plot.savefig("plot_cost_reduction_vs_redundant_occupation_{}.png".format(task_t) )
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_zerodelay_reped_vs_coded(loc, a):
  w_cancel = True
  K = 10 # 10 # 1050 # 15 # 400 # 10
  D, mu = 30, 0.5
  loc, a = 3, 4 # 1.2 # 1.5
  load_m = {'ro_0': 0.35, 'a_0': a} # None
  # l, u, a = 1, 100, 1.5
  task_t = "Pareto" # "SExp" # "TPareto" # "Google" # "Exp"
  
  if task_t == "Exp":
    task_t_in_latex = r'{} \sim {}(\mu={})'.format(r'\mathrm{Tasks}', r'\mathrm{Exp}', mu)
    task_t_rv = Exp(mu)
    task_dist_m = {'mu': mu}
  elif task_t == "SExp":
    task_t_in_latex = r'{} \sim {}(s={}, \mu={})'.format(r'\mathrm{Tasks}', r'\mathrm{SExp}', int(D/K), mu)
    task_t_rv = Exp(mu, D/K)
    task_dist_m = {'D': D, 'mu': mu}
  elif task_t == "Pareto":
    task_t_in_latex = r'{} \sim {}(s={}, \alpha={})'.format(r'\mathrm{Tasks}', r'\mathrm{Pareto}', loc, a)
    task_t_rv = Pareto(loc, a)
    task_dist_m = {'loc': loc, 'a': a}
  elif task_t == "TPareto":
    task_t_in_latex = r'{} \sim {}(l={}, u={}, \alpha={})'.format(r'\mathrm{Tasks}', r'\mathrm{Pareto}', l, u, a)
    task_t_rv = TPareto(l, u, a)
    task_dist_m = {"l": l, "u": u, 'a': a}
  elif task_t == "Google":
    task_t_in_latex = r'{} \sim {}'.format(r'\mathrm{Tasks}', r'\mathrm{Google}')
    task_t_rv = Google(K)
  d = 0
  mew = 3
  num_run = 10000*10 # * 4
  first_moment, second_moment, stdev, coeffvar = True, False, False, False
  def plot_(k=K, n=0, c=0, sim=False):
    x_l, x_sim_l, y_l, y_sim_l = [], [], [], []
    xerr_l, yerr_l = [], []
    color = next(dark_color)
    if c:
      for c_ in range(c+1):
        if task_t != "Google":
          E_T = E_T_k_c(task_t, task_dist_m, d, k, c_, load_m)
          E_C = E_C_k_c(task_t, task_dist_m, d, k, c_, w_cancel, load_m)
          E_T_2 = E_T_2_k_c(task_t, task_dist_m, k, c_, load_m)
          E_C_2 = E_C_2_k_c(task_t, task_dist_m, k, c_, load_m)
          stdev_T = 0 # math.sqrt(max(E_T_2 - E_T**2, 0) )
          stdev_C = 0 # math.sqrt(max(E_C_2 - E_C**2, 0) )
          coeffvar_T = 0 # stdev_T/E_T
          coeffvar_C = 0 # stdev_C/E_C
          
          if c_:
            xerr_l.append(stdev_T)
            yerr_l.append(stdev_C)
          else:
            xerr_l.append(0)
            yerr_l.append(0)
          if first_moment:
            x_l.append(E_T)
            y_l.append(E_C)
          elif second_moment:
            x_l.append(E_T_2)
            y_l.append(E_C_2)
          elif stdev:
            x_l.append(stdev_T)
            y_l.append(stdev_C)
          elif coeffvar:
            x_l.append(coeffvar_T)
            y_l.append(coeffvar_C)
          if c_ == 0:
            if first_moment:
              if task_t == "SExp":
                xy, xytext = (E_T, E_C), (E_T-0.7, E_C+9)
              elif task_t == "Pareto":
                # plot.axhline(y=E_C, color='k', alpha=0.4, ls='--')
                xy = (E_T, E_C)
                if a == 1.2: xytext = (E_T-20, E_C+3)
                elif a == 1.3: xytext = (E_T-12, E_C+6)
                elif a == 1.5: xytext = (E_T-6, E_C+3)
                elif a == 2: xytext = (E_T-1.2, E_C+3)
                elif a == 3: xytext = (E_T-1, E_C+10)
                else: xytext = (E_T-0.5, E_C+10) # (E_T+0.01, E_C+3)
              # print("*** AMK ***")
              plot.annotate('No\nredundancy\n $c=0$, $n={}$'.format(K), xy=xy, xytext=xytext, fontsize=14)
              # plot.annotate('No redundancy', xy=xy, xytext=xytext)
              # print('>> stdev_T= {}, stdev_C= {}'.format(stdev_T, stdev_C) )
              plot.errorbar([E_T], [E_C], xerr=[stdev_T], yerr=[stdev_C], color='black')
              # plot.errorbar([E_T], [E_C], xerr=[E_T-1], yerr=[stdev_C], color='black')
            elif coeffvar:
              if task_t == "SExp":
                xy, xytext = (coeffvar_T, coeffvar_C), (coeffvar_T-0.075, coeffvar_C-0.022)
              elif task_t == "Pareto":
                xy, xytext = (coeffvar_T, coeffvar_C), (coeffvar_T-0.5, coeffvar_C)
              plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=xy, xytext=xytext)
          elif c_ > 0:
            if task_t == "SExp":
              if first_moment:
                plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C), color=color, fontsize=16)
              elif coeffvar:
                plot.annotate(r'$c={}$'.format(c_), xy=(coeffvar_T, coeffvar_C), xytext=(coeffvar_T+0.003, coeffvar_C), color=color, fontsize=16)
            elif task_t == "Pareto":
              if first_moment:
                if a == 1.22:
                  plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+1, E_C+1), color=color, fontsize=16)
                elif a == 1.33:
                  plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+1, E_C+1), color=color, fontsize=16)
                elif a == 1.5:
                  plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.3, E_C+0.3), color=color, fontsize=16)
                elif a == 2:
                  plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.3, E_C+0.3), color=color, fontsize=16)
                elif a == 2.5:
                  plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.05, E_C+5), color=color, fontsize=16)
                elif a == 3:
                  plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.5, E_C+1), color=color, fontsize=16)
              elif coeffvar:
                plot.annotate(r'$c={}$'.format(c_), xy=(coeffvar_T, coeffvar_C), xytext=(coeffvar_T+0.003, coeffvar_C), color=color)
        if sim:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c_, num_run)
          T_l = stat_id__trial_sampleavg_l_m['T']
          E_T = sum(T_l)/len(T_l)
          stdev_T = math.sqrt(sum([(T - E_T)**2 for T in T_l] )/len(T_l) )
          E_T_2 = sum(stat_id__trial_sampleavg_l_m['T2'] )/len(stat_id__trial_sampleavg_l_m['T2'] )
          
          C_l = stat_id__trial_sampleavg_l_m['C_wc']
          E_C = sum(C_l)/len(C_l)
          stdev_C = math.sqrt(sum([(C - E_C)**2 for C in C_l] )/len(C_l) )
          E_C_2 = sum(stat_id__trial_sampleavg_l_m['C2_wc'] )/len(stat_id__trial_sampleavg_l_m['C2_wc'] )
          
          coeffvar_T = stdev_T/E_T
          coeffvar_C = stdev_C/E_C
          if task_t == "Google":
            E_T, E_C = E_T/60/60, E_C/60/60
          if first_moment:
            x_sim_l.append(E_T)
            y_sim_l.append(E_C)
          elif second_moment:
            x_sim_l.append(E_T_2)
            y_sim_l.append(E_C_2)
          elif stdev:
            x_sim_l.append(stdev_T)
            y_sim_l.append(stdev_C)
          elif coeffvar:
            x_sim_l.append(coeffvar_T)
            y_sim_l.append(coeffvar_C)
          
          if task_t == "Google":
            xy = (E_T, E_C)
            if c_ == 0:
              if k == 400: xytext = (E_T-0.4, E_C+2.5)
              elif k == 1050: xytext = (E_T-0.3, E_C+8)
              elif k == 15: xytext = (E_T-0.04, E_C+0.05)
              plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=xy, xytext=xytext)
            else:
              if k == 400: xytext = (E_T+0.025, E_C+1)
              elif k == 1050: xytext = (E_T+0.025, E_C+1)
              elif k == 15: xytext = (E_T+0.003, E_C+0.015)
              plot.annotate('$c={}$'.format(c_), color=color, xy=xy, xytext=xytext)
      if sim:
        plot.plot(x_sim_l, y_sim_l, label='Simulation, replication', color=color, marker=next(marker), ls=':', mew=0.1)
      else:
        plot.plot(x_l, y_l, label='Replication', color=color, marker='^', zorder=0, ms=10, mew=0.1, ls=':')
        # plot.errorbar(x_l, y_l, xerr=xerr_l, yerr=yerr_l, label='Replication', color=color, marker='^', zorder=0, ms=6, mew=0.1, ls=':')
    elif n:
      # if first_moment and load_m is not None:
      #   n_ = n_max_before_ETpain(load_m['ro_0'], load_m['a_0'], k)
      #   E_T = E_T_k_l_n(task_t, task_dist_m, d, k, k, n_, load_m)
      #   E_C = E_C_k_l_n(task_t, task_dist_m, d, k, k, n_, w_cancel, load_m)
      #   # plot.plot([E_T], [E_C], label=r'$n_{\max}$', color=next(dark_color), zorder=2, marker=next(marker), ls=':', ms=10, mew=3)
      #   # plot.plot([E_T], [E_C], label=r'$Max redundancy$', color=next(dark_color), zorder=2, marker=next(marker), ls=':', ms=10, mew=3)
      # counter = 0
      # for n_ in [*numpy.arange(k, k+4, 1), *numpy.linspace(k+4, 2*k-1, 15), *numpy.arange(2*k, n+1, k) ]:
      #   n_ = int(n_)
      for n_ in numpy.arange(k, n+1, 1):
        if task_t != "Google":
          E_T = E_T_k_l_n(task_t, task_dist_m, d, k, k, n_, load_m)
          E_C = E_C_k_l_n(task_t, task_dist_m, d, k, k, n_, w_cancel, load_m)
          E_T_2 = E_T_2_k_n(task_t, task_dist_m, k, n_, load_m)
          E_C_2 = E_C_2_k_n(task_t, task_dist_m, k, n_, load_m)
          stdev_T =  0 # math.sqrt(max(E_T_2 - E_T**2, 0) ) # 0
          stdev_C = 0 # math.sqrt(max(E_C_2 - E_C**2, 0) )
          coeffvar_T = 0 # stdev_T/E_T
          coeffvar_C = 0 # stdev_C/E_C
          
          if n_ != k:
            xerr_l.append(stdev_T)
            yerr_l.append(stdev_C)
          else:
            xerr_l.append(0)
            yerr_l.append(0)
          if first_moment:
            x_l.append(E_T)
            y_l.append(E_C)
          elif second_moment:
            x_l.append(E_T_2)
            y_l.append(E_C_2)
          elif stdev:
            x_l.append(stdev_T)
            y_l.append(stdev_C)
          elif coeffvar:
            x_l.append(coeffvar_T)
            y_l.append(coeffvar_C)
          # if n_ == k:
          #   plot.annotate('No \nredundancy', xy=(E_T, E_C), xytext=(E_T-75, E_C+30) )
        if sim:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n_, num_run)
          T_l = stat_id__trial_sampleavg_l_m['T']
          E_T = sum(T_l)/len(T_l) # /60/60
          stdev_T = math.sqrt(sum([(T - E_T)**2 for T in T_l] )/len(T_l) )
          E_T_2 = sum(stat_id__trial_sampleavg_l_m['T2'] )/len(stat_id__trial_sampleavg_l_m['T2'] )
          
          C_l = stat_id__trial_sampleavg_l_m['C_wc']
          E_C = sum(C_l)/len(C_l)
          stdev_C = math.sqrt(sum([(C - E_C)**2 for C in C_l] )/len(C_l) )
          E_C_2 = sum(stat_id__trial_sampleavg_l_m['C2_wc'] )/len(stat_id__trial_sampleavg_l_m['C2_wc'] )
          
          coeffvar_T = stdev_T/E_T
          coeffvar_C = stdev_C/E_C
          if task_t == "Google":
            E_T, E_C = E_T/60/60, E_C/60/60
          if first_moment:
            x_sim_l.append(E_T)
            y_sim_l.append(E_C)
          elif second_moment:
            x_sim_l.append(E_T_2)
            y_sim_l.append(E_C_2)
          elif stdev:
            x_sim_l.append(stdev_T)
            y_sim_l.append(stdev_C)
          elif coeffvar:
            x_sim_l.append(coeffvar_T)
            y_sim_l.append(coeffvar_C)
          if task_t == "Google":
            if n_ != k and n_ % k == 0:
              # plot.plot(E_T, E_C, 'x', color="blue", zorder=3, mew=0.1, ms=6)
              xy = (E_T, E_C)
              if k == 400: xytext = (E_T+0.02, E_C+1)
              elif k == 1050: xytext = (E_T+0.02, E_C+1)
              elif k == 15: xytext = (E_T+0.003, E_C+0.015)
              plot.annotate(r'$n={}$'.format(n_), xy=xy, xytext=xytext, color=color)
        if E_T is None:
          continue
        if task_t == "SExp":
          if first_moment:
            if n_ > K and n_ < k+3:
              plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C), color=color, fontsize=16)
            elif n_ != k and n_ % k == 0:
              plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
              plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C), color=color, fontsize=16)
              # plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), ha='left', va='center', xytext=(E_T+0.25, E_C+0.25), color=color)
          elif coeffvar:
            if n_ != k and n_ % k == 0:
              plot.plot(coeffvar_T, coeffvar_C, 'x', color="blue", zorder=3, ms=10, mew=3)
              plot.annotate(r'$n={}$'.format(n_), xy=(coeffvar_T, coeffvar_C), xytext=(coeffvar_T+0.003, coeffvar_C), color=color, fontsize=16)
        elif task_t == "Pareto":
          if first_moment:
            if a == 1.2:
              if n_ > K and n_ < k+3:
                # plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+2, E_C-2), color=color, fontsize=16)
                # elif n_ == k+3:
                #   plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                #   plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-4, E_C-6), color=color)
              elif n_ == 30:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-12, E_C-8), color=color, fontsize=16)
              elif n_ != k and n_ % k == 0:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-12, E_C-3), color=color, fontsize=16)
            elif a == 1.3:
              if n_ > K and n_ < k+3:
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C-5), color=color, fontsize=16)
              elif n_ != k and n_ % k == 0:
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-10, E_C), color=color, fontsize=16)
            elif a == 1.5:
              if n_ > K and n_ < k+3:
                # plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+2, E_C-7), color=color, fontsize=16)
              # elif n_ == k+3:
              #   plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
              #   plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-3.5, E_C-5), color=color)
              elif n_ == 40 or n_ == 60:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-4, E_C+1), color=color, fontsize=16)
              elif n_ != k and n_ % k == 0:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-4, E_C-2), color=color, fontsize=16)
            elif a == 2:
              if n_ > K and n_ < k+3:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T, E_C-5), color=color, fontsize=16)
                # elif n_ == k+3:
                #   plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                #   plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-1.5, E_C-5), color=color)
              elif n_ != k and n_ % k == 0:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-2, E_C-2), color=color, fontsize=16)
            elif a == 3:
              if n_ != k and n_ % k == 0:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-1.5, E_C), color=color, fontsize=16)
                # plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), ha='left', va='center', xytext=(E_T+0.25, E_C), color=color)
            elif a == 2.5:
              if n_ != k and n_ % k == 0:
                plot.plot(E_T, E_C, 'x', color="blue", zorder=3, ms=10, mew=3)
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T*0.78, E_C), color=color, fontsize=16)
                # plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), ha='left', va='center', xytext=(E_T+0.25, E_C), color=color)
          elif coeffvar:
            if n_ != k and n_ % k == 0:
              plot.plot(coeffvar_T, coeffvar_C, 'x', color="blue", zorder=3, ms=10, mew=3)
              plot.annotate(r'$n={}$'.format(n_), xy=(coeffvar_T, coeffvar_C), xytext=(coeffvar_T+0.003, coeffvar_C), color=color, fontsize=16)
      if sim:
        plot.plot(x_sim_l, y_sim_l, label='Simulation, coding', color=color, zorder=1, marker='.', ls=':', mew=1)
      else:
        plot.plot(x_l, y_l, label='Coding', color=color, zorder=1, marker='o', ms=8, mew=0.1, ls=':')
        # print("xerr_l= {}".format(pprint.pformat(xerr_l) ) )
        # print("yerr_l= {}".format(pprint.pformat(yerr_l) ) )
        # plot.errorbar(x_l, y_l, xerr=xerr_l, yerr=yerr_l, label='Coding', color=color, zorder=1, marker='o', ms=5, mew=0.1, ls=':')
  def plot_ET_min(k, n, c):
    x_l, y_l = [], []
    for c_ in range(c+1):
      E_T = E_T_k_c(task_t, task_dist_m, d, k, c_, load_m)
      E_C = E_C_k_c(task_t, task_dist_m, d, k, c_, w_cancel, load_m)
      x_l.append(E_T)
      y_l.append(E_C)
      
      if c_ == 0:
        plot.axhline(y=E_C, color='k', alpha=0.4, ls='--')
        xy = (E_T, E_C)
        if a == 1.5: xytext = (E_T-5, E_C+2)
        else: xytext = (E_T, E_C)
        # plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=xy, xytext=xytext)
        plot.annotate('No redundancy'.format(K), xy=xy, xytext=xytext, fontsize=14)
    color = next(dark_color)
    plot.plot(x_l, y_l, label='Replication', color=color, zorder=1, marker=next(marker), ms=10, mew=0.1, ls=':')
    c_max = 1
    E_Tmin = E_T_k_c(task_t, task_dist_m, d, k, c_max, load_m)
    E_C = E_C_k_c(task_t, task_dist_m, d, k, c_max, w_cancel, load_m)
    plot.annotate(r'$c={}$'.format(c_max), xy=(E_Tmin, E_C), xytext=(E_Tmin+0.3, E_C+0.3), color=color, fontsize=16)
    #
    x_l, y_l = [], []
    for n_ in numpy.arange(k, n+1, 1):
      E_T = E_T_k_l_n(task_t, task_dist_m, d, k, k, n_, load_m)
      E_C = E_C_k_l_n(task_t, task_dist_m, d, k, k, n_, w_cancel, load_m)
      x_l.append(E_T)
      y_l.append(E_C)
    color = next(dark_color)
    plot.plot(x_l, y_l, label='Coding', color=color, zorder=1, marker=next(marker), ms=10, mew=0.1, ls=':')
    n_max = 21
    E_Tmin = E_T_k_l_n(task_t, task_dist_m, d, k, k, n_max, load_m)
    E_C = E_C_k_l_n(task_t, task_dist_m, d, k, k, n_max, w_cancel, load_m)
    plot.annotate(r'$n={}$'.format(n_max), xy=(E_Tmin, E_C), xytext=(E_Tmin+0.5, E_C+0.5), color=color, fontsize=16)
    
    f = plot.gca()
    f.axes.xaxis.set_ticklabels([] )
    f.axes.yaxis.set_ticklabels([] )
  
  plot_(c=9)
  # plot_(c=4)
  plot_(n=10*K)
  
  # plot_ET_min(k=K, n=3*K, c=2)
  
  # plot_(c=5, sim=True)
  # plot_(n=6*K, sim=True)
  #
  # plot.xscale('log')
  # plot.yscale('log')
  
  plot.legend(framealpha=0.25, fontsize=14, numpoints=1) # , loc='lower right'
  
  plot.title(r'$k= {}, {}$'.format(K, task_t_in_latex), fontsize=18)
  if first_moment:
    # leg_T, leg_C = "$E[T]$", "$E[C]$"
    leg_T, leg_C = "Latency", "Cost"
  elif second_moment:
    leg_T, leg_C = "$\mathbb{E}[T^2]$", "$\mathbb{E}[C^2]$"
  elif stdev:
    leg_T, leg_C = "Standard deviation of $T$", "Standard deviation of $C$"
  elif coeffvar:
    leg_T, leg_C = "Coefficient of variance for $T$", "Coefficient of variance for $C$"
  
  plot.xlabel(r'{}'.format(leg_T), fontsize=20)
  plot.ylabel(r'{}'.format(leg_C), fontsize=20)
  plot.gcf().tight_layout()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.1, def_size[1]/1.1)
  prettify(plot.gca() )
  
  # plot.ylim(top=205)
  if task_t == "Pareto":
    plot.savefig("plot_zerodelay_reped_vs_coded_{}_k_{}_a_{}.png".format(task_t, K, a) )
  else:
    plot.savefig("plot_zerodelay_reped_vs_coded_{}_k_{}.png".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_reduct_in_ET_atnocost_wred():
  loc = 3
  d = 0
  def reduction_in_E_T_for_same_E_C_w_red(red_type, loc, a, k):
    E_T_wo_red, E_T = 0, 0
    if red_type == 'coded':
      E_T_wo_red = E_T_pareto_k_n(loc, a, d, k, n=k)
      E_C_wo_red = E_C_pareto_k_n_wrelaunch(loc, a, d, k, n=k)
      
      n_ = k+1
      while 1:
        E_C = E_C_pareto_k_n_wrelaunch(loc, a, d, k, n_)
        # print("n_= {}, E_C_wo_red= {}, E_C= {}".format(n_, E_C_wo_red, E_C) )
        if math.isnan(E_C):
          return reduction_in_E_T_w_red_atnocost__approx(red_type, loc, a, k)
        elif E_C >= E_C_wo_red:
          # print("breaking at n_= {}".format(n_) )
          break
        # print("n_= {}".format(n_) )
        n_ += 1
      E_T = E_T_pareto_k_n(loc, a, d, k, n_-1)
    elif red_type == 'reped':
      E_T_wo_red = E_T_pareto_k_c(loc, a, d, k, c=0)
      E_C_wo_red = E_C_pareto_k_c(loc, a, d, k, c=0)
      
      c_ = 1
      while 1:
        E_C = E_C_pareto_k_c(loc, a, d, k, c_)
        # print("c_= {}, E_C_wo_red= {}, E_C= {}".format(c_, E_C_wo_red, E_C) )
        if math.isnan(E_C) or math.isinf(E_C):
          return None
        elif E_C >= E_C_wo_red:
          # print("breaking at c_= {}".format(c_) )
          break
        c_ += 1
      E_T = E_T_pareto_k_c(loc, a, d, k, c_-1)
    # return E_T_wo_red - E_T
    return (E_T_wo_red - E_T)/E_T_wo_red
  
  def reduction_in_E_T_w_red_atnocost__approx(red_type, loc, a, k):
    if red_type == 'coded':
      E_T_wo_red = E_T_pareto_k_n(loc, a, d, k, n=k)
      # return E_T_wo_red - loc*a
      return max(E_T_wo_red - loc*a, 0)/E_T_wo_red
    elif red_type == 'reped':
      E_T_wo_red = E_T_pareto_k_c(loc, a, d, k, c=0)
      # E_T = loc*math.factorial(k)*G(1/a)/G(k+1/a)
      c_ = max(math.floor(1/(a-1) - 1), 0)
      E_T = E_T_pareto_k_c(loc, a, d, k, c_) # exact
      return max(E_T_wo_red - E_T, 0)/E_T_wo_red
  
  # def threshold_on_tail_for_latency_reduction_at_no_cost(red_type, a):
  #   if red_type == 'reped':
  #     return 1.5
  #   elif red_type == 'coded':
  #     return (a + math.factorial(k)*G(1-1/a)/G(k+1-1/a) )**(1/k)
  
  def plot_(red_type, k):
    x_l, y_l, y_approx_l = [], [], []
    for a in numpy.arange(1.05, 3.5, 0.05):
      x_l.append(a)
      y_l.append(reduction_in_E_T_for_same_E_C_w_red(red_type, loc, a, k) )
      # y_approx_l.append(reduction_in_E_T_w_red_atnocost__approx(red_type, loc, a, k) )
    # plot.axvline(x=threshold_on_tail_for_latency_reduction_at_no_cost(red_type),
    #             color=color, alpha=0.5, ls='--')
    legend = "Rep" if red_type == 'reped' else "Coding"
    plot.plot(x_l, y_l, label=r'{}, $k={}$'.format(legend, k), color=next(dark_color), marker=next(marker), zorder=0, ms=10, mew=0.1, ls=':')
    # plot.plot(x_l, y_approx_l, label=r'{},$k={}$, approx'.format(legend, k), color=next(dark_color), marker=next(marker), mew=0.1, zorder=1, ls='')
  
  # plot_('reped', k=1)
  # plot_('coded', k=1)
  
  # plot_('reped', k=2)
  # plot_('coded', k=2)
  
  plot_('reped', k=10)
  plot_('coded', k=10)
  plot_('reped', k=50)
  plot_('coded', k=50)
  
  plot.legend(fontsize=16, framealpha=0.25, numpoints=1)
  plot.title(r'${} \sim {}(s={}, \alpha)$'.format(r'\mathrm{Tasks}', r'\mathrm{Pareto}', loc), fontsize=18)
  plot.xlabel(r'Tail index $\alpha$', fontsize=18)
  # plot.ylabel('\n'.join(textwrap.wrap(r'Maximum relative latency reduction at no cost', 35) ), fontsize=18)
  plot.ylabel('\n'.join(textwrap.wrap(r'Maximum relative latency reduction at no cost', 25) ), fontsize=18)
  # fig = plot.gcf()
  # def_size = fig.get_size_inches()
  plot.gcf().set_size_inches(8.5, 4)
  prettify(plot.gca() )
  plot.savefig("plot_reduct_in_ET_atnocost_wred.png", bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done.")

# *************************************  (k, n/c, \Delta) with Relaunch  ************************************* #
def plot_arepeat_Pr_T_g_t_pareto_k_wrelaunch():
  K = 100
  loc, a = 3, 2
  
  def plot_(k, d, label=None):
    t_l, y_l, y_approx_l = [], [], []
    # for t in numpy.linspace(loc, loc*40, 100):
    for t in numpy.linspace(loc*20, loc*40, 30):
      t_l.append(t)
      y_l.append(Pr_T_g_t_pareto_k_wrelaunch(loc, a, d, k, t) )
      # y_approx_l.append(Pr_T_g_t_pareto_k_wrelaunch_approx(loc, a, d, k, t) )
    if label is None:
      label = r'$\Delta={}$'.format(d)
    plot.plot(t_l, y_l, label=label, color=next(dark_color), marker=next(marker), ls=':', mew=0.1)
    # plot.plot(t_l, y_approx_l, label=r'Approx, $\Delta={}$'.format(d), color=next(light_color), marker=next(marker), ls=':', mew=0.1)
  
  plot_(K, d=0)
  # plot_(K, d=6*loc)
  # plot_(K, d=45)
  d = Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, K)
  plot_(K, d, label='Minimum $E[T]$, $\Delta= {0:.2f}$'.format(d) )
  d = Delta_for_min_Pr_T_g_t_pareto_k_wrelaunch(loc, a, K)
  plot_(K, d, label='Minimum {}, $\Delta= {}$'.format(r'$Pr\{T > t\}$', '{0:.2f}'.format(d) ) )
  
  # def plot_(k, t):
  #   d_l, y_l = [], []
  #   for d in numpy.linspace(0, 2*t, 200):
  #     d_l.append(d)
  #     y_l.append(Pr_T_g_t_pareto_k_wrelaunch(loc, a, d, k, t) )
  #   plot.plot(d_l, y_l, label=r'$t={0:.2f}$'.format(t), color=next(dark_color), marker=next(marker), ls=':', mew=0.1)
  
  # plot_(K, t=10*loc)
  # plot_(K, t=20*loc)
  
  plot.legend()
  plot.title(r'$X \sim Pareto(\lambda={}, \alpha={}), k= {}$'.format(loc, a, K) )
  plot.xlabel(r'$t$', fontsize=13)
  # plot.xlabel(r'$\Delta$', fontsize=13)
  plot.ylabel(r'$Pr\{T > t\}$', fontsize=13)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_arepeat_Pr_T_g_t_pareto_k_wrelaunch_k_{}.png".format(K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_k_nc_wrelaunch():
  K = 100
  D, mu = 30, 1
  loc, a = 3, 2 # 100 # 1.2 # 2
  task_t = "Pareto"
  
  def plot_reduction_in_ET_wrelaunch(k, loc=loc):
    def max_a_lb():
      return math.log(k)/math.log(4)
    
    x_l, y_l = [], []
    for a in numpy.linspace(1.1, 9.5, 50):
      x_l.append(a)
      
      E_T_wrelaunch_min = float('Inf')
      # for d in numpy.linspace(loc+0.05, 10*loc, 1000):
      for d in numpy.linspace(loc+0.01, 5*loc, 100):
        E_T = E_T_pareto_k_n_wrelaunch(loc, a, d, k, n=k)
        if E_T < E_T_wrelaunch_min:
          E_T_wrelaunch_min = E_T
        else:
          # elif not math.isinf(E_T_wrelaunch_min):
          break
      ET_norelaunch = E_T_pareto_k_n(loc, a, d, k, n=k)
      y_l.append(max((ET_norelaunch - E_T_wrelaunch_min)/ET_norelaunch, 0) )
    c = next(dark_color)
    label = r'$k={}$'.format(k) # r'$\lambda={}$, $k={}$'.format(loc, k)
    plot.plot(x_l, y_l, label=label, color=c, marker=next(marker), ms=11, mew=0.1, ls=':')
    # plot.axvline(x=max_a_lb(), label=r'$\alpha_u$, {}'.format(label), color=c, ls='--')
    plot.axvline(x=max_a_lb(), color=c, ls='--', lw=2.5)
  
  E_T = True # False
  w_cancel = True # False
  num_srun = 10000
  def plot_(k, n=None, c=None, sim=False):
    task_t_rv = Pareto(loc, a)
    
    l = k
    x_l = []
    y_wrelaunch_l, y_wrelaunch_approx_l, y_sim_wrelaunch_l = [], [], []
    y_norelaunch_l = []
    for d in numpy.linspace(0, 20*loc, 50):
      x_l.append(d)
      
      if E_T:
        if n is not None:
          y_wrelaunch_l.append(E_T_pareto_k_n_wrelaunch(loc, a, d, k, n) )
          y_wrelaunch_approx_l.append(E_T_pareto_k_n_wrelaunch(loc, a, d, k, n) )
          if n == k: y_norelaunch_l.append(E_T_pareto_k_n(loc, a, d, k, n) )
        elif c is not None:
          y_wrelaunch_l.append(E_T_pareto_k_c_wrelaunch(loc, a, d, k, c) )
          # y_wrelaunch_approx_l.append(E_T_pareto_k_c_wrelaunch_approx(loc, a, d, k, c) )
          if c == 0: y_norelaunch_l.append(E_T_pareto_k_c(loc, a, d, k, c) )
      else:
        if n is not None:
          y_wrelaunch_l.append(E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel) )
          # y_wrelaunch_approx_l.append(E_C_pareto_k_n_approx(loc, a, d, k, n, w_cancel=w_cancel) )
        elif c is not None:
          y_wrelaunch_l.append(E_C_pareto_k_c_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel) )
          # y_wrelaunch_approx_l.append(E_C_pareto_k_c_approx(loc, a, d, k, n, w_cancel=w_cancel) )
      
      if sim:
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=num_srun, w_relaunch=True)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=num_srun, w_relaunch=True)
        if E_T:
          y_sim_wrelaunch_l.append(sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] ) )
        else:
          y_sim_wrelaunch_l.append(
            sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] ) if w_cancel else \
            sum(stat_id__trial_sampleavg_l_m['C'] )/len(stat_id__trial_sampleavg_l_m['C'] ) )
            
    color = next(dark_color)
    if n == k:
      d = Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, k)
      # plot.axvline(x=d, label=r'$\Delta^*$, n={}'.format(n), color=color)
      y = E_T_pareto_k_n_wrelaunch(loc, a, d, k, n)
      plot.plot(d, y, color='red', label=r'Approx optimal', marker='*', zorder=1, ms=10, mew=3)
    # print("y_wrelaunch_l= {}".format(y_wrelaunch_l) )
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    # label = 'With relaunch'
    # plot.plot(x_l, y_wrelaunch_l, label='Simulation, {}'.format(label), color=color, marker=next(marker), zorder=0, ms=10, mew=0.1, ls='--')
    # plot.plot(x_l, y_wrelaunch_approx_l, label='Approximate, {}'.format(label), color=next(light_color), marker='x', ms=8, mew=2, ls='')
    plot.plot(x_l, y_wrelaunch_l, label='Approximate, {}'.format(label), color=color, marker=next(marker), zorder=0, ms=10, mew=0.1, ls='--')
    if sim:
      plot.plot(x_l, y_sim_wrelaunch_l, label='Simulation, {}'.format(label), color=next(light_color), marker='x', ms=8, mew=2, ls='')
    
    if E_T:
      if n == k: plot.plot(x_l, y_norelaunch_l, label='No relaunch', color=next(dark_color), ms=10, mew=0.1, ls='--')
  
  sim = True # False
  K = 100
  # plot_reduction_in_ET_wrelaunch(K)
  # plot_reduction_in_ET_wrelaunch(10*K)
  # plot_reduction_in_ET_wrelaunch(100*K)
  # plot_reduction_in_ET_wrelaunch(1000*K)
  
  # plot_reduction_in_ET_wrelaunch(1, 10*K)
  # plot_reduction_in_ET_wrelaunch(loc, 10*K)
  # plot_reduction_in_ET_wrelaunch(10*loc, 10*K)
  # plot_reduction_in_ET_wrelaunch(100*loc, 10*K)
  
  # plot_(K, n=int(K*1.1), sim=sim)
  # plot_(K, n=int(K*1.4), sim=sim)
  # plot_(K, n=2*K, sim=sim)
  
  # plot_(K, n=K)
  # plot_(K, n=K+1, sim=sim)
  # plot_(K, n=K+2, sim=sim)
  # plot_(K, n=K+3, sim=sim)
  
  # plot_(K, c=0, sim=sim)
  plot_(K, c=1, sim=sim)
  plot_(K, c=2, sim=sim)
  # plot_(K, c=3, sim=sim)
  
  plot.legend(loc='lower right', fontsize=14, framealpha=0.25, numpoints=1)
  # plot.legend(loc='upper right', fontsize=14, numpoints=1)
  # plot.legend(loc='center left', bbox_to_anchor=(0.7, 0.5) )
  plot.title(r'${} \sim {}(s= {}, \alpha= {}), k= {}$'.format('\mathrm{Tasks}', '\mathrm{Pareto}', loc, a, K), fontsize=18)
  # plot.title(r'${} \sim {}(s= {}, \alpha)$'.format('\mathrm{Tasks}', '\mathrm{Pareto}', loc), fontsize=18)
  plot.xlabel(r'$\Delta$', fontsize=18)
  plot.ylabel('Latency' if E_T else 'Cost', fontsize=18)
  # plot.xlabel(r'Tail index $\alpha$', fontsize=18)
  # plot.ylabel('\n'.join(textwrap.wrap(r'Maximum relative reduction in latency with relaunch', 40) ), fontsize=18)
  # fig = plot.gcf()
  # fig.tight_layout()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  prettify(plot.gca() )
  plot.savefig("plot_k_nc_wrelaunch_{}_k_{}.png".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_k_nc_wrelaunch_EC_vs_ET():
  K = 100 # 15 # 1050 # 100
  loc, a = 3, 2 # 1.2
  w_cancel = True
  
  task_t = "Pareto" # "Google"
  if task_t == "Pareto": task_t_in_latex = r'\sim {}(s={}, \alpha={})'.format(r'\mathrm{Pareto}', loc, a)
  elif task_t == "Google": task_t_in_latex = r'\sim {}'.format(r'\mathrm{Google}',)
  
  def plot_wrelaunch_vs_worelaunch(k, n=None, c=None):
    task_dist_m = {'loc': loc, 'a': a}
    E_T_l, E_C_l = [], []
    E_T_wrelaunch_l, E_C_wrelaunch_l = [], []
    
    d = Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, k)
    if n is not None:
      for n_ in range(k, n+1):
        E_T_l.append(E_T_k_l_n(task_t, task_dist_m, 0, k, k, n_) )
        E_C_l.append(E_C_k_l_n(task_t, task_dist_m, 0, k, k, n_, w_cancel) )
        
        E_T_wrelaunch_l.append(E_T_pareto_k_n_wrelaunch(loc, a, d, k, n_) )
        E_C_wrelaunch_l.append(E_C_pareto_k_n_wrelaunch(loc, a, d, k, n_, w_cancel) )
    elif c is not None:
      for c_ in range(0, c+1):
        E_T_l.append(E_T_k_c(task_t, task_dist_m, 0, k, c_) )
        E_C_l.append(E_C_k_c(task_t, task_dist_m, 0, k, c_, w_cancel) )
        
        E_T_wrelaunch_l.append(E_T_pareto_k_c_wrelaunch(loc, a, d, k, c_) )
        E_C_wrelaunch_l.append(E_C_pareto_k_c_wrelaunch(loc, a, d, k, c_, w_cancel) )
    # plot.plot(E_T_l, E_C_l, label=r'w/o relaunch', color=next(dark_color), marker=next(marker), zorder=0, ls=':', mew=mew)
    # print("E_T_wrelaunch_l= {}".format(pprint.pformat(E_T_wrelaunch_l) ) )
    # print("E_C_wrelaunch_l= {}".format(pprint.pformat(E_C_wrelaunch_l) ) )
    plot.plot(E_T_wrelaunch_l, E_C_wrelaunch_l, label=r'w/ relaunch', color=next(dark_color), marker=next(marker), zorder=0, ls=':', mew=mew)
  
  d_ul = 10**5
  def plot_(k, n=None, c=None, sim=False):
    num_f = 5
    num_run = 10000 # 10000
    E_T_l, E_C_l = [], []
    E_T_sim_l, E_C_sim_l = [], []
    
    if task_t == "Pareto": task_t_rv = Pareto(loc, a)
    elif task_t == "Google": task_t_rv = Google(k)
    
    # for d in [*numpy.linspace(0, 5*loc, 50), *numpy.linspace(5*loc, 200*loc, 100)]:
    # for d in numpy.logspace(-3, 2, 100):
    for d in numpy.logspace(-2, 4, 50):
      if sim:
        E_T_sum, E_C_sum = 0, 0
        for i in range(num_f):
          if n is not None:
            stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n, num_run, w_relaunch=True)
          elif c is not None:
            stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run, w_relaunch=True)
          E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
          E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
          if task_t == "Google":
            E_T = E_T/60/60
            E_C = E_C/60/60
          E_T_sum += E_T
          E_C_sum += E_C
        E_T = E_T_sum/num_f
        E_C = E_C_sum/num_f
        E_T_sim_l.append(E_T)
        E_C_sim_l.append(E_C)
      if task_t != "Google":
        if n is not None:
          E_T_l.append(E_T_pareto_k_n_wrelaunch(loc, a, d, k, n) )
          E_C_l.append(E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel) )
        elif c is not None:
          E_T_l.append(E_T_pareto_k_c_wrelaunch(loc, a, d, k, c) )
          E_C_l.append(E_C_pareto_k_c_wrelaunch(loc, a, d, k, c, w_cancel) )
    color = next(dark_color)
    if n != k:
      label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    else:
      label = "" # "No redundancy with relaunch"
    plot.plot(E_T_l, E_C_l, label=label, color=color, marker=next(marker), zorder=0, ls=':', ms=10, mew=0.1)
    # plot.plot(E_T_sim_l, E_C_sim_l, label='Simulation, {}'.format(label), color=next(dark_color), marker=next(marker), zorder=0, ls=':', mew=0.1)
    
    if n is not None and n == k:
      if task_t != "Google":
        x = E_T_pareto_k_n_wrelaunch(loc, a, 0, k, n)
        y = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=w_cancel)
        if k == 100 and a == 2:
          # plot.annotate(r'$\Delta=0$ or $\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=12)
          plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+1, y+4), fontsize=16)
          plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x-1.5, y-25), fontsize=16)
          
          d = Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, k)
          x = E_T_pareto_k_n_wrelaunch(loc, a, d, k, n)
          y = E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel)
          plot.plot(x, y, color='red', label=r'Approximate optimal', marker='*', zorder=1, ms=20, mew=0.1, ls='None')
          plot.annotate(r'$\Delta={0:.2f}$'.format(d), xy=(x, y), xytext=(x-1, y+20), fontsize=16)
        else:
          # plot.plot(x, y, color='black', marker='*', ms=10, mew=3)
          plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x, y), fontsize=16)
      else: # task_t = "Google"
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, 0, k, k, k, num_run, w_relaunch=True)
        x = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )/60/60
        y = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )/60/60
        plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+0.005, y+2.5), fontsize=16)
        
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d_ul, k, k, k, num_run, w_relaunch=True)
        x = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )/60/60
        y = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )/60/60
        plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x-0.005, y-2.5), xytext=(x, y), fontsize=16)
  
  label_added = False
  def put_blackcross_atd0(k, n=None, c=None):
    nonlocal label_added
    if n is not None:
      x = E_T_pareto_k_n_wrelaunch(loc, a, 0, k, n)
      y = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=w_cancel)
    elif c is not None:
      x = E_T_pareto_k_c_wrelaunch(loc, a, 0, k, c)
      y = E_C_pareto_k_c_wrelaunch(loc, a, 0, k, c, w_cancel=w_cancel)
    # plot.plot(x, y, color='black', marker='o', zorder=1, ms=10, mew=0.1)
    # if not label_added:
    #   plot.plot(x, y, label=r'$\Delta=0$', color='black', marker='x', zorder=1, ms=11, mew=3)
    # else:
    #   plot.plot(x, y, color='black', marker='x', zorder=1, ms=11, mew=3)
    label = None
    if not label_added:
      label_added = True
      label = r'$\Delta=0$'
    plot.plot(x, y, label=label, color='black', marker='x', zorder=1, ms=11, mew=3, ls='None')
  def put_dtoinfty(k, n=None, c=None, d_ul_=None):
    nonlocal d_ul
    d_ul = d_ul if d_ul_ is None else d_ul_
    if n is not None:
      x = E_T_pareto_k_n_wrelaunch(loc, a, d_ul, k, n)
      y = E_C_pareto_k_n_wrelaunch(loc, a, d_ul, k, n, w_cancel=w_cancel)
      plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x-2.5, y+10), fontsize=16)
    elif c is not None:
      x = E_T_pareto_k_c_wrelaunch(loc, a, d_ul, k, c)
      y = E_C_pareto_k_c_wrelaunch(loc, a, d_ul, k, c, w_cancel=w_cancel)
      plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x-2.5, y+35), fontsize=16)
  
  plot_(K, n=K)
  plot.xlim(([24, 60] ) )
  
  # plot_(K, n=K+2)
  # plot_(K, n=K+1)
  # put_blackcross_atd0(k=K, n=K+2)
  # put_blackcross_atd0(k=K, n=K+1)
  # # plot.plot(x, y, color='black', label=r'$\Delta=0$', marker='o', zorder=1, ms=10, mew=0.1)
  # put_dtoinfty(k=K, n=K+1)
  
  # plot_(K, n=K+10)
  # plot_(K, n=K+20)
  # plot_(K, n=K+30)
  # put_blackcross_atd0(k=K, n=K+10)
  # put_blackcross_atd0(k=K, n=K+20)
  # put_blackcross_atd0(k=K, n=K+30)
  # put_dtoinfty(k=K, n=K+30)
  
  # plot_(K, c=1)
  # plot_(K, c=2)
  # plot_(K, c=3)
  # put_blackcross_atd0(k=K, c=1)
  # put_blackcross_atd0(k=K, c=2)
  # put_blackcross_atd0(k=K, c=3)
  # put_dtoinfty(k=K, c=3)
  
  # plot_wrelaunch_vs_worelaunch(K, n=2*K)
  # plot_wrelaunch_vs_worelaunch(K, c=50)
  
  plot.legend(loc='upper left', fontsize=14, framealpha=0.25, numpoints=1)
  # plot.legend(loc='upper right', fontsize=16, framealpha=0.25, numpoints=1)
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  # plot.xscale('log')
  # plot.yscale('log')
  plot.title(r'${}{}, k= {}$'.format(r'\mathrm{Tasks}', task_t_in_latex, K), fontsize=18)
  plot.xlabel(r'Latency', fontsize=18)
  plot.ylabel('Cost', fontsize=18)
  # fig = plot.gcf()
  # fig.tight_layout()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  prettify(plot.gca() )
  plot.savefig("plot_k_nc_wrelaunch_EC_vs_ET_{}_k_{}.pdf".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

# **********************************  (k, n, 0) with Relaunch at d  ****************************** #
def plot_k_ncd0_wrelaunch():
  K = 10 # 100
  loc, a = 3, 2 # 100 # 1.2 # 2
  task_t = "Pareto"
  
  ET = True # False
  w_cancel = True
  num_srun = 10000
  def plot_(k, n=None, c=None, sim=False):
    x_l = []
    y_wrelaunch_l, y_wrelaunch_approx_l, y_sim_wrelaunch_l = [], [], []
    y_norelaunch_l = []
    # for d in numpy.linspace(0, 20*loc, 50):
    for d in numpy.linspace(0, 10*loc, 50):
    # for d in numpy.linspace(0, 20*loc, 10):
      x_l.append(d)
      
      if ET:
        if n is not None:
          y_wrelaunch_l.append(E_T_pareto_k_nd0_wrelaunch(loc, a, d, k, n) )
          # y_wrelaunch_approx_l.append(E_T_pareto_k_n_approx(loc, a, d, k, n, w_relaunch=True) )
          if n == k: y_norelaunch_l.append(E_T_pareto_k_n(loc, a, d, k, n) )
        elif c is not None:
          y_wrelaunch_l.append(E_T_pareto_k_cd0_wrelaunch(loc, a, d, k, c) )
          # y_wrelaunch_approx_l.append(E_T_pareto_k_c_wrelaunch_approx(loc, a, d, k, c, w_relaunch=True) )
          if c == 0: y_norelaunch_l.append(E_T_pareto_k_c(loc, a, d, k, c) )
      else:
        if n is not None:
          y_wrelaunch_l.append(E_C_pareto_k_nd0_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel) )
          # y_wrelaunch_approx_l.append(E_C_pareto_k_n_approx(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
        elif c is not None:
          y_wrelaunch_l.append(E_C_pareto_k_cd0_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel) )
          # y_wrelaunch_approx_l.append(E_C_pareto_k_c_approx(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
      if sim:
        task_t_rv = Pareto(loc, a)
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_nd0_wrelaunch(task_t_rv, k, n, d, num_run=num_srun)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_cd0_wrelaunch(task_t_rv, k, c, d, num_run=num_srun)
        if ET:
          y_sim_wrelaunch_l.append(sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] ) )
        else:
          y_sim_wrelaunch_l.append(
            sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] ) if w_cancel else \
            sum(stat_id__trial_sampleavg_l_m['C'] )/len(stat_id__trial_sampleavg_l_m['C'] ) )
    # d = Delta_for_min_E_T_pareto_k_n0_wrelaunch(loc, a, k, n)
    # y = E_T_pareto_k_nd0_wrelaunch(loc, a, d, k, n) if ET else E_C_pareto_k_nd0_wrelaunch(loc, a, d, k, n, w_cancel)
    # plot.plot(d, y, color='red', label=r'Approx optimal', marker='*', zorder=1, ms=14, mew=0.1)
    
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    # label = 'With relaunch'
    plot.plot(x_l, y_wrelaunch_l, label=label, color=next(dark_color), marker=next(marker), zorder=0, ms=10, mew=0.1, ls='--')
    # plot.plot(x_l, y_wrelaunch_approx_l, label='approx, {}'.format(label), color=next(dark_color), marker=next(marker), ls='', mew=0.1)
    
    if sim:
      plot.plot(x_l, y_sim_wrelaunch_l, label=r'Simulation, {}'.format(label), color=next(light_color), marker=next(marker), ms=8, mew=0.1, ls='')
    if ET:
      if n == k: plot.plot(x_l, y_norelaunch_l, label='No relaunch', color=next(dark_color), ls='--', mew=0.1)
  
  sim = True # False
  # plot_(K, n=K, sim=sim)
  # plot_(K, n=K+1, sim=sim)
  # plot_(K, n=K+2, sim=sim)
  plot_(K, n=int(K*1.1), sim=sim)
  plot_(K, n=int(K*1.4), sim=sim)
  plot_(K, n=2*K, sim=sim)
  
  # plot_(K, n=K)
  # plot_(K, n=K+1, sim=sim)
  # plot_(K, n=K+2, sim=sim)
  # plot_(K, n=K+3, sim=sim)
  
  # plot_(K, c=0, sim=sim)
  # plot_(K, c=1, sim=sim)
  # plot_(K, c=2, sim=sim)
  # plot_(K, c=3, sim=sim)
  
  plot.legend(fontsize=14, framealpha=0.25, numpoints=1)
  # plot.legend(loc='center left', bbox_to_anchor=(0.7, 0.5) )
  plot.title(r'${} \sim {}(s= {}, \alpha= {}), k= {}$'.format(r'\mathrm{Tasks}', r'\mathrm{Pareto}', loc, a, K) )
  
  plot.xlabel(r'$\Delta$', fontsize=13)
  plot.ylabel('Latency' if ET else 'Cost', fontsize=13)
  # plot.xlabel(r'$\alpha$', fontsize=12)
  # plot.ylabel('\n'.join(textwrap.wrap(r'Maximum percetange reduction in $E[T]$ with relaunch', 40) ), fontsize=12)
  # fig = plot.gcf()
  # fig.tight_layout()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_k_ncd0_wrelaunch_{}_k_{}.png".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_k_ncd0_wrelaunch_EC_vs_ET():
  K = 100
  loc, a = 3, 1.2
  w_cancel = True
  
  task_t = "Pareto" # "Google"
  if task_t == "Pareto": task_t_in_latex = r'\sim {}(s={}, \alpha={})'.format(r'\mathrm{Pareto}', loc, a)
  elif task_t == "Google": task_t_in_latex = r'\sim \mathrm{Google}'
  
  def plot_wrelaunch_vs_worelaunch(k, n=None, c=None):
    task_dist_m = {'loc': loc, 'a': a}
    E_T_l, E_C_l = [], []
    E_T_wrelaunch_l, E_C_wrelaunch_l = [], []
    
    d = Delta_for_min_E_T_pareto_k_n0_wrelaunch(loc, a, k, n)
    if n is not None:
      for n_ in range(k, n+1):
        E_T_l.append(E_T_k_l_n(task_t, task_dist_m, 0, k, k, n_) )
        E_C_l.append(E_C_k_l_n(task_t, task_dist_m, 0, k, k, n_, w_cancel) )
        
        E_T_wrelaunch_l.append(E_T_pareto_k_nd0_wrelaunch(loc, a, d, k, n_) )
        E_C_wrelaunch_l.append(E_C_pareto_k_nd0_wrelaunch(loc, a, d, k, n_, w_cancel) )
    elif c is not None:
      for c_ in range(0, c+1):
        E_T_l.append(E_T_k_c(task_t, task_dist_m, 0, k, c_) )
        E_C_l.append(E_C_k_c(task_t, task_dist_m, 0, k, c_, w_cancel) )
        
        E_T_wrelaunch_l.append(E_T_pareto_k_c0_wrelaunch(loc, a, d, k, c_) )
        E_C_wrelaunch_l.append(E_C_pareto_k_c0_wrelaunch(loc, a, d, k, c_, w_cancel) )
    plot.plot(E_T_l, E_C_l, label=r'w/o relaunch', color=next(dark_color), marker=next(marker), zorder=0, ls=':', mew=mew)
    plot.plot(E_T_wrelaunch_l, E_C_wrelaunch_l, label=r'w/ relaunch', color=next(dark_color), marker=next(marker), zorder=0, ls=':', mew=mew)
  
  def plot_(k, n=None, c=None, sim=False):
    num_f = 5
    E_T_l, E_C_l = [], []
    E_T_sim_l, E_C_sim_l = [], []
    
    if task_t == "Pareto": task_t_rv = Pareto(loc, a)
    elif task_t == "Google": task_t_rv = Google(k)
    
    # for d in [*numpy.linspace(0, 5*loc, 50), *numpy.linspace(5*loc, 200*loc, 100)]:
    # for d in numpy.logspace(-3, 2, 100):
    d_ul = 10**5
    for d in numpy.logspace(-2, 4, 40):
      if sim:
        num_run = 10000
        E_T_sum, E_C_sum = 0, 0
        for i in range(num_f):
          if n is not None:
            stat_id__trial_sampleavg_l_m = sim_arepeat_k_nd0_wrelaunch(task_t_rv, k, n, d, num_run)
          elif c is not None:
            stat_id__trial_sampleavg_l_m = sim_arepeat_k_cd0_wrelaunch(task_t_rv, k, n, d, num_run)
          E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
          E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
          E_T_sum += E_T
          E_C_sum += E_C
        E_T = E_T_sum/num_f
        E_C = E_C_sum/num_f
        E_T_sim_l.append(E_T)
        E_C_sim_l.append(E_C)
      if task_t != "Google":
        if n is not None:
          E_T_l.append(E_T_pareto_k_nd0_wrelaunch(loc, a, d, k, n) )
          E_C_l.append(E_C_pareto_k_nd0_wrelaunch(loc, a, d, k, n, w_cancel) )
        elif c is not None:
          E_T_l.append(E_T_pareto_k_cd0_wrelaunch(loc, a, d, k, c) )
          E_C_l.append(E_C_pareto_k_cd0_wrelaunch(loc, a, d, k, c, w_cancel) )
    color = next(dark_color)
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    plot.plot(E_T_l, E_C_l, label=label, color=color, marker=next(marker), zorder=0, ls=':', mew=0.1)
    if sim:
      plot.plot(E_T_sim_l, E_C_sim_l, label='Simulation, {}'.format(label), color=next(dark_color), marker=next(marker), zorder=0, ls=':', mew=0.1)
    
      if task_t != "Google":
        x = E_T_pareto_k_nd0_wrelaunch(loc, a, 0, k, n)
        y = E_C_pareto_k_nd0_wrelaunch(loc, a, 0, k, n, w_cancel=w_cancel)
        if k == 100 and a == 2:
          plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+1, y+4), fontsize=12)
          plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x-1.5, y-25), fontsize=12)
          
          d = Delta_for_min_E_T_pareto_k_n0_wrelaunch(loc, a, k, n)
          x = E_T_pareto_k_nd0_wrelaunch(loc, a, d, k, n)
          y = E_C_pareto_k_nd0_wrelaunch(loc, a, d, k, n, w_cancel)
          plot.plot(x, y, color='red', label=r'Approx optimal', marker='*', zorder=1, ms=10, mew=3)
          plot.annotate(r'$\Delta={0:.2f}$'.format(d), xy=(x, y), xytext=(x-1.5, y+20), fontsize=12)
        else:
          # plot.plot(x, y, color='black', marker='*', ms=10, mew=3)
          plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x, y), fontsize=12)
      else: # task_t = "Google"
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, 0, k, k, k, num_run, w_relaunch=True)
        x = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
        y = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
        plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+0.005, y+2.5), fontsize=12)
        
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d_ul, k, k, k, num_run, w_relaunch=True)
        x = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
        y = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
        plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x-0.005, y-2.5), xytext=(x, y), fontsize=12)
    # elif n is not None and n > k:
    #   x = E_T_pareto_k_n_wrelaunch(loc, a, 0, k, n)
    #   y = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=w_cancel)
    #   if k == 100 and a == 2:
    #     plot.plot(x, y, color='black', marker='o', zorder=1, ms=10, mew=0.1)
    #     # plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+1, y+5), color=color, fontsize=13)
    #   x = E_T_pareto_k_n_wrelaunch(loc, a, d_ul, k, n)
    #   y = E_C_pareto_k_n_wrelaunch(loc, a, d_ul, k, n, w_cancel=w_cancel)
    #   if n == 101 and k == 100 and a == 2:
    #     plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=13)
    # elif c is not None:
    #   x = E_T_pareto_k_c_wrelaunch(loc, a, 0, k, c)
    #   y = E_C_pareto_k_c_wrelaunch(loc, a, 0, k, c, w_cancel=w_cancel)
    #   plot.plot(x, y, color='black', marker='o', zorder=1, ms=6, mew=0.1)
  # plot_(K, n=K)
  plot_(K, n=K+1)
  # plot_(K, n=K+2)
  
  # plot_wrelaunch_vs_worelaunch(K, n=2*K)
  # plot_wrelaunch_vs_worelaunch(K, c=50)
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  # plot.xscale('log')
  # plot.yscale('log')
  plot.title(r'${}{}, k= {}$'.format(r'\mathrm{Tasks}', task_t_in_latex, K) )
  plot.xlabel(r'$E[T]$', fontsize=13)
  plot.ylabel(r'$E[C]$', fontsize=13)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_k_ncd0_wrelaunch_EC_vs_ET_{}_k_{}.png".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

# *************************************  (k, n/c, \Delta) with Relaunch  ************************************* #
def plot_k_nc_retainl_atd():
  K = 15
  D, mu = 30, 1
  loc, a = 3, 1.5 # 4 # 1.2 # 2
  task_t = "Pareto"
  
  E_T = False # True
  def plot_(k, n=None, c=None, sim=False):
    task_t_rv = Pareto(loc, a)
    
    x_l, x_sim_l = [], []
    y_l, y_approx_l, y_sim_l = [], [], []
    
    # y_ref_l = []
    for d in numpy.linspace(loc+0.05, 20*loc, 100):
      x_l.append(d)
      
      if n is not None:
        if E_T:
          y_l.append(E_T_pareto_k_n_retainl_atd(loc, a, k, n, d) )
          # y_approx_l.append(approx_E_T_pareto_k_n_retainl_atd(loc, a, k, n, d) )
        else:
          y_l.append(E_C_pareto_k_n_retainl_atd(loc, a, k, n, d) )
          
          # y_ref_l.append(E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n) )
      elif c is not None:
        pass
      if sim:
        # if n is not None:
        #   stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=num_srun, w_relaunch=True)
        # elif c is not None:
        #   stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=num_srun, w_relaunch=True)
        # if E_T:
        #   y_sim_wrelaunch_l.append(sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] ) )
        pass
    label = 'c= {}'.format(c) if c is not None else 'n= {}'.format(n)
    plot.plot(x_l, y_l, label=label, color=next(dark_color), marker=next(marker), ls='--', zorder=0, mew=0.1)
    # plot.plot(x_l, y_approx_l, label='approx, {}'.format(label), color=next(dark_color), marker=next(marker), ls='', mew=0.1)
    # plot.plot(x_l, y_ref_l, label='ref, n= {}'.format(n), color=next(dark_color), marker=next(marker), ls='--', zorder=0, mew=0.1)
    if sim:
      plot.plot(x_l, y_sim_l, label=r'sim, {}'.format(label), color=next(light_color), marker=next(marker), ls='', mew=0.1)
  
  sim = False # True
  plot_(K, n=K+20, sim=sim)
  # plot_(K, n=K, sim=sim)
  
  plot.legend()
  plot.title(r'$X \sim Pareto(\lambda= {}, \alpha= {}), k= {}$'.format(loc, a, K) )
  
  plot.xlabel(r'$\Delta$', fontsize=12)
  plot.ylabel(r'$E[T]$' if E_T else r'$E[C]$', fontsize=12)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_k_nc_retainl_atd_{}_k_{}.png".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_k_nc_retainl_atd__C_vs_T():
  K = 15
  D, mu = 30, 1
  loc, a = 3, 2 # 1.5 # 2
  task_t = "Pareto"
  
  def plot_varying_d(k, n):
    x_l, y_l = [], []
    
    d = 1
    for n_ in range(k+1, n+1, 1):
      d_ = (n-k)*d/(n_-k)
      E_T = E_T_pareto_k_n_retainl_atd(loc, a, k, n_, d_)
      E_C = E_C_pareto_k_n_retainl_atd(loc, a, k, n_, d_)
      
      plot.plot([E_T], [E_C], label=r'$n$= {}, $\Delta$= {}'.format(n_, d_), color=next(dark_color), marker=next(marker), ls='--', zorder=2, mew=1)
      
      # x_l.append(E_T)
      # y_l.append(E_C)
    # plot.plot(x_l, y_l, label=r'n= {}, $\Delta_0$= {}'.format(n, d), color=next(dark_color), marker=next(marker), ls='--', zorder=0, mew=1)
  
  def plot_(k, n, sim=False):
    # task_t_rv = Pareto(loc, a)
    
    x_l, x_sim_l = [], []
    y_l, y_sim_l = [], []
    
    # E_T = E_X_n_k_pareto(loc, a, k, k)
    # E_C = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, k)
    # plot.plot([E_T], [E_C], label="n=k", color='red', marker=next(marker), zorder=1, mew=3)
    
    # E_T = E_X_n_k_pareto(loc, a, n, k)
    # E_C = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n)
    # plot.plot([E_T], [E_C], label="n= {}".format(n), color='blue', marker=next(marker), zorder=1, mew=3)
    
    col = next(dark_color)
    mar = next(marker)
    for n_ in range(k, n+1, 1):
      E_T = E_X_n_k_pareto(loc, a, n_, k)
      E_C = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n_)
      # plot.plot([E_T], [E_C], label="n= {}".format(n_), color=next(dark_color), marker=next(marker), zorder=1, mew=3)
      plot.plot([E_T], [E_C], color=col, marker=mar, zorder=1, mew=1)
    
    # for d in numpy.linspace(loc+0.05, 20*loc, 1000):
    for d in numpy.logspace(0, 4, 100):
      x_l.append(E_T_pareto_k_n_retainl_atd(loc, a, k, n, d) )
      y_l.append(E_C_pareto_k_n_retainl_atd(loc, a, k, n, d) )
      if sim:
        pass
    plot.plot(x_l, y_l, label='n= {}'.format(n), color=next(dark_color), marker=next(marker), ls='--', zorder=0, mew=1)
    if sim:
      plot.plot(x_sim_l, y_sim_l, label=r'sim, n= {}'.format(n), color=next(light_color), marker=next(marker), ls='', mew=0.1)
  
  sim = False # True
  plot_(K, n=K+20, sim=sim)
  plot_varying_d(K, n=K+20)
  
  plot.legend(fontsize=14, framealpha=0.25, numpoints=1)
  plot.title(r'$X \sim Pareto(\lambda= {}, \alpha= {}), k= {}$'.format(loc, a, K) )
  
  plot.xlabel(r'$E[T]$', fontsize=12)
  plot.ylabel(r'$E[C]$', fontsize=12)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_k_nc_retainl_atd__C_vs_T_{}_k_{}.png".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

if __name__ == "__main__":
  # plot_arepeat_Pr_T_g_t()
  
  # plot_arepeat_k_l_n__E_C_vs_E_T()
  # plot_arepeat_E_T_k_n_vs_k_l_n()
  # plot_arepeat_E_T_k_l_n()
  
  # plot_dE_T_shiftedexp_k_n_dk()
  # plot_arepeat_shiftedexp_k_n()
  
  # plot_EC_vs_ET_wdelay(w_cancel=True)
  # plot_zerodelay_reped_vs_coded(loc=3, a=2)
  # plot_zerodelay_reped_vs_coded(loc=3, a=1.2)
  # plot_zerodelay_reped_vs_coded(loc=3, a=1.5)
  
  # plot_k_nc()
  # plot_k_nc_E_C_vs_E_T()
  
  # plot_k_nc_wrelaunch()
  # plot_k_ncd0_wrelaunch()
  plot_k_nc_wrelaunch_EC_vs_ET()
  # plot_k_ncd0_wrelaunch_EC_vs_ET()
  
  # plot_arepeat_Pr_T_g_t_pareto_k_wrelaunch()
  
  # plot_reduct_in_ET_atnocost_wred()
  
  # plot_k_nc_retainl_atd()
  # plot_k_nc_retainl_atd__C_vs_T()
  
  # plot_cost_reduction_vs_redundant_occupation()
 
