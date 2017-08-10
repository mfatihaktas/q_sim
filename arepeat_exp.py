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

def plot_arepeat_k_nc():
  K = 10
  N = 14
  D, mu = 30, 1
  loc, a = 3, 2
  task_t = "Pareto" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim SExp(D/k, \mu), D={}, \mu={}'.format(D, mu)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  
  def plot_(k, n=None, c=None, sim=False):
    if task_t == "Exp": task_t_rv = Exp(mu)
    elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
    elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    
    d_l = []
    y_l, y_approx_l, y_sim_l = [], [], []
    u_l = 3*H(K)/mu + 1
    for d in numpy.arange(0, u_l, 0.5):
      d_l.append(d)
      if task_t != "Google":
        if n is not None:
          y_l.append(E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n) )
          # y_l.append(E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n, w_cancel=True) )
        elif c is not None:
          y_l.append(E_T_k_c(task_t, D, mu, loc, a, d, k, c) )
          # y_l.append(E_C_k_l_n(task_t, D, mu, loc, a, d, k, c, w_cancel=True) )
      if sim:
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=10000)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=10000)
        # y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
        y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] ) )
        # y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C'] )/len(stat_id__trial_sampleavg_l_m['E_C'] ) )
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    # plot.plot(d_l, y_l, label=label, color=next(dark_color), linestyle=':', mew=2)
    # plot.plot(d_l, y_approx_l, label=label, color=next(light_color), marker=next(marker), linestyle='', mew=2)
    plot.plot(d_l, y_sim_l, label='Sim, {}'.format(label), color=next(light_color), marker=next(marker), linestyle='', mew=2)
  
  sim = False
  k_step = 2
  for n in [K, K+1, *range(K+k_step, N, k_step)]:
    plot_(K, n=n, sim=sim)
  plot_(K, n=20, sim=sim)
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  plot.xlabel(r'$\Delta$ (s)')
  # plot.ylabel(r'Expected latency $E[T]$ (s)')
  plot.ylabel(r'Expected cost $E[C^c]$ (s)')
  plot.savefig("plot_arepeat_k_nc_{}_k_{}.png".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )
  
def plot_arepeat_k_nc_E_C_vs_E_T():
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
        E_T_l.append(E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n) )
        E_C_l.append(E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n, w_cancel=w_cancel) )
      if sim:
        if task_t == "Exp": task_t_rv = Exp(mu)
        elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
        elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
        elif task_t == "Google": task_t_rv = Google(k)
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=10000)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=10000)
        E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )
        E_C = sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] )/60/60
        if task_t == "Google":
          E_T = E_T/60/60
          E_C = E_C/60/60
        E_T_sim_l.append(E_T)
        E_C_sim_l.append(E_C)
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    # plot.plot(E_T_l, E_C_l, label=label, color=next(dark_color), linestyle=':', mew=2)
    plot.plot(E_T_sim_l, E_C_sim_l, label='Sim, {}'.format(label), marker=next(marker), color=next(dark_color), linestyle=':', mew=2)
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
  plot.savefig("plot_arepeat_k_nc_E_C_vs_E_T_{}_k_{}.png".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

# ***************************************  REP vs. CODING  *************************************** #
def plot_reped_vs_coded(w_cancel=True):
  K = 10
  D = 30
  mu = 0.5
  loc, a = 3, 2
  task_t = "SExp" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim D/k + Exp(\mu), D={}, k={}, \mu={}'.format(D, K, mu)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  ann_color = 'black'
  def plot_d_extremes(c=0, n=0):
    l, k, d = K, K, 0
    E_T_l, E_C_l = [], []
    if c:
      for c_ in range(1,c+1):
        E_T = E_T_k_c(task_t, D, mu, loc, a, d, k, c_)
        E_C = E_C_k_c(task_t, D, mu, loc, a, d, k, c_, w_cancel=w_cancel)
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        if c_ == c:
          plot.annotate(r'$\Delta=0$', ha='center', va='center', xy=(E_T, E_C), xytext=(E_T, E_C-10), color=ann_color, fontsize=20)
      plot.plot(E_T_l, E_C_l, color=ann_color, alpha=0.6, linestyle='--')
    if n:
      E_T_l.clear()
      E_C_l.clear()
      for n_ in range(k+1, n+1):
        E_T = E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n_)
        E_C = E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n_, w_cancel=w_cancel)
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        if n_ == n:
          plot.annotate(r'$\Delta=0$', ha='center', va='center', xy=(E_T, E_C), xytext=(E_T, E_C-5), color=ann_color, fontsize=20)
      plot.plot(E_T_l, E_C_l, color=ann_color, alpha=0.6, linestyle='--')
  
  def plot_(k=K, n=0, c=0):
    l = k
    if task_t == "Exp": task_t_rv = Exp(mu)
    elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
    elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    
    E_T_l, E_T_sim_l, E_C_l, E_C_sim_l = [], [], [], []
    num_run = 100000
    u_l = 25 + 1
    color = next(dark_color)
    if c:
      for d in numpy.arange(0, u_l, 0.5):
        # stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=num_run)
        # E_T_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
        # key = 'E_C_wc' if w_cancel else 'E_C'
        # E_C_sim_l.append(sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] ) )
        
        E_T = E_T_k_c(task_t, D, mu, loc, a, d, k, c)
        E_C = E_C_k_c(task_t, D, mu, loc, a, d, k, c, w_cancel=w_cancel)
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        if d == 0:
          plot.annotate(r'$c={}$'.format(c), ha='center', va='center', xy=(E_T, E_C), xytext=(E_T-0.3, E_C), color=color, fontsize=12)
      plot.plot(E_T_l, E_C_l, color=color, marker=next(marker), ms=10, mew=2, zorder=2, lw=2, linestyle='-')
      # plot.plot(E_T_l, E_C_l, label=r'Rep,$c={}$'.format(c), color=next(dark_color), marker=next(marker), ms=10, mew=2, zorder=2, lw=2, linestyle='-')
      # plot.plot(E_T_sim_l, E_C_sim_l, label=r'Rep,$c={}$'.format(c), color=next(dark_color), marker=next(marker), zorder=2, mew=5, linestyle=':')
    elif n:
      for d in numpy.arange(0, u_l, 0.5):
        # stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n, num_run=num_run)
        # E_T_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
        # key = 'E_C_wc' if w_cancel else 'E_C'
        # E_C_sim_l.append(sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] ) )
        E_T = E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n)
        E_C = E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n, w_cancel=w_cancel)
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        if d == 0:
          plot.annotate(r'$n={}$'.format(n), ha='center', va='center', xy=(E_T, E_C), xytext=(E_T-0.4, E_C), color=color, fontsize=12)
      # label=r'Code,$n={}$'.format(n)
      plot.plot(E_T_l, E_C_l, color=color, marker=next(marker), ms=9, mew=2, zorder=1, lw=2, linestyle=':')
      # plot.plot(E_T_sim_l, E_C_sim_l, label=r'Code,$n={}$'.format(n), color=next(dark_color), marker=next(marker), ms=8, zorder=1, mew=2, linestyle=':')
  plot_(c=1)
  plot_(c=2)
  plot_(c=3)
  plot_d_extremes(c=3)
  
  # plot_(n=K+1)
  # plot_(n=K+2)
  # for n in range(K+5, 3*K+1, 5):
  #   plot_(n=n)
  # plot_d_extremes(n=3*K)
  #
  if task_t == "SExp":
    E_T_nored = E_T_k_c(task_t, D, mu, loc, a, d=0, k=K, c=0)
    E_C_nored = E_C_k_c(task_t, D, mu, loc, a, d=0, k=K, c=0, w_cancel=w_cancel)
    plot.plot([E_T_nored], [E_C_nored], 'o', zorder=3, mew=4, color=ann_color)
    plot.annotate(r'$\Delta \to \infty$', xy=(E_T_nored, E_C_nored), ha='center', va='center', xytext=(E_T_nored+0.5, E_C_nored), color=ann_color, fontsize=20)
    axes = plot.gca()
    # axes.set_xlim([3.2, 9.2] )
    axes.set_xlim([3.9, 9.2] ) # for only rep
    # axes.set_ylim([47, 115] )
  # legend = plot.legend()
  # # legend = plot.legend(loc='center left', bbox_to_anchor=(0.85, 0.7), fontsize=10)
  # frame = legend.get_frame()
  # frame.frameon=True
  # frame.set_edgecolor('black')
  plot.title(r'${}$'.format(task_t_in_latex), fontsize=12)
  plot.xlabel(r'Expected Latency $E[T]$ (s)', fontsize=12)
  # w_cancel_text = "w/ cancel $E[C^c]$" if w_cancel else "w/o cancel $E[C]$"
  w_cancel_text = r'$E[C]$'
  plot.ylabel(r'Expected Cost {} (s)'.format(w_cancel_text), fontsize=12)
  fig = plot.gcf()
  # plot.savefig("plot_reped_vs_coded_{}_k_{}.pdf".format(task_t, K, w_cancel), bbox_extra_artists=(legend,), bbox_inches='tight', dpi=fig.dpi)
  plot.savefig("plot_reped_vs_coded_{}_k_{}.pdf".format(task_t, K, w_cancel), bbox_inches='tight', dpi=fig.dpi)
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_reped_vs_coded_zerodelay(loc, a):
  w_cancel=True
  K = 15 # 400 # 1050 # 15 # 400 # 10
  D, mu = 30, 0.5
  loc, a = 3, 1.5 # 2 # 3
  num_run = 10000
  task_t = "Pareto" # "Google" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim SExp(D/k, \mu={}), D={}'.format(mu, D)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  elif task_t == "Google": task_t_in_latex = r'X \sim Google'
  mew = 3
  def plot_(k=K, n=0, c=0, sim=False):
    l = k
    if task_t == "Exp": task_t_rv = Exp(mu)
    elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
    elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    elif task_t == "Google": task_t_rv = Google(k)
    
    E_T_l, E_T_sim_l, E_C_l, E_C_sim_l = [], [], [], []
    d = 0
    color = next(dark_color)
    if c:
      for c_ in range(c+1):
        # E_T = E_T_k_c(task_t, D, mu, loc, a, d, k, c_)
        # E_C = E_C_k_c(task_t, D, mu, loc, a, d, k, c_, w_cancel=w_cancel)
        E_T = E_T_2_pareto_k_c(loc, a, k, c_)
        E_C = E_C_2_pareto_k_c(loc, a, k, c_)
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        if c_ > 0:
          if task_t == "SExp":
            plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C), color=color)
          elif task_t == "Pareto":
            if a == 1.2:
              plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+1, E_C+1), color=color)
            elif a == 1.5:
              plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.3, E_C+0.3), color=color)
            elif a == 2:
              plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.3, E_C+0.3), color=color)
        if sim:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c_, num_run)
          E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )/60/60
          key = 'E_C_wc' if w_cancel else 'E_C'
          E_C = sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] )/60/60
          E_T_sim_l.append(E_T)
          E_C_sim_l.append(E_C)
          if task_t == "Google":
            if k == 400:
              if c_ == 0:
                plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T, E_C), xytext=(E_T-0.2, E_C+2.5) )
              else:
                plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.025, E_C+1), color=color)
            elif k == 1050:
              if c_ == 0:
                plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T, E_C), xytext=(E_T-0.2, E_C+8) )
              else:
                plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.025, E_C+1), color=color)
            elif k == 15:
              if c_ == 0:
                plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T, E_C), xytext=(E_T-0.0275, E_C+0.05) )
              else:
                plot.annotate(r'$c={}$'.format(c_), xy=(E_T+0.003, E_C), xytext=(E_T, E_C), color=color)
      plot.plot(E_T_l, E_C_l, label='Replication', color=color, marker=next(marker), zorder=0, mew=2, linestyle=':')
      # plot.plot(E_T_sim_l, E_C_sim_l, label='Simulation, replication', color=color, marker=next(marker), linestyle=':', mew=2)
    elif n:
      counter = 0
      # for n_ in [*numpy.arange(k, k+50, 2), *numpy.arange(k+50, n+1, 100) ]:
      for n_ in numpy.arange(k, n+1, 1):
        # E_T = E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n_)
        # E_C = E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n_, w_cancel=w_cancel)
        E_T = E_T_2_pareto_k_n(loc, a, k, n_)
        E_C = E_C_2_pareto_k_n(loc, a, k, n_)
        E_T_l.append(E_T)
        E_C_l.append(E_C)
        
        if sim:
          counter += 1
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n_, num_run)
          E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )/60/60
          key = 'E_C_wc' if w_cancel else 'E_C'
          E_C = sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] )/60/60
          E_T_sim_l.append(E_T)
          E_C_sim_l.append(E_C)
          if task_t == "Google":
            if counter % 10 == 0:
              plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
              if k == 400:
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.02, E_C+1), color=color)
              elif k == 1050:
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.02, E_C+1), color=color)
              elif k == 15:
                plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.003, E_C), color=color)
        
        # if task_t == "SExp":
        #   if n_ > K and n_ < k+3:
        #     plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C), color=color)
        #   elif n_ != k and n_ % k == 0:
        #     plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #     plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C), color=color)
        # elif task_t == "Pareto":
        #   if a == 1.2:
        #     if n_ > K and n_ < k+3:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+2, E_C-2), color=color)
        #       # elif n_ == k+3:
        #       #   plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       #   plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-4, E_C-6), color=color)
        #     elif n_ == 30:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-15, E_C-8), color=color)
        #     elif n_ != k and n_ % k == 0:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-15, E_C-3), color=color)
        #   elif a == 1.5:
        #     if n_ > K and n_ < k+3:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T, E_C-5), color=color)
        #       # elif n_ == k+3:
        #       #   plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       #   plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-3.5, E_C-5), color=color)
        #     elif n_ == 40 or n_ == 60:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-4, E_C+1), color=color)
        #     elif n_ != k and n_ % k == 0:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-4, E_C-2), color=color)
        #   elif a == 2:
        #     if n_ > K and n_ < k+3:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T, E_C-5), color=color)
        #       # elif n_ == k+3:
        #       #   plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       #   plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-1.5, E_C-5), color=color)
        #     elif n_ != k and n_ % k == 0:
        #       plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
        #       plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-2, E_C-2), color=color)
      plot.plot(E_T_l, E_C_l, label='Coding', color=color, zorder=1, marker=next(marker), mew=1, linestyle=':')
      # plot.plot(E_T_sim_l, E_C_sim_l, label='Simulation, coding', color=color, zorder=1, marker=next(marker), linestyle=':', mew=1)
  plot_(c=5)
  plot_(n=6*K)
  #
  # E_T_nored = E_T_k_c(task_t, D, mu, loc, a, d=0, k=K, c=0)
  # E_C_nored = E_C_k_c(task_t, D, mu, loc, a, d=0, k=K, c=0, w_cancel=w_cancel)
  # if task_t == "SExp":
  #   plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-0.6, E_C_nored+9) )
  # elif task_t == "Pareto":
  #   plot.axhline(y=E_C_nored, color='k', alpha=0.4, linestyle='--')
  #   if a == 1.2:
  #     plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-20, E_C_nored+3) )
  #   elif a == 1.5:
  #     plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-6, E_C_nored+3) )
  #   elif a == 2:
  #     plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-1.2, E_C_nored+3) )
  plot.legend()
  # plot.legend(loc='lower right')
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  # plot.xlabel(r'Expected Latency $E[T]$ (s)', fontsize=12)
  # plot.xlabel(r'Expected latency $E[T]$ (hour)', fontsize=12)
  plot.xlabel(r'Expected Latency $E[T^2]$', fontsize=12)
  # plot.ylabel(r'Expected cost $E[C^c]$' if w_cancel else r'$E[C]$', fontsize=12)
  plot.ylabel(r'Expected cost $E[C^2]$', fontsize=12)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  if task_t == "Pareto":
    plot.savefig("plot_reped_vs_coded_zerodelay_{}_k_{}_a_{}.png".format(task_t, K, a) )
  else:
    plot.savefig("plot_reped_vs_coded_zerodelay_{}_k_{}.pdf".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_reduct_in_E_T_atnocost_w_zerodelay_red():
  loc = 3
  d = 0
  def reduction_in_E_T_for_same_E_C_w_red(red_type, loc, a, k):
    E_T_wo_red, E_T = 0, 0
    if red_type == 'coded':
      E_T_wo_red = E_T_pareto_k_n(loc, a, d, k, n=k)
      E_C_wo_red = E_C_pareto_k_n(loc, a, d, k, n=k)
      
      n_ = k+1
      while 1:
        E_C = E_C_pareto_k_n(loc, a, d, k, n_)
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
    #             color=color, alpha=0.5, linestyle='--')
    legend = "Replication" if red_type == 'reped' else "Coding"
    plot.plot(x_l, y_l, label=r'{},$k={}$'.format(legend, k), color=next(dark_color), marker=next(marker), mew=2, zorder=0, linestyle=':')
    # plot.plot(x_l, y_approx_l, label=r'{},$k={}$, approx'.format(legend, k), color=next(dark_color), marker=next(marker), mew=2, zorder=1, linestyle='')
  
  plot_('reped', k=10)
  # plot_('coded', k=10)
  
  # plot_('reped', k=50)
  # plot_('coded', k=50)
  
  plot.legend()
  plot.title(r'$X \sim Pareto(\lambda={}, \alpha)$'.format(loc) )
  plot.xlabel(r'Tail index $\alpha$', fontsize=12)
  plot.ylabel('\n'.join(textwrap.wrap(r'Maximum percetange reduction in $E[T]$ at no added cost', 40) ),
              fontsize=12)
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_reduct_in_E_T_atnocost_w_zerodelay_red.pdf", bbox_inches='tight')
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
    plot.plot(t_l, y_l, label=label, color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
    # plot.plot(t_l, y_approx_l, label=r'Approx, $\Delta={}$'.format(d), color=next(light_color), marker=next(marker), linestyle=':', mew=2)
  
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
  #   plot.plot(d_l, y_l, label=r'$t={0:.2f}$'.format(t), color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  
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
  plot.savefig("plot_arepeat_Pr_T_g_t_pareto_k_wrelaunch_k_{}.pdf".format(K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_arepeat_k_nc_wrelaunch():
  K = 100
  D, mu = 30, 1
  loc, a = 3, 2 # 100 # 1.2 # 2
  task_t = "Pareto"
  
  def plot_reduction_in_E_T_wrelaunch(k, loc=loc):
    def max_a_lb():
      return math.log(k+1)/math.log(4)
    
    x_l, y_l = [], []
    for a in numpy.linspace(1.05, 8, 50):
      x_l.append(a)
      
      E_T_wrelaunch_min = float('Inf')
      for d in numpy.linspace(loc+0.05, 10*loc, 1000):
        E_T = E_T_pareto_k_n(loc, a, d, k, n=k, w_relaunch=True)
        if E_T < E_T_wrelaunch_min:
          E_T_wrelaunch_min = E_T
        else:
          # elif not math.isinf(E_T_wrelaunch_min):
          break
      E_T_norelaunch = E_T_pareto_k_n(loc, a, d, k, n=k, w_relaunch=False)
      y_l.append(max((E_T_norelaunch - E_T_wrelaunch_min)/E_T_norelaunch, 0) )
    c = next(dark_color)
    label = r'$k={}$'.format(k) # r'$\lambda={}$, $k={}$'.format(loc, k)
    plot.plot(x_l, y_l, label=label, color=c, marker=next(marker), linestyle=':', mew=2)
    plot.axvline(x=max_a_lb(), label=r'$\alpha_u$, {}'.format(label), color=c, linestyle='--')
  
  E_T = True # False # True
  w_cancel = True # False
  def plot_(k, n=None, c=None, sim=False):
    task_t_rv = Pareto(a, loc=loc)
    
    l = k
    x_l = []
    y_wrelaunch_l, y_wrelaunch_approx_l, y_sim_wrelaunch_l = [], [], []
    y_norelaunch_l = []
    for d in numpy.linspace(0, 20*loc, 50):
      x_l.append(d)
      
      if E_T:
        if n is not None:
          y_wrelaunch_l.append(E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=True) )
          # y_wrelaunch_approx_l.append(E_T_pareto_k_n_approx(loc, a, d, k, n, w_relaunch=True) )
          if n == k: y_norelaunch_l.append(E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=False) )
        elif c is not None:
          y_wrelaunch_l.append(E_T_pareto_k_c_wrelaunch(loc, a, d, k, c, w_relaunch=True) )
          # y_wrelaunch_approx_l.append(E_T_pareto_k_c_wrelaunch_approx(loc, a, d, k, c, w_relaunch=True) )
          if c == 0: y_norelaunch_l.append(E_T_pareto_k_c_wrelaunch(loc, a, d, k, c, w_relaunch=False) )
      else:
        if n is not None:
          y_wrelaunch_l.append(E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel) )
          # y_wrelaunch_approx_l.append(E_C_pareto_k_n_approx(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
        elif c is not None:
          y_wrelaunch_l.append(E_C_pareto_k_c_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel) )
          # y_wrelaunch_approx_l.append(E_C_pareto_k_c_approx(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
      
      if sim:
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=10000, w_relaunch=True)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=10000, w_relaunch=True)
        if E_T:
          y_sim_wrelaunch_l.append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
        else:
          if w_cancel:
            y_sim_wrelaunch_l.append(sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] ) )
          else:
            y_sim_wrelaunch_l.append(sum(stat_id__trial_sampleavg_l_m['E_C'] )/len(stat_id__trial_sampleavg_l_m['E_C'] ) )
    c = next(dark_color)
    if n == k:
      d = Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, k)
      # plot.axvline(x=d, label=r'$\Delta^*$, n={}'.format(n), color=c)
      y = E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=True)
      plot.plot(d, y, color='red', label=r'Approx optimal', marker='*', zorder=1, ms=10, mew=3)
    # print("y_wrelaunch_l= {}".format(y_wrelaunch_l) )
    # label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    label = 'With relaunch'
    plot.plot(x_l, y_wrelaunch_l, label=label, color=c, marker=next(marker), linestyle='--', zorder=0, mew=2)
    # plot.plot(x_l, y_wrelaunch_approx_l, label='approx, {}'.format(label), color=next(dark_color), marker=next(marker), linestyle='', mew=2)
    
    if sim:
      plot.plot(x_l, y_sim_wrelaunch_l, label=r'sim, {}'.format(label), color=next(light_color), marker=next(marker), linestyle='', mew=2)
    if E_T:
      if n == k: plot.plot(x_l, y_norelaunch_l, label='No relaunch', color=next(dark_color), linestyle='--', mew=2)
  
  sim = False # True
  
  # plot_reduction_in_E_T_wrelaunch(K)
  # plot_reduction_in_E_T_wrelaunch(10*K)
  # plot_reduction_in_E_T_wrelaunch(100*K)
  # plot_reduction_in_E_T_wrelaunch(1000*K)
  
  # plot_reduction_in_E_T_wrelaunch(1, 10*K)
  # plot_reduction_in_E_T_wrelaunch(loc, 10*K)
  # plot_reduction_in_E_T_wrelaunch(10*loc, 10*K)
  # plot_reduction_in_E_T_wrelaunch(100*loc, 10*K)
  
  plot_(K, n=K)
  # plot_(2*K, n=2*K)
  # plot_(10*K, n=10*K)
  # plot_(20*K, n=20*K)
  
  # plot_(K, n=K)
  # plot_(K, n=K+1, sim=sim)
  # plot_(K, n=K+2, sim=sim)
  # plot_(K, n=K+3, sim=sim)
  
  # plot_(K, c=0, sim=sim)
  # plot_(K, c=1, sim=sim)
  # plot_(K, c=2, sim=sim)
  # plot_(K, c=3, sim=sim)
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.7, 0.5) )
  plot.title(r'$X \sim Pareto(\lambda= {}, \alpha= {}), k= {}$'.format(loc, a, K) )
  # plot.title(r'$X \sim Pareto(\lambda= {}, \alpha)$'.format(loc) )
  
  plot.xlabel(r'Relaunch delay $\Delta$ (s)', fontsize=13)
  y_label = r'Expected latency $E[T]$' if E_T else r'Expected cost $E[C]$'
  plot.ylabel(r'{} (s)'.format(y_label), fontsize=13)
  # plot.xlabel(r'$\alpha$', fontsize=12)
  # plot.ylabel('\n'.join(textwrap.wrap(r'Maximum percetange reduction in $E[T]$ by task relaunch', 40) ), fontsize=12)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_arepeat_k_nc_wrelaunch_{}_k_{}.pdf".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_arepeat_k_nc_wrelaunch_E_C_vs_E_T():
  K = 15 # 1050 # 100
  loc, a = 3, 2
  w_cancel = True
  
  task_t = "Pareto" # "Google"
  if task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  elif task_t == "Google": task_t_in_latex = r'X \sim Google'
  
  def plot_(k, n=None, c=None, sim=False):
    num_f = 5
    num_run = 10000 # 10000
    E_T_l, E_C_l = [], []
    E_T_sim_l, E_C_sim_l = [], []
    
    if task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    elif task_t == "Google": task_t_rv = Google(k)
    
    # for d in [*numpy.linspace(0, 5*loc, 50), *numpy.linspace(5*loc, 200*loc, 100)]:
    # for d in numpy.logspace(-3, 2, 100):
    d_ul = 10**5
    for d in numpy.logspace(-2, 4, 40):
      if sim:
        E_T_sum, E_C_sum = 0, 0
        for i in range(num_f):
          if n is not None:
            stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n, num_run, w_relaunch=True)
          elif c is not None:
            stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run, w_relaunch=True)
          E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )
          E_C = sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] )
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
          E_C_l.append(E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel) )
        elif c is not None:
          E_T_l.append(E_T_pareto_k_c_wrelaunch(loc, a, d, k, c) )
          E_C_l.append(E_C_pareto_k_c_wrelaunch(loc, a, d, k, c, w_cancel=w_cancel) )
    color = next(dark_color)
    if n != k:
      label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    else:
      label = "No redundancy with relaunch"
    plot.plot(E_T_l, E_C_l, label=label, color=color, marker=next(marker), zorder=0, linestyle=':', mew=2)
    # plot.plot(E_T_sim_l, E_C_sim_l, label='Simulation, {}'.format(label), color=next(dark_color), marker=next(marker), zorder=0, linestyle=':', mew=2)
    
    if n is not None and n == k:
      if task_t != "Google":
        x = E_T_pareto_k_n_wrelaunch(loc, a, 0, k, n)
        y = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=w_cancel)
        if k == 100 and a == 2:
          # plot.annotate(r'$\Delta=0$ or $\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=12)
          plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+1, y+4), fontsize=12)
          plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x-1.5, y-25), fontsize=12)
          
          d = Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, k)
          x = E_T_pareto_k_n_wrelaunch(loc, a, d, k, n)
          y = E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel)
          plot.plot(x, y, color='red', label=r'Approx optimal', marker='*', zorder=1, ms=10, mew=3)
          plot.annotate(r'$\Delta={0:.2f}$'.format(d), xy=(x, y), xytext=(x-1.5, y+20), fontsize=12)
        else:
          # plot.plot(x, y, color='black', marker='*', ms=8, mew=3)
          plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x, y), fontsize=12)
      else: # task_t = "Google"
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, 0, k, k, k, num_run, w_relaunch=True)
        x = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )/60/60
        y = sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] )/60/60
        plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+0.005, y+2.5), fontsize=12)
        
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d_ul, k, k, k, num_run, w_relaunch=True)
        x = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )/60/60
        y = sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] )/60/60
        plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x-0.005, y-2.5), xytext=(x, y), fontsize=12)
    elif n is not None and n > k:
      x = E_T_pareto_k_n_wrelaunch(loc, a, 0, k, n)
      y = E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=w_cancel)
      if k == 100 and a == 2:
        plot.plot(x, y, color='black', marker='o', zorder=1, ms=8, mew=2)
        # plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+1, y+5), color=color, fontsize=13)
      x = E_T_pareto_k_n_wrelaunch(loc, a, d_ul, k, n)
      y = E_C_pareto_k_n_wrelaunch(loc, a, d_ul, k, n, w_cancel=w_cancel)
      if n == 101 and k == 100 and a == 2:
        plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=13)
    elif c is not None:
      x = E_T_pareto_k_c_wrelaunch(loc, a, 0, k, c)
      y = E_C_pareto_k_c_wrelaunch(loc, a, 0, k, c, w_cancel=w_cancel)
      plot.plot(x, y, color='black', marker='o', zorder=1, ms=6, mew=2)
  """
  # plot_(K, n=K)
  
  # plot_(K, n=K+2)
  # plot_(K, n=K+1)
  
  plot_(K, n=K+10)
  plot_(K, n=K+20)
  plot_(K, n=K+30)
  k, n = K, K+30
  x = E_T_pareto_k_n(loc, a, 0, k, n, w_relaunch=True)
  y = E_C_pareto_k_n(loc, a, 0, k, n, w_cancel=w_cancel, w_relaunch=True)
  plot.plot(x, y, color='black', label=r'$\Delta=0$', marker='o', zorder=1, ms=8, mew=2)
  x = E_T_pareto_k_n(loc, a, 10**5, k, n, w_relaunch=True)
  y = E_C_pareto_k_n(loc, a, 10**5, k, n, w_cancel=w_cancel, w_relaunch=True)
  plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=13)
  plot_(K, c=1)
  plot_(K, c=2)
  plot_(K, c=3)
  k, c = K, 3
  x = E_T_pareto_k_c_wrelaunch(loc, a, 0, k, c, w_relaunch=True)
  y = E_C_pareto_k_c_wrelaunch(loc, a, 0, k, c, w_cancel=w_cancel, w_relaunch=True)
  plot.plot(x, y, color='black', label=r'$\Delta=0$', marker='o', zorder=1, ms=6, mew=2)
  x = E_T_pareto_k_c_wrelaunch(loc, a, 10**2, k, c, w_relaunch=True)
  y = E_C_pareto_k_c_wrelaunch(loc, a, 10**2, k, c, w_cancel=w_cancel, w_relaunch=True)
  plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=13)
  """
  plot_(K, n=K, sim=False)
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  unit = '(s)' if task_t != "Google" else '(hours)'
  plot.xlabel(r'Expected Latency $E[T]$ {}'.format(unit), fontsize=13)
  y = r'$E[C^c]$' if w_cancel else r'$E[C]$'
  plot.ylabel('Expected Cost {} {}'.format(y, unit), fontsize=13)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_arepeat_k_nc_wrelaunch_E_C_vs_E_T_{}_k_{}.pdf".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

if __name__ == "__main__":
  # plot_arepeat_Pr_T_g_t()
  
  # plot_arepeat_k_l_n__E_C_vs_E_T()
  # plot_arepeat_E_T_k_n_vs_k_l_n()
  # plot_arepeat_E_T_k_l_n()
  
  # plot_dE_T_shiftedexp_k_n_dk()
  # plot_arepeat_shiftedexp_k_n()
  
  # plot_reped_vs_coded(w_cancel=True)
  plot_reped_vs_coded_zerodelay(loc=3, a=2)
  # plot_reped_vs_coded_zerodelay(loc=3, a=1.2)
  # plot_reped_vs_coded_zerodelay(loc=3, a=1.5)
  
  # plot_arepeat_k_nc()
  # plot_arepeat_k_nc_E_C_vs_E_T()
  # plot_arepeat_k_nc_wrelaunch()
  # plot_arepeat_k_nc_wrelaunch_E_C_vs_E_T()
  # plot_arepeat_Pr_T_g_t_pareto_k_wrelaunch()
  
  # plot_arepeat_k_c()
  # plot_arepeat_k_c_E_C_vs_E_T()
  
