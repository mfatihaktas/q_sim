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

# ##################  Send n any k is enough, X_i ~ Exp(mu), each packet drops ~ Exp(gamma)  ################ #
def plot_send_n_w_drop():
  mu = 1
  gamma = 0.2
  k = 10
  N = 20
  
  plot_Pr_succ = False # True
  
  n__y_l_m, n__y_approx_l_m = {}, {}
  n__x_l_m = {}
  
  k_step = 2
  for n in range(k, N, k_step):
    n__y_l_m[n] = []
    n__y_approx_l_m[n] = []
    n__x_l_m[n] = []
    for gamma_ in numpy.arange(gamma, 2.9, 0.1):
      n__x_l_m[n].append(gamma_)
      if plot_Pr_succ:
        n__y_l_m[n].append(Pr_succ_n_k_w_drop(mu, gamma_, n, k) )
      else:
        n__y_l_m[n].append(E_T_n_k_w_drop_given_succ(mu, gamma_, n, k) )
        n__y_approx_l_m[n].append(E_T_n_k_w_drop_given_succ_approx(mu, gamma_, n, k) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  color = iter(cm.rainbow(numpy.linspace(0, 1, int((N-k)/k_step)+1) ) )
  for n in range(k, N, k_step):
    m = next(marker)
    c = next(color)
    plot.plot(n__x_l_m[n], n__y_l_m[n], label=r'n:{}'.format(n), color=c, marker=m, linestyle='', mew=2)
    plot.plot(n__x_l_m[n], n__y_approx_l_m[n], label=r'~ n:{}'.format(n), color=c, alpha=0.7, marker=m, linestyle='', mew=2)
  # plot.legend()
  plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.xlabel(r'$\gamma$')
  plot.title(r'$\mu$= {}, k= {}'.format(mu, k) )
  if plot_Pr_succ:
    plot.ylabel(r'$Pr\{succ\}$')
  else:
    plot.ylabel(r'$E[T|succ]$ (s)')
  plot.savefig("plot_send_n_w_drop_k_{}.png".format(k) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(k) )

# ##################  X_i ~ G, (l=k, k, n=k+1, \Delta)  ################ #
def plot_arepeat_E_T_vs_E_C_Gred():
  k = 7
  red = 1
  def plot_(k, rv, rv_in_latex):
    x_l_m, y_l_m, y2_l_m = [], [], []
    x_sim_l_m, y_sim_l_m, y2_sim_l_m = [], [], []
    
    for d in numpy.arange(0, 5*rv.mean(), 1):
      if red == 1:
        x_l_m.append(d)
        # x_l_m.append(E_T_G_1red(rv, d, k) )
        y_l_m.append(E_C_G_1red(rv, d, k, w_cancel=True) )
        y2_l_m.append(E_C_G_1red(rv, d, k, w_cancel=False) )
      
      stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(rv, d, k, l=k, n=k+red, num_run=100000)
      E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )
      E_C_wc = sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] )
      E_C = sum(stat_id__trial_sampleavg_l_m['E_C'] )/len(stat_id__trial_sampleavg_l_m['E_C'] )
      x_sim_l_m.append(d)
      # x_sim_l_m.append(E_T)
      y_sim_l_m.append(E_C_wc)
      y2_sim_l_m.append(E_C)
    if red == 1:
      plot.plot(x_l_m, y_l_m, label=r'w/ cancel,{}'.format(rv_in_latex), color=next(dark_color), marker=next(marker), linestyle='-', mew=2)
      plot.plot(x_l_m, y2_l_m, label=r'no cancel,{}'.format(rv_in_latex), color=next(dark_color), marker=next(marker), linestyle='-', mew=2)
    plot.plot(x_sim_l_m, y_sim_l_m, label=r'w/ cancel,{}'.format(rv_in_latex), color=next(light_color), alpha=1, marker=next(marker), linestyle='', mew=2)
    plot.plot(x_sim_l_m, y2_sim_l_m, label=r'no cancel,{}'.format(rv_in_latex), color=next(light_color), alpha=1, marker=next(marker), linestyle='', mew=2)
  
  # plot_(7, Pareto(a=3, loc=6), r"$Pareto(\alpha=3,\lambda=6)$")
  # task_t_rv, name = Pareto(a=3, loc=6), "$Pareto$"
  # task_t_rv, task_t_in_latex = Pareto(a=3, loc=6), r"$Pareto(\alpha=3,\lambda=6)$"
  # task_t_rv, task_t_in_latex = Exp(mu=1/3, D=6), r"$ShiftedExp(D=6,\mu=1/3)$"
  # task_t_rv, task_t_in_latex = Exp(mu=1/9), r"$Exp(\mu=1/3)$"
  # plot_(7, task_t_rv, task_t_in_latex)
  
  # plot_(k, Pareto(a=3, loc=6), r"$Pareto(\alpha=3,\lambda=6)$")
  # plot_(k, Exp(mu=1/3, D=6), r"$ShiftedExp(D=6,\mu=1/3)$")
  # plot_(k, Exp(mu=1/9), r"$Exp(\mu=1/3)$")
  plot_(k, Pareto(a=3, loc=6), r"$Pareto$")
  plot_(k, Exp(mu=1/3, D=6), r"$SExp$")
  plot_(k, Exp(mu=1/9), r"$Exp$")
  
  # plot.legend()
  plot.legend(loc='center left', bbox_to_anchor=(0.75, 0.5) )
  plot.xlabel(r'$\Delta$ (s)')
  # plot.xlabel(r'$E[T]$ (s)')
  plot.ylabel(r'$E[C]$ (s)')
  plot.title(r'$k= {}, l= {}, n= {}$'.format(k, k, k+red) )
  # plot.title(r'$X_i \sim $ {}'.format(task_t_in_latex) )
  # plot.savefig("plot_arepeat_E_T_vs_E_C_Gred_{}.png".format(task_t_rv) )
  plot.savefig("plot_arepeat_E_T_vs_E_C_G_red_{}_k_{}.png".format(red, k) )
  plot.gcf().clear()
  log(WARNING, "done.")
  
def plot_arepeat_E_T_G_1red(w_cancel=True):
  K = 10
  # task_t_rv = Exp(mu=2/3)
  task_t_rv = Pareto(a=3, loc=6)
  
  k__x_l_m, k__y_l_m, k__y_approx_l_m = {}, {}, {}
  k__x_sim_l_m, k__y_sim_l_m = {}, {}
  k_step = 2
  for k in range(2, K, k_step):
  # for k in range(K, K+1, k_step):
    k__x_l_m[k] = []
    k__y_l_m[k] = []
    # k__y_approx_l_m[k] = []
    k__x_sim_l_m[k] = []
    k__y_sim_l_m[k] = []
    # for d in numpy.arange(0, K*task_t_rv.mean(), 0.1):
    for d in numpy.arange(0, 2*task_t_rv.mean(), 0.1):
      k__x_l_m[k].append(d)
      k__y_l_m[k].append(E_T_G_1red(task_t_rv, d, k) )
      # k__y_approx_l_m[k].append(E_T_G_1red_approx(task_t_rv, d, k) )
      
      stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l=k, n=k+1, num_run=1000)
      E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )
      # C_key = 'E_C' if not w_cancel else 'E_C_wc'
      # E_C = sum(stat_id__trial_sampleavg_l_m[C_key] )/len(stat_id__trial_sampleavg_l_m[C_key] )
      k__x_sim_l_m[k].append(d)
      k__y_sim_l_m[k].append(E_T)
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  color = iter(cm.rainbow(numpy.linspace(0, 1, int(K/k_step)+1) ) )
  for k in range(2, K, k_step):
  # for k in range(K, K+1, k_step):
    m = next(marker)
    c = next(color)
    plot.plot(k__x_l_m[k], k__y_l_m[k], label=r'k:{}'.format(k), color=c, marker=m, linestyle='', mew=2)
    # plot.plot(k__x_l_m[k], k__y_approx_l_m[k], label=r'~ k:{}'.format(k), color=c, alpha=0.5, marker=m, linestyle='', mew=2)
    plot.plot(k__x_sim_l_m[k], k__y_sim_l_m[k], label=r'sim k:{}'.format(k), color=c, alpha=0.5, marker=m, linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$\Delta$ (s)')
  plot.ylabel(r'$E[T]$ (s)')
  plot.title(r'$X_i \sim$ {}'.format(task_t_rv) )
  plot.savefig("plot_arepeat_E_T_G_1red_{}.png".format(task_t_rv) )
  plot.gcf().clear()
  log(WARNING, "done; K= {}, task_t_rv= {}".format(K, task_t_rv) )

# ##################  X_i ~ Exp(mu), (l, k, n, \Delta)  ################ #
def plot_arepeat_k_l_n__E_C_vs_E_T():
  mu = 1
  k = 10
  N = 20
  task_t_rv = Exp(mu)
  
  l__x_l_m, l__y_l_m = {}, {}
  l__x_sim_l_m, l__y_sim_l_m = {}, {}
  
  l_step = 2
  for l in range(k, N, l_step):
    # l__x_sim_l_m[l] = []
    # l__y_sim_l_m[l] = []
    l__y_l_m[l] = []
    l__x_l_m[l] = []
    for d in numpy.arange(0.1, 3*H(k)/mu, 0.1):
      # stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, N, num_run=100000)
      # E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )
      # E_C = sum(stat_id__trial_sampleavg_l_m['E_C'] )/len(stat_id__trial_sampleavg_l_m['E_C'] )
      # l__x_sim_l_m[l].append(E_T)
      # l__y_sim_l_m[l].append(E_C)
      
      l__x_l_m[l].append(E_T_exp_k_l_n(mu, d, k, l, N) )
      l__y_l_m[l].append(E_C_exp_k_l_n(mu, d, k, l, N, w_cancel=False) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  color = iter(cm.rainbow(numpy.linspace(0, 1, int((N-k)/l_step)+1) ) )
  for l in range(k, N, l_step):
    plot.plot(l__x_l_m[l], l__y_l_m[l], label=r'l={}'.format(l), color=next(dark_color), marker=next(marker), linestyle='-', mew=2)
    # plot.plot(l__x_sim_l_m[l], l__y_sim_l_m[l], label=r'l:{}'.format(l), color=next(light_color), alpha=0.7, marker=next(marker), linestyle='', mew=2)
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.xlabel(r'$E[T]$ (s)')
  plot.ylabel(r'$E[C]$ (s)')
  plot.title(r'$X \sim Exp(\mu={}), k= {}, n= {}$'.format(mu, k, N) )
  plot.savefig("plot_arepeat_k_l_n__E_C_vs_E_T__n_{}.png".format(N) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}, N= {}".format(k, N) )

def plot_arepeat_E_T_k_n_vs_k_l_n():
  mu = 1
  k = 10
  N = 20
  
  l__y_l_m = {}
  l__x_l_m = {}
  
  l_step = 2
  for l in range(k, N, l_step):
    l__y_l_m[l] = []
    l__x_l_m[l] = []
    for d in numpy.arange(0.1, 2*H(k)/mu, 0.1):
      l__x_l_m[l].append(d)
      # if l == k:
      #   # l__y_l_m[l].append(E_T_exp_k_l_n(mu, d, k, l, N) )
      #   l__y_l_m[l].append(E_T_exp_k_n(mu, d, k, N) )
      # else:
      l__y_l_m[l].append(E_T_exp_k_l_n(mu, d, k, l, N) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  color = iter(cm.rainbow(numpy.linspace(0, 1, int((N-k)/l_step)+1) ) )
  for l in range(k, N, l_step):
    m = next(marker)
    c = next(color)
    plot.plot(l__x_l_m[l], l__y_l_m[l], label=r'l:{}'.format(l), color=c, marker=m, linestyle='', mew=2)
  # plot.legend()
  plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.xlabel(r'$\Delta$ (s)')
  plot.ylabel(r'$E[T]$ (s)')
  plot.title(r'$\mu$= {}, k= {}, N= {}'.format(mu, k, N) )
  plot.savefig("plot_arepeat_E_T_k_n_vs_k_l_n__N_{}.png".format(N) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}, N= {}".format(k, N) )

def plot_arepeat_E_T_k_l_n():
  mu = 1
  l = 11 #  91 # 11
  k = 10 # 90 # 10
  N = 20 # 120 # 20
  d = 1
  task_t_rv = Exp(mu)
  
  n__y_sim_l_m, n__y_l_m, n__y_approx_l_m = {}, {}, {}
  n__x_l_m = {}
  
  l_step = int((N-l)/4)
  for n in range(l, N, l_step):
    n__y_sim_l_m[n] = []
    n__y_l_m[n] = []
    n__y_approx_l_m[n] = []
    n__x_l_m[n] = []
    for d in numpy.arange(0.1, 3/mu*(H(l)-H(l-k) ), 0.1):
      n__x_l_m[n].append(d)
      
      # stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=10000)
      # n__y_sim_l_m[n].append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
      n__y_l_m[n].append(E_T_exp_k_l_n(mu, d, k, l, n) )
      n__y_approx_l_m[n].append(E_T_exp_k_l_n_approx(mu, d, k, l, n) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  color = iter(cm.rainbow(numpy.linspace(0, 1, int((N-l)/l_step)+1) ) )
  for n in range(l, N, l_step):
    m = next(marker)
    c = next(color)
    # plot.plot(n__x_l_m[n], n__y_sim_l_m[n], label=r'sim, n:{}'.format(n), color=c, alpha=0.5, marker=m, linestyle='', mew=2)
    plot.plot(n__x_l_m[n], n__y_l_m[n], label=r'n:{}'.format(n), color=c, marker=m, linestyle='', mew=2)
    plot.plot(n__x_l_m[n], n__y_approx_l_m[n], label=r'~ n:{}'.format(n), color=c, alpha=0.5, marker=m, linestyle='', mew=2)
  # plot.legend()
  plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.xlabel(r'$\Delta$ (s)')
  plot.title(r'$X_i \sim Exp(\mu={})$, k= {}, l= {}'.format(mu, k, l) )
  plot.ylabel(r'$E[T]$ (s)')
  plot.savefig("plot_arepeat_k_l_{}.png".format(l) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}, l= {}".format(k, l) )

# ############################  (l=k, k, n, \Delta)  ########################## #
def plot_arepeat_conf_k_n():
  mu = 1
  k = 10
  n = 20 # 15
  # d = H(k)/mu * 0.3
  L = 2
  
  d_l, prob_T_k_n_geq_L_l, Pr_T_k_n_geq_L_approx_l = [], [], []
  # for d in numpy.arange(0.1, 3*H(k)/mu, 0.1):
  for d in numpy.arange(0.1, L, 0.1):
    d_l.append(d)
    
    prob_T_k_n_geq_L_l.append(prob_T_k_n_geq_t(mu, d, k, n, L) )
    Pr_T_k_n_geq_L_approx_l.append(prob_T_k_n_geq_t_approx(mu, d, k, n, L) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(d_l, prob_T_k_n_geq_L_l, color='black', label=r'$Pr\{T \geq L\}$', marker=next(marker), linestyle='', mew=2)
  plot.plot(d_l, Pr_T_k_n_geq_L_approx_l, color='brown', label=r'$Pr\{T_{alt} \geq L\}$', marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$\Delta$')
  plot.ylabel(r'$Pr\{T \geq L\}$')
  plot.title(r'mu= {}, k= {}, n= {}, L= {}'.format(mu, k, n, L) )
  plot.savefig("plot_arepeat_conf_k_{}_n_{}.png".format(k, n) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}, n= {}".format(k, n) )

def plot_arepeat_dist_k_n():
  mu = 1
  k = 10
  n = 20
  d = H(k)/mu * 0.3
  
  t_l, prob_T_k_n_geq_t_l, prob_T_k_n_geq_t_approx_l = [], [], []
  # for d in numpy.arange(0.1, 3*H(k)/mu, 0.1):
  for t in numpy.arange(0, 3*d, 0.1):
    t_l.append(t)
    
    prob_T_k_n_geq_t_l.append(prob_T_k_n_geq_t(mu, d, k, n, t) )
    prob_T_k_n_geq_t_approx_l.append(prob_T_k_n_geq_t_approx(mu, d, k, n, t) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(t_l, prob_T_k_n_geq_t_l, color='black', label=r'$Pr\{T \geq t\}$', marker=next(marker), linestyle='', mew=2)
  plot.plot(t_l, prob_T_k_n_geq_t_approx_l, color='brown', label=r'$Pr\{T_{alt} \geq t\}$', marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$t$')
  plot.ylabel(r'$Pr\{T \geq t\}$')
  plot.title(r'mu= {}, k= {}, n= {}, $\Delta$= {}'.format(mu, k, n, d) )
  plot.savefig("plot_arepeat_dist_k_{}_n_{}.png".format(k, n) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}, n= {}".format(k, n) )

def plot_dE_T_shiftedexp_k_n_dk():
  D = 1
  mu = 1
  k = 3 # 10
  N = k + 1 # 20 # 20
  d = 1
  
  n__k_l_m, n__d_E_T_dk_l_m = {}, {}
  
  k_step = 2
  for n in range(k, N, k_step):
    n__d_E_T_dk_l_m[n] = []
    n__k_l_m[n] = []
    for k_ in numpy.arange(1, n+1, 1):
      n__k_l_m[n].append(k_)
      n__d_E_T_dk_l_m[n].append(d_E_T_shiftedexp_k_n_dk(D, mu, d, k_, n) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  color = iter(cm.rainbow(numpy.linspace(0, 1, int(N/k_step) ) ) )
  for n in range(k, N, k_step):
    m = next(marker)
    c = next(color)
    plot.plot(n__k_l_m[n], n__d_E_T_dk_l_m[n], label=r'n:{}, $d(E[T])/dk$'.format(n), color=c, marker=m, linestyle='', mew=2)
  # plot.legend()
  plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.xlabel(r'$k$')
  plot.title(r'D= {}, $\mu$= {}, $\Delta$= {}'.format(D, mu, d) )
  plot.ylabel("E[T] (s)")
  plot.savefig("plot_dE_T_shiftedexp_k_n_dk__d_{}.png".format(d) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(k) )

def plot_arepeat_shiftedexp_k_n():
  D = 1
  mu = 1
  k = 10
  N = 15 # 20
  d = 1
  
  d_l = []
  n__E_T_l_m, n__E_T_approx_l_m, n__E_T_sim_l_m = {}, {}, {}
  n__k_l_m = {}
  
  varying_d = True # False
  varying_k = False # True
  
  first_loop = True
  k_step = 2
  for n in range(k, N, k_step):
    n__E_T_l_m[n] = []
    n__E_T_approx_l_m[n] = []
    n__E_T_sim_l_m[n] = []
    if varying_d:
      for d in numpy.arange(0, 2*H(k)/mu, 0.5):
        if first_loop:
          d_l.append(d)
        n__E_T_l_m[n].append(E_T_shiftedexp_k_n(D, mu, d, k, n) )
        # n__E_T_approx_l_m[n].append(E_T_exp_k_n_approx(mu, d, k, n) )
        
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(Exp(mu, D/k), d, k, k, n, num_run=10000)
        n__E_T_sim_l_m[n].append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
      first_loop = False
    elif varying_k:
      n__k_l_m[n] = []
      for k_ in numpy.arange(1, n+1, 1):
        n__k_l_m[n].append(k_)
        n__E_T_l_m[n].append(E_T_shiftedexp_k_n(D, mu, d, k_, n) )
  for n in range(k, N, k_step):
    if varying_d:
      plot.plot(d_l, n__E_T_l_m[n], label=r'$n={}$'.format(n), color=next(dark_color), linestyle='-', mew=2)
      # plot.plot(d_l, n__E_T_approx_l_m[n], label=r'$n={}$'.format(n), color=next(light_color), marker=next(marker), linestyle='', mew=2)
      plot.plot(d_l, n__E_T_sim_l_m[n], label=r'$n={}$'.format(n), color=next(light_color), marker=next(marker), linestyle='', mew=2)
    elif varying_k:
      plot.plot(n__k_l_m[n], n__E_T_l_m[n], label=r'$n={}$'.format(n), color=next(dark_color), linestyle='-', mew=2)
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  if varying_d:
    plot.xlabel(r'$\Delta$ (s)')
    plot.title(r'$X \sim SExp(D\k, \mu), D={}, \mu={}, k= {}$'.format(D, mu, k) )
  elif varying_k:
    plot.xlabel(r'$k$')
    plot.title(r'$X \sim SExp(D\k, \mu), D={}, \mu={}, \Delta= {}$'.format(D, mu, d) )
  plot.ylabel(r'$E[T]$ (s)')
  if varying_d:
    plot.savefig("plot_arepeat_shiftedexp_k_n__k_{}.png".format(k) )
  elif varying_k:
    plot.savefig("plot_arepeat_shiftedexp_k_n__d_{}.png".format(d) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(k) )

def plot_arepeat_k_n():
  K = 10
  N = 14
  D = 30
  mu = 1
  a, loc = 2, 3
  task_t = "SExp" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim SExp(D/k, \mu), D={}, \mu={}'.format(D, mu)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\alpha={}, \lambda={})'.format(a, loc)
  
  def plot_(k, n=None, c=None, sim=False):
    if task_t == "Exp": task_t_rv = Exp(mu)
    elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
    elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    
    d_l = []
    y_l, y_approx_l, y_sim_l = [], [], []
    u_l = 3*H(K)/mu + 1
    for d in numpy.arange(0, u_l, 0.5):
      d_l.append(d)
      # if n is not None:
      #   y_l.append(E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n) )
      #   # y_l.append(E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n, w_cancel=True) )
      # elif c is not None:
      #   y_l.append(E_T_k_c(task_t, D, mu, loc, a, d, k, c) )
      #   # y_l.append(E_C_k_l_n(task_t, D, mu, loc, a, d, k, c, w_cancel=True) )
      if sim:
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=10000)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=10000)
        # y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
        y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] ) )
        # y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C'] )/len(stat_id__trial_sampleavg_l_m['E_C'] ) )
    # plot.plot(d_l, y_l, label=r'$n={}$'.format(n), color=next(dark_color), linestyle='-', mew=2)
    # plot.plot(d_l, y_approx_l, label=r'$n={}$'.format(n), color=next(light_color), marker=next(marker), linestyle='', mew=2)
    plot.plot(d_l, y_sim_l, label=r'$n={}$'.format(n), color=next(light_color), marker=next(marker), linestyle='', mew=2)
  
  sim = False
  k_step = 2
  for n in [K, K+1, *range(K+k_step, N, k_step)]:
    plot_(K, n=n, sim=sim)
  plot_(K, n=20, sim=sim)
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  plot.xlabel(r'$\Delta$ (s)')
  # plot.ylabel(r'$E[T]$ (s)')
  plot.ylabel(r'$E[C]$ (s)')
  plot.savefig("plot_arepeat_k_n_{}_k_{}.png".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )
  
def plot_arepeat_k_nc_E_C_vs_E_T():
  K = 400 # 1050 # 1000 # 400 # 10
  N = 12
  D = 30
  mu = 1
  a, loc = 2, 3
  
  task_t = "Google" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim SExp(D/k, \mu), D={}, \mu={}'.format(D, mu)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\alpha={}, \lambda={})'.format(a, loc)
  elif task_t == "Google": task_t_in_latex = r'X \sim Google'
  
  def plot_(k, n=None, c=None, sim=False):
    if task_t == "Exp": task_t_rv = Exp(mu)
    elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
    elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    elif task_t == "Google": task_t_rv = Google(k)
    l = k
    E_T_l, E_T_sim_l, E_T_approx_l = [], [], []
    E_Cwcancel_l, E_Cwcancel_sim_l, E_C_l, E_C_sim_l = [], [], [], []
    
    u_l = 10*task_t_rv.mean()
    for d in numpy.linspace(0, u_l, 20):
      # E_Cwcancel_l.append(E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n, w_cancel=True) )
      # E_C_l.append(E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n, w_cancel=False) )
      # E_T_l.append(E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n) )
      
      if sim:
        if n is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=10000)
        elif c is not None:
          stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=10000)
        E_T_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
        E_Cwcancel_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] ) )
      # E_C_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C'] )/len(stat_id__trial_sampleavg_l_m['E_C'] ) )
    # plot.plot(E_T_l, E_Cwcancel_l, label=r'w/ c,n={}'.format(n), color=next(dark_color), linestyle='-', mew=2)
    plot.plot(E_T_sim_l, E_Cwcancel_sim_l, label=r'w/ c,n={}'.format(n), marker=next(marker), color=next(dark_color), linestyle=':', mew=2)
    # plot.plot(E_T_l, E_C_l, label=r'no-c,n={}'.format(n), color=next(dark_color), linestyle='--', mew=2)
    # plot.plot(E_T_sim_l, E_C_sim_l, label=r'no-c,n={}'.format(n), marker=next(marker), color=next(dark_color), linestyle=':', mew=2)
  sim = True
  plot_(K, n=math.floor(K*1.1), sim=sim)
  plot_(K, c=1, sim=sim)
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.8, 0.6) )
  plot.xlabel(r'$E[T]$ (s)')
  # plot.ylabel(r'$E[C]$ (s)')
  plot.ylabel(r'$E[C^c]$ (s)')
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  fig = plot.gcf()
  fig.tight_layout()
  plot.savefig("plot_arepeat_k_nc_E_C_vs_E_T_{}_k_{}.png".format(task_t, K) )
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
          y_wrelaunch_l.append(E_C_pareto_k_n(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
          # y_wrelaunch_approx_l.append(E_C_pareto_k_n_approx(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
        elif c is not None:
          y_wrelaunch_l.append(E_C_pareto_k_c(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
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
      d = Delta_for_min_E_T_pareto_k_n_wrelaunch(loc, a, k, n)
      # plot.axvline(x=d, label=r'$\Delta^*$, n={}'.format(n), color=c)
      y = E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=True)
      plot.plot(d, y, color='red', label=r'Approx minimum', marker='*', zorder=1, ms=10, mew=3)
    # print("y_wrelaunch_l= {}".format(y_wrelaunch_l) )
    # label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    label = 'w/ relaunch'
    plot.plot(x_l, y_wrelaunch_l, label=label, color=c, marker=next(marker), linestyle='--', zorder=0, mew=2)
    # plot.plot(x_l, y_wrelaunch_approx_l, label='approx, {}'.format(label), color=next(dark_color), marker=next(marker), linestyle='', mew=2)
    
    if sim:
      plot.plot(x_l, y_sim_wrelaunch_l, label=r'sim, {}'.format(label), color=next(light_color), marker=next(marker), linestyle='', mew=2)
    if E_T:
      if n == k: plot.plot(x_l, y_norelaunch_l, label='no-relaunch', color=next(dark_color), linestyle='--', mew=2)
  
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
  # plot.title(r'$X \sim Pareto(\lambda= {}, \alpha= {}), k= {}$'.format(loc, a, K) )
  plot.title(r'$X \sim Pareto(\lambda= {}, \alpha)$'.format(loc) )
  
  # plot.xlabel(r'$\Delta$ (s)', fontsize=13)
  # y_label = r'$E[T]$' if E_T else r'$E[C]$'
  # plot.ylabel(r'{} (s)'.format(y_label), fontsize=13)
  plot.xlabel(r'$\alpha$', fontsize=12)
  plot.ylabel('\n'.join(textwrap.wrap(r'Maximum percetange reduction in $E[T]$ by task relaunch', 40) ), fontsize=12)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_arepeat_k_nc_wrelaunch_{}_k_{}.pdf".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

def plot_arepeat_k_nc_wrelaunch_E_C_vs_E_T():
  K = 100
  loc, a = 3, 2
  task_t = "Pareto"
  w_cancel = True # False # True
  
  def plot_(k, n=None, c=None):
    E_T_l, E_C_l = [], []
    # for d in [*numpy.linspace(0, 5*loc, 50), *numpy.linspace(5*loc, 200*loc, 100)]:
    for d in numpy.logspace(-3, 2, 100):
      if n is not None:
        E_T_l.append(E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=True) )
        E_C_l.append(E_C_pareto_k_n(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True) )
      elif c is not None:
        E_T_l.append(E_T_pareto_k_c_wrelaunch(loc, a, d, k, c, w_relaunch=True) )
        E_C_l.append(E_C_pareto_k_c_wrelaunch(loc, a, d, k, c, w_cancel=w_cancel, w_relaunch=True) )
    color = next(dark_color)
    if n is not None and n == k:
      x = E_T_pareto_k_n(loc, a, 0, k, n, w_relaunch=True)
      y = E_C_pareto_k_n(loc, a, 0, k, n, w_cancel=w_cancel, w_relaunch=True)
      if k == 100 and a == 2:
        # plot.annotate(r'$\Delta=0$ or $\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=12)
        plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+1, y+4), fontsize=12)
        plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x-1.5, y-25), fontsize=12)
      else:
        # plot.plot(x, y, color='black', marker='*', ms=8, mew=3)
        plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x, y), fontsize=12)
      
      x = E_T_pareto_k_n(loc, a, 20*loc, k, n, w_relaunch=True)
      y = E_C_pareto_k_n(loc, a, 20*loc, k, n, w_cancel=w_cancel, w_relaunch=True)
      # if k == 100 and a == 2:
      #   plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+1, y), fontsize=12)
      # else:
      #   # plot.plot(x, y, color='red', marker='*', ms=8, mew=3)
      #   plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x, y), fontsize=12)
      if k == 100 and a == 2:
        d = Delta_for_min_E_T_pareto_k_n_wrelaunch(loc, a, k, n)
        x = E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=True)
        y = E_C_pareto_k_n(loc, a, d, k, n, w_cancel=w_cancel, w_relaunch=True)
        plot.plot(x, y, color='red', label=r'Approx minimum', marker='*', zorder=1, ms=10, mew=3)
        plot.annotate(r'$\Delta={0:.2f}$'.format(d), xy=(x, y), xytext=(x-1.5, y+20), fontsize=12)
    elif n is not None and n > k:
      x = E_T_pareto_k_n(loc, a, 0, k, n, w_relaunch=True)
      y = E_C_pareto_k_n(loc, a, 0, k, n, w_cancel=w_cancel, w_relaunch=True)
      if k == 100 and a == 2:
        plot.plot(x, y, color='black', marker='o', zorder=1, ms=8, mew=2)
        # plot.annotate(r'$\Delta=0$', xy=(x, y), xytext=(x+1, y+5), color=color, fontsize=13)
      x = E_T_pareto_k_n(loc, a, 10**5, k, n, w_relaunch=True)
      y = E_C_pareto_k_n(loc, a, 10**5, k, n, w_cancel=w_cancel, w_relaunch=True)
      if n == 101 and k == 100 and a == 2:
        plot.annotate(r'$\Delta \rightarrow \infty$', xy=(x, y), xytext=(x+0.5, y), fontsize=13)
    elif c is not None:
      x = E_T_pareto_k_c_wrelaunch(loc, a, 0, k, c, w_relaunch=True)
      y = E_C_pareto_k_c_wrelaunch(loc, a, 0, k, c, w_cancel=w_cancel, w_relaunch=True)
      plot.plot(x, y, color='black', marker='o', zorder=1, ms=6, mew=2)
    label = r'$n={}$'.format(n) if n is not None else r'$c={}$'.format(c)
    plot.plot(E_T_l, E_C_l, label=label, color=color, marker=next(marker), zorder=0, linestyle=':', mew=2)
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
  """
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
  
  plot.legend()
  # plot.legend(loc='center left', bbox_to_anchor=(0.9, 0.5) )
  plot.title(r'$X \sim Pareto(\lambda= {}, \alpha= {}), k= {}$'.format(loc, a, K) )
  plot.xlabel(r'Expected Latency $E[T]$ (s)', fontsize=13)
  y = r'$E[C^c]$' if w_cancel else r'$E[C]$'
  plot.ylabel('Expected Cost {} (s)'.format(y), fontsize=13)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_arepeat_k_nc_wrelaunch_E_C_vs_E_T_{}_k_{}.pdf".format(task_t, K), bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

# #########################  Reped vs. Coded  ##################### #
def plot_reped_vs_coded(w_cancel):
  K = 10
  D = 30
  mu = 0.5
  a, loc = 2, 3
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

def plot_reped_vs_coded_nodelay(loc, a):
  w_cancel=True
  K = 400 # 10
  D = 30
  mu = 0.5
  # loc, a = 3, 2 # 1.5
  num_run = 10 # 10000
  task_t = "Google" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim SExp(D={}/k, \mu={})'.format(D, mu)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  elif task_t == "Google": task_t_in_latex = r'X \sim Google'
  mew = 3
  def plot_(k=K, n=0, c=0):
    l = k
    if task_t == "Exp": task_t_rv = Exp(mu)
    elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
    elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    elif task_t == "Google": task_t_rv = Google(k)
    
    E_T_l, E_T_sim_l, E_C_l, E_C_sim_l = [], [], [], []
    d = 0
    if c:
      color = next(dark_color)
      for c_ in range(c+1):
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c_, num_run)
        E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )/60/60
        key = 'E_C_wc' if w_cancel else 'E_C'
        E_C = sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] )/60/60
        E_T_sim_l.append(E_T)
        E_C_sim_l.append(E_C)
        if task_t == "Google":
          if c_ == 0:
            plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T, E_C), xytext=(E_T-10, E_C+200) )
          else:
            plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+10, E_C+100), color=color)
        
        # E_T = E_T_k_c(task_t, D, mu, loc, a, d, k, c_)
        # E_C = E_C_k_c(task_t, D, mu, loc, a, d, k, c_, w_cancel=w_cancel)
        # E_T_l.append(E_T)
        # E_C_l.append(E_C)
        # if c_ > 0:
        #   if task_t == "SExp":
        #     plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.1, E_C), color=color)
        #   elif task_t == "Pareto":
        #     if a == 1.2:
        #       plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+1, E_C+1), color=color)
        #     elif a == 1.5:
        #       plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.3, E_C+0.3), color=color)
        #     elif a == 2:
        #       plot.annotate(r'$c={}$'.format(c_), xy=(E_T, E_C), xytext=(E_T+0.3, E_C+0.3), color=color)
      # plot.plot(E_T_l, E_C_l, label='Replication', color=color, marker=next(marker), zorder=0, mew=2, linestyle=':')
      plot.plot(E_T_sim_l, E_C_sim_l, label='Simulation, replication', color=color, marker=next(marker), linestyle=':', mew=2)
    elif n:
      color = next(dark_color)
      for n_ in range(k, n+1):
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, k, n_, num_run)
        E_T = sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] )/60/60
        key = 'E_C_wc' if w_cancel else 'E_C'
        E_C = sum(stat_id__trial_sampleavg_l_m[key] )/len(stat_id__trial_sampleavg_l_m[key] )/60/60
        E_T_sim_l.append(E_T)
        E_C_sim_l.append(E_C)
        if task_t == "Google":
          if n_ > K and n_ < k+3:
            plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T-10, E_C-200), color=color)
          elif n_ != k and n_ % k == 0:
            plot.plot(E_T, E_C, 'x', color="blue", mew=mew, zorder=3)
            plot.annotate(r'$n={}$'.format(n_), xy=(E_T, E_C), xytext=(E_T+10, E_C+100), color=color)
        
        # E_T = E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n_)
        # E_C = E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n_, w_cancel=w_cancel)
        # E_T_l.append(E_T)
        # E_C_l.append(E_C)
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
      # plot.plot(E_T_l, E_C_l, label='Coding', color=color, zorder=1, marker=next(marker), mew=1, linestyle=':')
      plot.plot(E_T_sim_l, E_C_sim_l, label='Simulation, coding', color=color, zorder=1, marker=next(marker), linestyle=':', mew=1)
  plot_(c=5)
  plot_(n=6*K)
  #
  E_T_nored = E_T_k_c(task_t, D, mu, loc, a, d=0, k=K, c=0)
  E_C_nored = E_C_k_c(task_t, D, mu, loc, a, d=0, k=K, c=0, w_cancel=w_cancel)
  if task_t == "SExp":
    plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-0.6, E_C_nored+9) )
  elif task_t == "Pareto":
    plot.axhline(y=E_C_nored, color='k', alpha=0.4, linestyle='--')
    if a == 1.2:
      plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-20, E_C_nored+3) )
    elif a == 1.5:
      plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-6, E_C_nored+3) )
    elif a == 2:
      plot.annotate('No redundancy \n $c=0$, $n={}$'.format(K), xy=(E_T_nored, E_C_nored), xytext=(E_T_nored-1.2, E_C_nored+3) )
  plot.legend()
  # plot.legend(loc='lower right')
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  # plot.xlabel(r'Expected Latency $E[T]$ (s)', fontsize=12)
  plot.xlabel(r'Expected latency $E[T]$ (hour)', fontsize=12)
  plot.ylabel(r'Expected cost $E[C^c]$ (hour)' if w_cancel else r'$E[C]$', fontsize=12)
  fig = plot.gcf()
  fig.tight_layout()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  # plot.savefig("plot_reped_vs_coded_nodelay_{}_k_{}__a_{}.png".format(task_t, K, a) )
  if task_t == "Pareto":
    plot.savefig("plot_reped_vs_coded_nodelay_{}_k_{}_a_{}_.pdf".format(task_t, K, a) )
  else:
    plot.savefig("plot_reped_vs_coded_nodelay_{}_k_{}.pdf".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

if __name__ == "__main__":
  # plot_send_n_w_drop()
  # plot_arepeat_dist_k_n()
  # plot_arepeat_conf_k_n()
  # plot_arepeat_E_T_G_1red()
  # plot_arepeat_E_T_vs_E_C_Gred()
  
  # plot_arepeat_k_l_n__E_C_vs_E_T()
  # plot_arepeat_E_T_k_n_vs_k_l_n()
  # plot_arepeat_E_T_k_l_n()
  
  # plot_dE_T_shiftedexp_k_n_dk()
  # plot_arepeat_shiftedexp_k_n()
  
  # plot_reped_vs_coded(w_cancel=True)
  plot_reped_vs_coded_nodelay(loc=3, a=2)
  # plot_reped_vs_coded_nodelay(loc=3, a=1.2)
  # plot_reped_vs_coded_nodelay(loc=3, a=1.5)
  
  # plot_arepeat_k_n()
  # plot_arepeat_k_nc_E_C_vs_E_T()
  # plot_arepeat_k_nc_wrelaunch()
  # plot_arepeat_k_nc_wrelaunch_E_C_vs_E_T()
  
  # plot_arepeat_k_c()
  # plot_arepeat_k_c_E_C_vs_E_T()
  