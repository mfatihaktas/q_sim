# ##################  Send n any k is enough, X_i ~ Exp(mu), each packet drops ~ Exp(gamma)  ############## #
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
  plot.title(r'$\mu$= {}, k= {}'.format(mu, k) )
  plot.xlabel(r'$\gamma$')
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

def plot_arepeat_Pr_T_g_t():
  K = 10
  D, mu = 30, 1
  loc, a = 3, 2
  task_t = "Exp" # "Exp" # "SExp" # "Pareto"
  task_t_rv, task_t_in_latex = None, None
  
  if task_t == "Exp": task_t_in_latex = r'X \sim Exp(\mu={})'.format(mu)
  elif task_t == "SExp": task_t_in_latex = r'X \sim SExp(D/k, \mu), D={}, \mu={}'.format(D, mu)
  elif task_t == "Pareto": task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha={})'.format(loc, a)
  
  def plot_(k, d, n=None, c=None, sim=False):
    if task_t == "Exp": task_t_rv = Exp(mu)
    elif task_t == "SExp": task_t_rv = Exp(mu, D/k)
    elif task_t == "Pareto": task_t_rv = Pareto(a, loc=loc)
    
    t_l = []
    y_l, y_approx_l, y_sim_l = [], [], []
    for t in numpy.linspace(0, 3*H(K)/mu, 100):
      t_l.append(t)
      if n is not None:
        y_l.append(Pr_T_g_t_k_n(task_t, D, mu, loc, a, d, k, n, t) )
        # y_approx_l.append(Pr_T_g_t_Exp_k_n_approx(mu, d, k, n, t) )
      elif c is not None:
        y_l.append(Pr_T_g_t_k_c(task_t, D, mu, loc, a, d, k, c, t) )
      # if sim:
      #   if n is not None:
      #     stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(task_t_rv, d, k, l, n, num_run=10000)
      #   elif c is not None:
      #     stat_id__trial_sampleavg_l_m = sim_arepeat_k_c(task_t_rv, d, k, c, num_run=10000)
      #   # y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_T'] )/len(stat_id__trial_sampleavg_l_m['E_T'] ) )
      #   y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C_wc'] )/len(stat_id__trial_sampleavg_l_m['E_C_wc'] ) )
      #   # y_sim_l.append(sum(stat_id__trial_sampleavg_l_m['E_C'] )/len(stat_id__trial_sampleavg_l_m['E_C'] ) )
    label = r'$n={}$, $\Delta={}$'.format(n, d) if n is not None else r'$c={}$'.format(c)
    # print("y_l= {}".format(y_l) )
    plot.plot(t_l, y_l, label=label, color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
    # plot.plot(t_l, y_approx_l, label=r'Approx, {}'.format(label), color=next(light_color), marker=next(marker), linestyle=':', mew=2)
    # plot.plot(d_l, y_sim_l, label=r'$n={}$'.format(n), color=next(light_color), marker=next(marker), linestyle='', mew=2)
  
  sim = False
  plot_(K, d=0, n=K+1, sim=sim)
  plot_(K, d=1/mu, n=K+1, sim=sim)
  plot_(K, d=2/mu, n=K+1, sim=sim)
  plot_(K, d=3/mu, n=K+1, sim=sim)
  plot_(K, d=4/mu, n=K+1, sim=sim)
  
  plot.legend()
  plot.title(r'${}, k= {}$'.format(task_t_in_latex, K) )
  plot.xlabel(r'$t$')
  plot.ylabel(r'$Pr\{T \geq t\}$')
  plot.savefig("plot_arepeat_Pr_T_g_t_{}_k_{}.png".format(task_t, K) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(K) )

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