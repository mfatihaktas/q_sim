import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
from cycler import cycler

from auto_repeat_models import *

# ------------------------  Checks on Pei's work  --------------------- #
def plot_autorepeat_conf_k_n():
  mu = 1
  k = 10
  n = 20 # 15
  # d = H(k)/mu * 0.3
  L = 2
  
  d_l, prob_T_k_n_geq_L_l, prob_T_k_n_geq_L_alt_l = [], [], []
  # E_T_k_n_11_l, E_T_k_n_r_1_5_alt_l, E_T_k_n_r_1_l, E_T_k_n_r_1_alt_l = [], [], [], []
  # for d in numpy.arange(0.1, 3*H(k)/mu, 0.1):
  for d in numpy.arange(0.1, L, 0.1):
    d_l.append(d)
    
    prob_T_k_n_geq_L_l.append(prob_T_k_n_geq_t(mu, d, k, n, L) )
    prob_T_k_n_geq_L_alt_l.append(prob_T_k_n_geq_t_alt(mu, d, k, n, L) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(d_l, prob_T_k_n_geq_L_l, color='black', label=r'$Pr\{T \geq L\}$', marker=next(marker), linestyle='', mew=2)
  plot.plot(d_l, prob_T_k_n_geq_L_alt_l, color='brown', label=r'$Pr\{T_{alt} \geq L\}$', marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$d$')
  plot.ylabel(r'$Pr\{T \geq L\}$')
  plot.title(r'mu= {}, k= {}, n= {}, L= {}'.format(mu, k, n, L) )
  plot.savefig("plot_autorepeat_conf_k_{}_n_{}.png".format(k, n) )
  log(WARNING, "done; k= {}, n= {}".format(k, n) )

def plot_autorepeat_dist_k_n():
  mu = 1
  k = 10
  n = 10
  d = H(k)/mu * 0.3
  
  t_l, prob_T_k_n_geq_t_l, prob_T_k_n_geq_t_alt_l = [], [], []
  # E_T_k_n_11_l, E_T_k_n_r_1_5_alt_l, E_T_k_n_r_1_l, E_T_k_n_r_1_alt_l = [], [], [], []
  # for d in numpy.arange(0.1, 3*H(k)/mu, 0.1):
  for t in numpy.arange(0, 3*d, 0.1):
    t_l.append(t)
    
    prob_T_k_n_geq_t_l.append(prob_T_k_n_geq_t(mu, d, k, n, t) )
    prob_T_k_n_geq_t_alt_l.append(prob_T_k_n_geq_t_alt(mu, d, k, n, t) )
    # E_T_k_n_11_l.append(E_T_k_n(mu, d, k, int(k*1.5) ) )
    # E_T_k_n_r_1_5_alt_l.append(E_T_k_n_alt(mu, d, k, int(k*1.5) ) )
    # E_T_k_n_r_1_l.append(E_T_k_n(mu, d, k, k) )
    # E_T_k_n_r_1_alt_l.append(E_T_k_n_alt(mu, d, k, k) )
    # E_T_k_n_lb_l.append(E_T_k_n_lb(mu, d, k, n) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(t_l, prob_T_k_n_geq_t_l, color='black', label=r'$Pr\{T \geq t\}$', marker=next(marker), linestyle='', mew=2)
  plot.plot(t_l, prob_T_k_n_geq_t_alt_l, color='brown', label=r'$Pr\{T_{alt} \geq t\}$', marker=next(marker), linestyle='', mew=2)
  # plot.plot(d_l, E_T_k_n_11_l, color='magenta', label=r'$E[T]$, r:1.5', marker=next(marker), linestyle='', mew=2)
  # plot.plot(d_l, E_T_k_n_r_1_5_alt_l, color='pink', label=r'$E[T_{alt}]$, r:1.5', marker=next(marker), linestyle='', mew=2)
  # plot.plot(d_l, E_T_k_n_r_1_l, color='red', label=r'$E[T]$, r:1', marker=next(marker), linestyle='', mew=2)
  # plot.plot(d_l, E_T_k_n_lb_l, color='green', label=r'$E[\hat{T}_{LB}]$', marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$t$')
  plot.ylabel(r'$Pr\{T \geq t\}$')
  plot.title(r'mu= {}, k= {}, n= {}, d= {}'.format(mu, k, n, d) )
  plot.savefig("plot_autorepeat_dist_k_{}_n_{}.png".format(k, n) )
  log(WARNING, "done; k= {}, n= {}".format(k, n) )
  
def plot_autorepeat_k_n():
  mu = 1
  k = 10
  N = 20
  
  d_l = []
  n__E_T_l_m, n__E_T_model_l_m = {}, {}
  
  
  first_loop = True
  k_step = 1
  for n in range(k, N, k_step):
    n__E_T_l_m[n] = []
    n__E_T_model_l_m[n] = []
    for d in numpy.arange(0.1, 2*H(k)/mu, 0.1):
      if first_loop:
        d_l.append(d)
      n__E_T_l_m[n].append(E_T_k_n(mu, d, k, n) )
      # n__E_T_model_l_m[n].append(E_T_k_n_alt(mu, d, k, n) )
    first_loop = False
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )  
  color = iter(cm.rainbow(numpy.linspace(0, 1, int((N-k)/k_step) ) ) )
  plot.plot(d_l, d_l, label=r'$y=d$', color="black", marker='.', linestyle='', mew=2)
  for n in range(k, N, k_step):
    m = next(marker)
    c = next(color)
    plot.plot(d_l, n__E_T_l_m[n], label=r'$E[T]$, n:{}'.format(n), color=c, marker=m, linestyle='', mew=2)
    # plot.plot(d_l, n__E_T_model_l_m[n], label=r'$E[\hat{T}]$', color=c, alpha=0.7, marker=m, linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$d$')
  plot.ylabel("E[T] (s)")
  plot.title(r'$\mu$= {}, k= {}'.format(mu, k) )
  plot.savefig("plot_autorepeat_k_{}.png".format(k) )
  log(WARNING, "done; k= {}".format(k) )

def plot_prob_N_d():
  k = 10
  n = 15
  d = 0.6 # delta
  def Pr_N_d(r):
    if r > k:
      return 0
    return binomial(k, r)*math.exp(-d*(k-r) )*(1-math.exp(-d) )**r
  r_l, prob_l = [], []
  for r in range(1, k+1):
    r_l.append(r)
    prob_l.append(Pr_N_d(r) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  print("r_l= {}".format(pprint.pformat(r_l) ) )
  plot.plot(r_l, prob_l, color='red', label=r'$Pr\{N_{\delta}=r\}$', marker=next(marker), linestyle='', mew=2)
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'k= {}, n= {}, d= {}'.format(k, n, d) )
  plot.savefig("plot_prob_N_d_k_{}_n_{}_d_{}.png".format(k, n, d) )
  log(WARNING, "done; k= {}, n= {}, d= {}".format(k, n, d) )

def plot_binomial_dist__approx():
  n = 15
  p = 0.4
  def comp_dist(k):
    if k > n:
      return 0
    sum_ = 0
    for i in range(k, n+1):
      sum_ += binomial(n, i) * p**i * (1-p)**(n-i)
    return sum_
  def dist(k):
    if k > n:
      return 0
    sum_ = 0
    for i in range(0, k+1):
      sum_ += binomial(n, i) * p**i * (1-p)**(n-i)
    return sum_
  def chernoff_bound_on_upper_tail(k):
    p_ = k/n
    print("p_= {}".format(p_) )
    return math.exp(n*((p_*math.log(p/p_) ) + (1-p_)*math.log((1-p)/(1-p_) ) ) )
  
  k_l, dist_l, approx_dist_l = [], [], []
  # for k in range(0, n+2):
  for k in range(int(p*n), n):
    k_l.append(k)
    dist_l.append(comp_dist(k) )
    approx_dist_l.append(chernoff_bound_on_upper_tail(k) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  # print("k_l= {}".format(pprint.pformat(k_l) ) )
  plot.plot(k_l, approx_dist_l, color='red', label=r'$Pr\{N \leq k\}_{UB}$', marker=next(marker), linestyle='', mew=2)
  plot.plot(k_l, dist_l, color='black', label=r'$Pr\{N \leq k\}$', marker=next(marker), linestyle='', mew=2)
  plot.xlabel(r'$k$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, p= {}'.format(n, p) )
  plot.savefig("plot_binomial_dist__approx_n_{}.png".format(n) )
  log(WARNING, "done; n= {}, p= {}".format(n, p) )

def plot_den():
  k = 10
  n = 20 # 15
  def exp(r):
    if n-r+1 == 0:
      return None
    return (k-r)*(k-r+1)*(r+1)/(n-r+1)
  r_l, exp_l = [], []
  for r in range(1, k):
    r_l.append(r)
    exp_l.append(exp(r) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(r_l, exp_l, color='red', label='', marker=next(marker), linestyle='', mew=2)
  plot.xlabel(r'$r$')
  plot.ylabel("exp")
  plot.title(r'k= {}, n= {}'.format(k, n) )
  plot.savefig("plot_debinomial_{}_n_{}.png".format(k, n) )
  log(WARNING, "done; k= {}, n= {}".format(k, n) )

if __name__ == "__main__":
  # plot_prob_N_d()
  # plot_den()
  # plot_binomial_dist__approx()
  plot_autorepeat_k_n()
  # plot_autorepeat_dist_k_n()
  # plot_autorepeat_conf_k_n()
