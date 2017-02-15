import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
from cycler import cycler

from auto_repeat_models import *

# ------------------------  Checks on Pei's work  --------------------- #
def plot_autorepeat_k_n():
  mu = 1
  k = 10
  n = 20
  # d = 5
  
  d_l, E_T_k_n_l, E_T_k_n_alt_l, E_T_k_n_lb_l = [], [], [], []
  E_T_k_n_r_1_5_l, E_T_k_n_r_1_l = [], []
  for d in numpy.arange(0.1, 3*H(k)/mu, 0.1):
    d_l.append(d)
    
    E_T_k_n_l.append(E_T_k_n(mu, d, k, n) )
    E_T_k_n_alt_l.append(E_T_k_n_alt(mu, d, k, n) )
    E_T_k_n_r_1_5_l.append(E_T_k_n(mu, d, k, int(k*1.5) ) )
    E_T_k_n_r_1_l.append(E_T_k_n(mu, d, k, k) )
    
    E_T_k_n_lb_l.append(E_T_k_n_lb(mu, d, k, n) )
  
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(d_l, E_T_k_n_l, color='black', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
  plot.plot(d_l, E_T_k_n_alt_l, color='brown', label=r'$E[T_{alt}]$', marker=next(marker), linestyle='', mew=2)
  # plot.plot(d_l, E_T_k_n_r_1_5_l, color='magenta', label=r'$E[T]$, r:1.5', marker=next(marker), linestyle='', mew=2)
  # plot.plot(d_l, E_T_k_n_r_1_l, color='red', label=r'$E[T]$, r:1', marker=next(marker), linestyle='', mew=2)
  # plot.plot(d_l, E_T_k_n_lb_l, color='green', label=r'$E[\hat{T}_{LB}]$', marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$d$')
  plot.ylabel("E[T] (s)")
  plot.title(r'mu= {}, k= {}, n= {}'.format(mu, k, n) )
  plot.savefig("plot_autorepeat_k_{}_n_{}.png".format(k, n) )
  log(WARNING, "done; k= {}, n= {}".format(k, n) )

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
