import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy
from patch import *

"""
  S = a*X
  File size: X ~ Pareto
  Slowdown: a ~ Dolly -- for tractable analysis --> Bern(p_s) in {L, U}
"""
def Pr_slow_simplex(t, p_s):
  return p_s*(1-(1-p_s)**2)**t

def Pr_slow_mds(n, k, p_s):
  return I(p_s, n-k+1, k)

def plot_Pr_slow():
  def plot_(k):
    t = 2**(k-1) - 1
    n = 2**k - 1
    
    p_s_l = []
    pr_slow_simplex_l, pr_slow_mds_l = [], []
    for p_s in numpy.linspace(0, 1, 20):
      p_s_l.append(p_s)
      pr_slow_simplex_l.append(Pr_slow_simplex(t, p_s) )
      pr_slow_mds_l.append(Pr_slow_mds(n, k, p_s) )
    plot.plot(p_s_l, pr_slow_simplex_l, label=r'Simplex, $k= {}$, $t= {}$'.format(k, t), color=next(dark_color), marker=next(marker) )
    plot.plot(p_s_l, pr_slow_mds_l, label=r'MDS, $k= {}$, $n= {}$'.format(k, n), color=next(dark_color), marker=next(marker) )
  
  # plot_(k=2)
  # plot_(k=3)
  plot_(k=4)
  plot_(k=5)
  plot_(k=6)
  
  plot.legend()
  plot.xlabel(r'$p_s$')
  plot.ylabel(r'Probability of Slowdown')
  plot.savefig("plot_Pr_slow.png")

"""
  A batch of k tasks executed with Qing or PS on a single node
  With PS, X' ~ X/k
"""
def plot_E_T():
  # k tasks served at a q with server ~ Pareto(loc, a)
  def E_T_w_q(k, loc, a):
    return k*loc*a/(a-1)
  # k tasks served at a server ~ Pareto(loc, a)
  def E_T_w_ps(k, loc, a):
    return k*loc*G(k+1)*G(1-1/a)/G(k+1-1/a)
  # k tasks served at a server ~ Pareto(loc, a) by m tasks at a time
  def E_T_w_q_ps(k, loc, a, m):
    if k % m:
      log(ERROR, "k % m= {}".format(k % m) )
      return 1
    return k/m * m*loc*G(m+1)*G(1-1/a)/G(m+1-1/a)
  
  loc = 3
  a = 2
  k = 100
  
  x_l, E_T_w_q_l, E_T_w_ps_l = [], [], []
  # title = r'$\lambda= {}$, $\alpha= {}$'.format(loc, a)
  # xlabel = r'$k$'
  # for k in range(1, 20):
  #   x_l.append(k)
  
  # title = r'$k= {}$, $\alpha= {}$'.format(k, a)
  # xlabel = r'$\lambda$'
  # for loc in numpy.linspace(0.5, 10, 50):
  #   x_l.append(loc)
  
  title = r'$k= {}$, $\lambda= {}$'.format(k, loc)
  xlabel = r'$\alpha$'
  for a in numpy.linspace(2, 10, 50):
    x_l.append(a)
    E_T_w_q_l.append(E_T_w_q(k, loc, a) )
    E_T_w_ps_l.append(E_T_w_ps(k, loc, a) )
  plot.plot(x_l, E_T_w_q_l, label=r'q', color=next(dark_color), marker=next(marker) )
  plot.plot(x_l, E_T_w_ps_l, label=r'ps', color=next(dark_color), marker=next(marker) )
  
  def plot_E_T_w_q_ps(m):
    x_l, y_l = [], []
    for a in numpy.linspace(2, 10, 50):
      x_l.append(a)
      y_l.append(E_T_w_q_ps(k, loc, a, m) )
    plot.plot(x_l, y_l, label=r'q-ps, m= {}'.format(m), color=next(dark_color), marker=next(marker) )
  
  plot_E_T_w_q_ps(m=2)
  plot_E_T_w_q_ps(m=5)
  
  plot.legend()
  plot.title(title)
  plot.xlabel(xlabel)
  plot.ylabel(r'$E[T]$')
  plot.savefig("plot_E_T.png")

if __name__ == "__main__":
  # plot_Pr_slow()
  plot_E_T()

  # plot_Pr_slow()
  plot_E_T()
