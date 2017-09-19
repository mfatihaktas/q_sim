import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy
from patch import *
from rvs import *

"""
  Small and large files arrive with rate l_s, l_l
  Small f size; S ~ U[S_l, S_u]
  Large f size; L ~ U[L_l, L_u]
  Sching:
  divide-to-all: All f's are divided across k servers
  divide-to-some: Small f's are divided across m while large f's divided across k-m
"""
def plot_sching_small_large_fs():
  pass

"""
  k servers divided into two groups of m and k-m servers
  1st group stores files of size < M while second stores of size > M
  Files are distributed uniformly across all servers within each group.
  File sizes ~ Pareto(l, g, a)
"""
E_T_MAX = 100
def plot_distribute_Pareto_w_threshold():
  k = 100
  l, u, a = 1, 10000, 1.5
  
  F = TPareto(l, u, a)
  def E_T(m, M, ar):
    if m == 0:
      return PK(1/k*F.moment(1), 1/k**2*F.moment(2), ar)
    
    def S_moment(i):
      def S_cdf(x):
        if x < l: return 0
        elif x >= M: return 1
        else: return F.cdf(x)/F.cdf(M)
      return 1/m**i * mpmath.quad(lambda x: i*x**(i-1) * (1 - S_cdf(x) ), [0, mpmath.inf] )
    
    def L_moment(i):
      def L_tail(x):
        if x < M: return 1
        elif x >= u: return 0
        else: return F.tail(x)/F.tail(M)
      return 1/(k-m)**i * mpmath.quad(lambda x: i*x**(i-1) * L_tail(x), [0, mpmath.inf] )
    
    p_S = F.cdf(M)
    p_L = 1 - p_S
    print("m= {}, M= {}: p_S= {}".format(m, M, p_S) )
    E_T_S = PK(S_moment(1), S_moment(2), ar*p_S)
    E_T_L = PK(L_moment(1), L_moment(2), ar*p_L)
    if E_T_S is None or E_T_L is None:
      return None
    return p_S*E_T_S + p_L*E_T_L
  
  ar_ub = 1/(F.moment(1)/k)
  def plot_(m, M=0):
    ar_l, E_T_l = [], []
    for ar in numpy.linspace(0.05, ar_ub, 20):
      ar_l.append(ar)
      E_T_l.append(E_T(m, M, ar) )
    plot.plot(ar_l, E_T_l, label=r'$m= {}, M= {}$'.format(m, M), color=next(dark_color), marker=next(marker) )
  
  plot_(m=0)
  # for m in range(1, 10, 2):
  for m in [99]:
    # for M in numpy.linspace(u/10, u/1.5, 3):
    for M in [*numpy.linspace(1.05, 1.5, 10), u/2]:
      plot_(m, M)
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel(r'$E[T]$')
  plot.title(r'Pareto(l={}, u={}, a={}), k= {}'.format(l, u, a, k) )
  plot.savefig("plot_distribute_Pareto_w_threshold_k_{}.png".format(k) )
  
if __name__ == "__main__":
  plot_distribute_Pareto_w_threshold()
  