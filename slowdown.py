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

if __name__ == "__main__":
  plot_Pr_slow()
