import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import pprint, math, numpy

from commonly_used import *

def ar_ub(n, r, k, dist_m):
  EV = EXm_n_k(1, r, k, dist_m)
  return float(1/EV)*n/r

def ET(ar, n, r, k, dist_m):
  ar = ar/(n/r)
  EV = EXm_n_k(1, r, k, dist_m)
  if ar*EV >= 1: return None
  EV2 = EXm_n_k(2, r, k, dist_m)
  ET = EV + ar/2 * EV2/(1 - ar*EV)
  
  if ET < 0 or ET >= 100: return None
  return ET

def plot_split_r1():
  n = 20
  # dist_m = {'dist': 'Exp', 'mu': 1}
  # dist_m = {'dist': 'SExp', 'D': 0.1, 'mu': 1}
  dist_m = {'dist': 'Pareto', 'l': 2, 'a': 2}
  log(WARNING, "n= {}, dist_m={}".format(n, dist_m) )
  
  def plot_(r, k):
    dist_m_ = scale_dist(dist_m, k)
    ub = ar_ub(n, r, k, dist_m_)
    print("r= {}, ar_ub= {}".format(r, ub) )
    
    x_l, y_l = [], []
    for ar in numpy.linspace(0.05, ub, 20):
      x_l.append(ar)
      y_l.append(ET(ar, n, r, k, dist_m_) )
    plot.plot(x_l, y_l, color=next(dark_color), label='r= {}, k= {}'.format(r, k), marker=next(marker), linestyle=':', mew=2)
  
  # plot_(r=1, k=1)
  # plot_(r=2, k=2)
  # plot_(r=4, k=4)
  # plot_(r=5, k=5)
  # plot_(r=10, k=10)
  # plot_(r=20, k=20)
  
  plot_(r=10, k=10)
  plot_(r=10, k=5)
  plot_(r=10, k=2)
  plot_(r=10, k=1)
  
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel(r'$E[T]$')
  plot.title(r'$n= {}$, $V \sim {}$'.format(n, dist_m) )
  plot.savefig("plot_split_r1_n_{}.png".format(n) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_split_r1()