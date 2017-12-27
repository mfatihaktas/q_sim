import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot

import pprint, numpy, scipy
from scipy.stats import norm
from scipy.integrate import dblquad

from patch import *

# *************************************  Binomial Order Stats  *********************************** #
# X ~ Bin(N, p)
def Pr_Bin_X_leq_x(N, p, x):
  return I(1-p, N-x, 1+x)

# From "Order statistics arising from independent binomial populations"
def Pr_Xj_minus_Xi(N, p, n, j, i, x):
  # Does not work correctly for some reason
  s = 0
  for k in range(N-x+1):
    if x == 0:
      # int_ = scipy.integrate.dblquad(lambda v, w: w**(i-1)*(v-w)**(j-i-1)*(1-v)**(n-j), \
      #                               Pr_Bin_X_leq_x(N, p, k-1), Pr_Bin_X_leq_x(N, p, k), \
      #                               lambda w: w, lambda w: Pr_Bin_X_leq_x(N, p, k) )
      pass
    else:
      # int_ = scipy.integrate.dblquad(lambda v, w: w**(i-1)*(v-w)**(j-i-1)*(1-v)**(n-j), \
      #                               Pr_Bin_X_leq_x(N, p, k-1), Pr_Bin_X_leq_x(N, p, k), \
      #                               lambda w: Pr_Bin_X_leq_x(N, p, k+x-1), lambda w: Pr_Bin_X_leq_x(N, p, k+x) )
      def bounds_w():
        return [Pr_Bin_X_leq_x(N, p, k-1), Pr_Bin_X_leq_x(N, p, k) ]
      def bounds_v(w):
        return [Pr_Bin_X_leq_x(N, p, k+x-1), Pr_Bin_X_leq_x(N, p, k+x) ]
      int_ = scipy.integrate.nquad(lambda v, w: w**(i-1)*(v-w)**(j-i-1)*(1-v)**(n-j), [bounds_v, bounds_w] )
    s += G(n+1)/G(i)/G(j-i)/G(n-j+1) * int_[0]
  return s

def Pr_Xn_minus_X1(N, p, n, x):
  s = 0
  if x == 0:
    for k in range(N+1):
      s += (Pr_Bin_X_leq_x(N, p, k) - Pr_Bin_X_leq_x(N, p, k-1) )**n
  else:
    for k in range(N-x+1):
      s += (Pr_Bin_X_leq_x(N, p, k+x) - Pr_Bin_X_leq_x(N, p, k-1) )**n \
         - (Pr_Bin_X_leq_x(N, p, k+x) - Pr_Bin_X_leq_x(N, p, k) )**n \
         - (Pr_Bin_X_leq_x(N, p, k+x-1) - Pr_Bin_X_leq_x(N, p, k-1) )**n \
         + (Pr_Bin_X_leq_x(N, p, k+x-1) - Pr_Bin_X_leq_x(N, p, k) )**n
  return s
  
def plot_Pr_Xj_minus_Xi():
  N = 100
  p = 0.5
  n = 10
  
  # def plot_(j, i):
  #   x_l, Pr_l = [], []
    
  #   for x in range(1, 10):
  #     x_l.append(x)
  #     Pr_l.append(Pr_Xj_minus_Xi(N, p, n, j, i, x) )
  #   plot.plot(x_l, Pr_l, label=r'$j= {}$, $i= {}$'.format(j, i), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot_(10, 6)
  
  x_l, Pr_l = [], []
  for x in range(1, 10):
    x_l.append(x)
    Pr_l.append(Pr_Xn_minus_X1(N, p, n, x) )
  plot.plot(x_l, Pr_l, label=r'$j= {}$, $i= {}$'.format(n, 1), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend()
  # plot.xscale('log')
  # plot.yscale('log')
  plot.xlabel(r'$x$', fontsize=13)
  plot.ylabel(r'$p(X_j - X_i = x)$', fontsize=13)
  plot.title(r'$N= {}$, $p= {}$, $n= {}$'.format(N, p, n) )
  plot.savefig("plot_Pr_Xj_minus_Xi_n_{}.png".format(n) )
  plot.gcf().clear()
  log(WARNING, "done; n= {}".format(n) )

def EXj_minus_EXi(N, p, n, j, i):
  s = 0
  for x in range(N):
    f = Pr_Bin_X_leq_x(N, p, x)
    s += I(f, i, n-i+1) - I(f, j, n-j+1)
  return s

def approx_EXj_minus_EXi(N, p, n, j, i):
  # return math.sqrt(N*p*(1-p)) * (math.sqrt((j-1)/(n-j+1) ) - math.sqrt((i-1)/(n-i+1) ) )
  return math.sqrt(N*p*(1-p)) * (scipy.stats.norm.ppf(j/(n+1) ) - scipy.stats.norm.ppf(i/(n+1) ) )

def plot_EXj_minus_EXi():
  p = 0.5
  n = 10
  
  def plot_(j, i):
    N_l, EXj_minus_EXi_l, approx_EXj_minus_EXi_l = [], [], []
    
    # for N in numpy.logspace(2, 6, 10):
    for N in numpy.logspace(2, 3, 10):
      N = int(N)
      N_l.append(N)
      
      EXj_minus_EXi_l.append(EXj_minus_EXi(N, p, n, j, i) )
      approx_EXj_minus_EXi_l.append(approx_EXj_minus_EXi(N, p, n, j, i) )
    label = r'$E[X_{}-X_{}]$'.format("{}".format(j), "{}".format(i) )
    plot.plot(N_l, EXj_minus_EXi_l, label=label, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(N_l, approx_EXj_minus_EXi_l, label=r'Approx {}'.format(label), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot_(j=n, i=6)
  plot_(j=n, i=8)
  
  plot.legend()
  plot.xscale('log')
  plot.yscale('log')
  plot.xlabel(r'$N$', fontsize=13)
  plot.ylabel(r'', fontsize=13)
  plot.title(r'$n= {}$, $p= {}$'.format(n, p) )
  plot.savefig("plot_EXj_minus_EXi_n_{}.png".format(n) )
  plot.gcf().clear()
  log(WARNING, "done.")

if __name__ == "__main__":
  # plot_EXj_minus_EXi()
  plot_Pr_Xj_minus_Xi()
  