import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot

import pprint, numpy, scipy
from scipy.stats import norm

from patch import *

# *************************************  Binomial Order Stats  *********************************** #
# X ~ Bin(N, p)
def Pr_Bin_X_lep_x(N, p, x):
  return I(1-p, N-x, 1+x)

def EXj_minus_EXi(N, p, n, j, i):
  s = 0
  for x in range(N):
    f = Pr_Bin_X_lep_x(N, p, x)
    s += I(f, i, n-i+1) - I(f, j, n-j+1)
  
  return s

def approx_EXj_minus_EXi(N, p, n, j, i):
  # return math.sqrt(N*p*(1-p)) * (math.sqrt((j-1)/(n-j+1) ) - math.sqrt((i-1)/(n-i+1) ) )
  return math.sqrt(N*p*(1-p)) * (scipy.stats.norm.ppf(j/(n+1) ) - scipy.stats.norm.ppf(i/(n+1) ) )

def plot_EXj_minus_EXi():
  p = 0.5
  n = 100 # 10
  
  N_l = []
  EXn_minus_EX1_l, EXn_minus_EX2_l, EXn_minus_EXnm1_l = [], [], []
  approx_EXn_minus_EX2_l, approx_EXn_minus_EXnm1_l = [], []
  for N in numpy.logspace(2, 6, 10):
    N = int(N)
    N_l.append(N)
    
    # EXn_minus_EX1_l.append(EXj_minus_EXi(N, p, n, n, 1) )
    EXn_minus_EX2_l.append(EXj_minus_EXi(N, p, n, n, 2) )
    approx_EXn_minus_EX2_l.append(approx_EXj_minus_EXi(N, p, n, n, 2) )
    EXn_minus_EXnm1_l.append(EXj_minus_EXi(N, p, n, n, n-1) )
    approx_EXn_minus_EXnm1_l.append(approx_EXj_minus_EXi(N, p, n, n, n-1) )
  # plot.plot(N_l, EXn_minus_EX1_l, label=r'$E[X_n-X_1]$', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(N_l, EXn_minus_EX2_l, label=r'$E[X_n-X_2]$', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  print("approx_EXn_minus_EX2_l= {}".format(pprint.pformat(approx_EXn_minus_EX2_l) ) )
  plot.plot(N_l, approx_EXn_minus_EX2_l, label=r'Approx $E[X_n-X_2]$', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(N_l, EXn_minus_EXnm1_l, label=r'$E[X_n-X_{n-1}]$', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  print("approx_EXn_minus_EXnm1_l= {}".format(pprint.pformat(approx_EXn_minus_EXnm1_l) ) )
  plot.plot(N_l, approx_EXn_minus_EXnm1_l, label=r'Approx $E[X_n-X_{n-1}]$', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
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
  plot_EXj_minus_EXi()
  