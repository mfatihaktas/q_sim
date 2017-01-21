import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy, simpy, getopt, itertools
from math import factorial
from numpy import linalg
from patch import *

def E_S_replicate_to_n(f, D, n):
  q = 1-numpy.e**(-f*D)
  return 1/(1-q**n) * D

def E_S_dist_to_n(f, D, n):
  q = 1-numpy.e**(-f*D)
  return (-1/numpy.log(q)*harmonic_sum(n) + 0.5) * D/n

def plot_det():
  f = 0.5 # 0.1
  D = 1
  
  E_S_rep_to_n_l, E_S_dist_to_n_l = [], []
  n_l = []
  for n in range(1, 10):
    n_l.append(n)
    
    E_S_rep_to_n_l.append(E_S_replicate_to_n(f, D, n) )
    E_S_dist_to_n_l.append(E_S_dist_to_n(f, D, n) )
  color = iter(cm.rainbow(numpy.linspace(0, 10, 10) ) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(n_l, E_S_rep_to_n_l, label=r'$E[S], rep-to-n$', color=next(color), marker=next(marker), linestyle='', mew=2)
  plot.plot(n_l, E_S_dist_to_n_l, label=r'$E[S], dist-to-n$', color=next(color), marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$n$')
  plot.ylabel("E[T] (s)")
  plot.title(r'$f$= {}, $D$= {}'.format(f, D) )
  plot.savefig("plot_det.png")
  log(WARNING, "done; f= {}, D= {}".format(f, D) )

if __name__ == "__main__":
  plot_det()