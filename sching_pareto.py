import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot

from patch import *

def E_T(D, k, a):
  return D/k * G(k+1)/G(k+1-1/a)*G(1-1/a)
  # return D * G(k+1)/G(k+1-1/a)*G(1-1/a)

def plot_E_T():
  D = 1
  a = 2
  
  def plot_(a):
    k_l, E_T_l = [], []
    for k in range(1, 31):
      k_l.append(k)
      E_T_l.append(E_T(D, k, a) )
    plot.plot(k_l, E_T_l, label=r'$\alpha={}$'.format(a), color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  plot_(a=1.4)
  plot_(a=1.5)
  plot_(a=2)
  plot_(a=2.5)
  plot_(a=5)
  plot_(a=25)
  
  plot.legend()
  plot.xlabel(r'$k$')
  plot.ylabel(r'$E[T]$')
  plot.title(r'$X \sim Pareto(D={}/k, \alpha)$'.format(D) )
  plot.savefig("plot_E_T.png" )
  plot.gcf().clear()
  log(WARNING, "done.")

if __name__ == "__main__":
  plot_E_T()
