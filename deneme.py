import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams['ps.useafm'] = True
# matplotlib.rcParams['pdf.use14corefonts'] = True
# matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy

from patch import *

def MG1_tail_for_Pareto_serv(ro, b, a, t):
  def R(x, a):
    s = 0
    for r in range(1, 1000):
      den = 1
      for i in range(1, r+1):
        den *= a - i
      s += x**r/den
    pi = math.pi
    return 1 + s - pi*x**a*math.exp(-x)/G(a)*mpmath.cot(pi*a)
  
  def I(x, a):
    return x**a*math.exp(-x)/G(a)
  
  def H1(x, a, ro):
    return (1 - ro*R(x, a-1) )**2 + (math.pi*ro*I(x, a-1) )**2
  
  return ro*(1-ro)/G(a-1) * \
         mpmath.quad(lambda x: x**(a-2)/H1(x, a, ro) * math.exp(-(1+t/b)*x), [0, mpmath.inf] )

def plot_MG1_tail():
  b, a = 1, 2.11
  
  def plot_(ro):
    t_l, tail_l = [], []
    for t in numpy.logspace(0, 3, 10):
      t_l.append(t)
      tail_l.append(MG1_tail_for_Pareto_serv(ro, b, a, t) )
    plot.plot(t_l, tail_l, label=r'$\rho= {}$'.format(ro), color=next(dark_color), marker=next(marker), linestyle=':')
  
  plot_(ro=0.1)
  plot_(ro=0.3)
  plot_(ro=0.6)
  plot_(ro=0.9)
  
  plot.legend()
  plot.xscale('log')
  plot.yscale('log')
  
  plot.xlabel(r'$t$', fontsize=12)
  plot.ylabel(r'$p(T > t)$', fontsize=12)
  plot.title('$b= {}$, $a= {}$'.format(b, a) )
  
  plot.savefig("plot_MG1_tail.png", bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_ratios_of_Gammas():
  a = 1 # 10 # 1.05 # 3
  k = 10
  # """
  def func(c):
    return G((c+1)*k+1)/G(c*k+1) * G(c*k+1-1/a)/G((c+1)*k+1-1/a) * G(k+1-1/a/(c+1) )/G(k+1)/G(1-1/a/(c+1) )
    # return ((c+1)*k + 1)**(1/a) / (c*k)**(1/a) / k**(1/a/(c+1) ) # / G(1-1/a/(c+1) )
  c_l, f_l = [], []
  for c in range(1, 100):
    c_l.append(c)
    f_l.append(func(c) )
  plot.plot(c_l, f_l, color=next(dark_color), marker=next(marker), linestyle=':')
  plot.xlabel(r'$c$', fontsize=12)
  """
  def func(n, _a=a):
    # return G(n+1)/G(n-k+1) * G(n-k+1-1/a)/G(n+1-1/a) / (1 + k/(a*n-k) )
    return G(n+1)/G(n-k+1) * G(n-k+1-1/_a)/G(n+1-1/_a) - k/(_a*n-k)
    # return (n-k+1)/(n+1)*(n+1-1/a)/(n-k+1-1/a) - a*n/(a*n-k)
  # n_l, f_l = [], []
  # for n in range(k+1, 10*k):
  #   n_l.append(n)
  #   f_l.append(func(n) )
  # plot.plot(n_l, f_l, color=next(dark_color), marker=next(marker), linestyle=':')
  # plot.xlabel(r'$n$', fontsize=12)
  
  # a_l, f_l = [], []
  # for a in numpy.linspace(1, 100, 100):
  #   a_l.append(a)
  #   f_l.append(func(2*k, a) )
  # plot.plot(a_l, f_l, color=next(dark_color), marker=next(marker), linestyle=':')
  # plot.xlabel(r'$\alpha$', fontsize=12)
  
  def func(n, i):
    return (n+1-1/a)/(n-k+1-1/a)*(i+n+1-1/a)/(i+n-k+1-1/a) \
           - (n+1)/(n-k+1)*(i+n+1)/(i+n-k+1)*(a*n)/(a*n-k)
  i_l, f_l = [], []
  for i in range(1, 1000):
    i_l.append(i)
    f_l.append(func(2*k, i) )
  plot.plot(i_l, f_l, color=next(dark_color), marker=next(marker), linestyle=':')
  plot.xlabel(r'$i$', fontsize=12)
  """
  # plot.ylabel(r'', fontsize=12)
  plot.savefig("plot_ratios_of_Gammas.png", bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done.")

if __name__ == "__main__":
  plot_ratios_of_Gammas()
  # plot_MG1_tail()
  