import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams['ps.useafm'] = True
# matplotlib.rcParams['pdf.use14corefonts'] = True
# matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy

from rvs import *
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

def reptoall(t):
  def compare_fl_festl(ro_l):
    f_l = [0]*(t+1)
    f_l[0] = 1/(1 + sum([numpy.prod(ro_l[:j+1] ) for j in range(1, t+1) ] ) )
    f_l[1:] = [f_l[0]*numpy.prod(ro_l[:j+1] ) for j in range(1, t+1) ]
    print("sum(f_l)= {}, f_l=\n  {}".format(sum(f_l), f_l) )
    
    ro = max(ro_l[1:] )
    fest_l = [0]*(t+1)
    fest_l[0] = (1 - ro)/(1 - ro**(t+1) )
    fest_l[1:] = [fest_l[0]*ro**j for j in range(1, t+1) ]
    print("sum(fest_l)= {}, fest_l=\n  {}".format(sum(f_l), f_l) )
    
    print(">>> fest_l - f_l= {}".format(list(numpy.array(fest_l) - numpy.array(f_l) ) ) )
  
  for _ in range(10):
    l = numpy.random.randint(100, size=t)
    ro_l = [1, *(l/sum(l) ) ]
    print("ro_l= {}".format(ro_l) )
    compare_fl_festl(ro_l)

def rep_wcancel():
  def Pr_Tgt(V, X, c, t):
    Pr_V_l_X = mpmath.quad(lambda x: V.cdf(x)*X.pdf(x), [0, 10000*10] )
    return V.tail(t) * (V.tail(t)*Pr_V_l_X + 1 - Pr_V_l_X)**c
  
  def momenti(i, V, X, c):
    return mpmath.quad(lambda t: i*t**(i-1) * Pr_Tgt(V, X, c, t), [0, 10000*10] ) # mpmath.inf
  
  def steadystate_Pr_Tgt(V, X, c, t):
    def Pr_Tgt_(Pr_Vgt, X, c, t):
      return  Pr_Vgt * (Pr_Vgt*X.tail(t) + X.cdf(t) )**c
    tail_pr = V.tail(t)
    for i in range(20):
      print("i= {}, tail_pr= {}".format(i, tail_pr) )
      tail_pr = Pr_Tgt_(tail_pr, X, c, t)
  
  V = Pareto(1, 2)
  def plot_tail(c, ar):
    X = Exp(ar)
    
    t_l, Pr_Tgt_l = [], []
    for t in numpy.linspace(0.05, 20, 100):
      t_l.append(t)
      Pr_Tgt_l.append(Pr_Tgt(V, X, c, t) )
    plot.plot(t_l, Pr_Tgt_l, label=r'$c= {}$, $\lambda= {}$'.format(c, ar), color=next(dark_color), linestyle='-')
    plot.xlabel(r'$t$', fontsize=14)
    plot.ylabel(r'$Pr\{T > t\}$', fontsize=14)
  
  def plot_ETi(c, i):
    ar_l, ETi_l = [], []
    for ar in numpy.logspace(-10, 0.5, 50):
      X = Exp(ar)
      ar_l.append(ar)
      ETi_l.append(momenti(i, V, X, c) )
    plot.plot(ar_l, ETi_l, label=r'$c= {}$'.format(c), color=next(dark_color), linestyle='-')
    plot.xlabel(r'$\lambda$', fontsize=14)
    plot.ylabel(r'$E[T^{}]$'.format(i), fontsize=14)
    
  c = 2
  plot_tail(c, ar=0.1)
  plot_tail(c, ar=1)
  plot_tail(c, ar=10)
  plot_tail(c, ar=100)
  
  # X = Exp(0.1)
  # steadystate_Pr_Tgt(V, X, c=1, t=2)
  
  # i = 3
  # plot_ETi(c=1, i=i)
  # plot_ETi(c=2, i=i)
  # plot_ETi(c=3, i=i)
  
  plot.legend()
  plot.title(r'$V \sim {}$, $X \sim Exp(\lambda)$'.format(V) )
  plot.savefig("rep_wcancel.png", bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done.")

def waitingtime_repwcancel():
  def laplace(X, r):
    return mpmath.quad(lambda x: math.exp(-r*x) * X.pdf(x), [0, mpmath.inf] ) # 10000*10
  
  # V = Exp(1)
  # V = Pareto(1, 3)
  V = Pareto(0.1, 3)
  V21 = X_n_k(V, 2, 1)
  EV = moment_ith(1, V)
  EV21 = moment_ith(1, V21)
  def solvefor_War(ar):
    X = Exp(ar)
    V_21 = X_n_k(V, 2, 1)
    a = laplace(V_21, ar)
    b = laplace(V, ar)
    
    ro = ar*EV
    eq = lambda W: (a-b)*W**2 + (b + ar*(EV21 - EV) )*W + ar*EV-1
    roots = scipy.optimize.brentq(eq, 0.0001, 2)
    print("ar= {}, roots= {}".format(ar, roots) )
  
  for ar in numpy.linspace(0.05, 1/EV-0.05, 10):
    solvefor_War(ar)

if __name__ == "__main__":
  # plot_ratios_of_Gammas()
  # plot_MG1_tail()
  # reptoall(t=5)
  # rep_wcancel()
  
  waitingtime_repwcancel()
  