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
from arepeat_models import a_wred_

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
  
  plot.savefig("plot_MG1_tail.pdf", bbox_inches='tight')
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
  plot.savefig("plot_ratios_of_Gammas.pdf", bbox_inches='tight')
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
  plot.savefig("rep_wcancel.pdf", bbox_inches='tight')
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

def EC_k_n():
  def EC(n, k, l, a):
    return l*n/(a-1) * (a - G(n)/G(n-k)*G(n-k+1-1/a)/G(n+1-1/a) )
  
  def ECkp1_m_ECk(k, l, a):
    return a - G(k+2)*G(2-1/a)/G(k+2-1/a) + G(k+1)*G(1-1/a)/G(k+1-1/a)
  
  l, a = 1, 1.5
  # k = 10
  # EC0 = EC(k, k, l, a)
  # for n in range(k, 3*k):
  #   EC_m_EC0 = EC(n, k, l, a) - EC0
  #   print("n-k= {}, EC_m_EC0= {}".format(n-k, EC_m_EC0) )
  
  print("G(0)= {}".format(G(0) ) )
  
  for k in range(10, 20):
    c = ECkp1_m_ECk(k, l, a)
    print("k= {}, ECkp1_m_ECk= {}".format(k, c) )

def plots_for_proposalslides():
  # """
  def plot_(S, color, label=None):
    s_l = numpy.linspace(1, 40, 1000) # 1, 10
    Pr_S_g_s_l = [S.tail(s) for s in s_l]
    label = 'Mean= {}'.format(S.mean() ) if label is None else label
    plot.plot(s_l, Pr_S_g_s_l, marker='None', label=label, color=color, ls='-', lw=2)
  plot.xscale('log')
  plot.yscale('log')
  
  '''
  plot_(S=Pareto(1, 2), color='black')
  plot_(S=Pareto(1.5, 4), color='orange')
  
  plot.legend(loc='lower right', fontsize=14, framealpha=0.25)
  plot.xscale('log')
  plot.yscale('log')
  
  # plot.title(r'${}$'.format(task_t_in_latex), fontsize=24)
  plot.xlabel(r'Execution time $x$', fontsize=20)
  plot.ylabel(r'$\Pr\{X > x\}$', fontsize=24)
  prettify(plot.gca() )
  plot.gcf().set_size_inches(3, 4)
  plot.savefig("plot_low_vs_high_variability.pdf", bbox_inches='tight') # , dpi=fig.dpi
  '''
  # '''
  plot.ylim(top=1.5)
  plot_(S=Exp(1, 1), color='red', label='1 + Exp(1)') # r'Exp: $e^{-\mu x}$'
  plot_(S=Pareto(1, 2), color='black', label='Pareto(1, 2)') # r'Pareto: $x^{-\alpha}$'
  
  # plot_(S=Pareto(1, 3), color='black', label=r'$\alpha_i=3$')
  # plot_(S=Pareto(1, 2), color='orange', label=r'$\alpha_j=2$')
  
  plot.legend(loc='lower left', fontsize=14, framealpha=0.25)
  plot.xlabel(r'$s$', fontsize=20)
  plot.ylabel(r'$\Pr\{Sl > s\}$', fontsize=20)
  prettify(plot.gca() )
  # plot.gcf().set_size_inches(3, 4)
  plot.gcf().set_size_inches(4, 3)
  plot.savefig("plot_tail_exp_vs_pareto.pdf", bbox_inches='tight')
  # '''
  """
  r_l = []
  suffcond_rep_l, suffcond_coding_l = [], []
  # for r in range(2, 20):
  # for r in np.linspace(1.5, 20, 100):
  r = 1.4
  while r < 6:
    r_l.append(r)
    r_ = r + 0.1
    suffcond_rep_l.append(r_/r)
    suffcond_coding_l.append(math.log(r/(r-1)) / math.log(r_/(r_-1)) )
    r = r_
  
  plot.plot(r_l, suffcond_rep_l, label='Replication', color='green', marker='o', ls=':', lw=2)
  plot.plot(r_l, suffcond_coding_l, label='Coding', color='red', marker='x', ls=':', lw=2)
  
  plot.legend(loc='upper right', fontsize=14, framealpha=0.25)
  plot.xlabel(r'$r$', fontsize=20)
  plot.ylabel(r'$\alpha_i/\alpha_j < $', fontsize=20)
  prettify(plot.gca() )
  plot.gcf().set_size_inches(4, 3)
  plot.savefig("plot_suffcond_on_tailindexratio.pdf", bbox_inches='tight')
  """
  """
  a_0 = 4
  r_l, a_l = [], []
  for r in numpy.linspace(1, 10, 100):
    r_l.append(r)
    a_l.append(a_wred_(a_0, r) )
  plot.plot(r_l, a_l, color='green', ls='-', lw=2)
  
  plot.ylim(top=a_0+0.5)
  plot.legend(loc='upper right', fontsize=14, framealpha=0.25)
  plot.xlabel(r'$r$', fontsize=20)
  plot.ylabel(r'$\alpha$', fontsize=20)
  prettify(plot.gca() )
  plot.gcf().set_size_inches(4, 3)
  plot.savefig("plot_tailindex_vs_r.pdf", bbox_inches='tight')
  """
  plot.gcf().clear()
  log(WARNING, "done.")

if __name__ == "__main__":
  # plot_ratios_of_Gammas()
  # plot_MG1_tail()
  # reptoall(t=5)
  # rep_wcancel()
  
  # waitingtime_repwcancel()
  # EC_k_n()
  
  plots_for_proposalslides()
