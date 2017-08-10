import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

import numpy
import matplotlib.pyplot as plot
from scipy.optimize import fsolve

from patch import *

def E_T_mixednet_ub(n, k, l, qlambda_l=[] ):
  if len(qlambda_l):
    E_T_l = []
    for i,l in enumerate(qlambda_l):
      qlambda_l_ = list(qlambda_l)
      qlambda_l_.remove(l)
      # print("l= {}, qlambda_l_= {}".format(l, qlambda_l_) )
      mu = sum(qlambda_l_[0:k] )
      E_T_l.append(1/(mu-l) )
    log(WARNING, "n= {}, k= {}, qlambda_l= {}\n\t E_T_l= {}".format(n, k, qlambda_l_, E_T_l) )
    return E_T_l
  else:
    E_S = 1/l * (H(n-1) - H(n-k) )
    E_S_2 = 1/l**2 * (H_2(n-1) - H_2(n-k) ) + E_S**2
    E_T = E_S + l*E_S_2/2/(1-l*E_S)
    log(WARNING, "n= {}, k= {}, l= {}\n\t E_T= {}".format(n, k, l, E_T) )
    if E_T < 0: return None
    return E_T

def E_T_mixednet_lb(n, k, l):
  E_S = 1/(n-k+1)/l
  E_S_2 = 2/((n-k+1)*l)**2
  E_T = E_S + l*E_S_2/2/(1-l*E_S)
  log(WARNING, "n= {}, k= {}, l= {}\n\t E_T= {}".format(n, k, l, E_T) )
  if E_T < 0: return None
  return E_T

def E_T_mixednet_approx(n, k, l):
  # Using MC
  # pbusy = Pr_busy_mixednet_approx(n, k)
  # p = 1/(1 + pbusy*(n-k) )
  # mu = p*(n-k+1)*l
  # ro = l/mu
  # ro_1 = (1-p)*ro
  # E_N = ro_1/(1-ro)/(1-ro+ro_1)
  # return E_N/l
  
  # pbusy = (1/(n-k+1) )**(1/k)
  pbusy = (1/(n-k+1) )**(1/(k-1))
  p = pbusy**(k-1)
  print("pbusy= {}, p= {}".format(pbusy, p) )
  # p = pbusy**(k-2)
  mu = p*(n-k+1)*l
  # return (k-1)/n * 1/(mu-l)
  return 1/(mu-l)

def Pr_busy_mixednet_approx(n=100, k=None):
  def eq(pbusy, data):
    k = data
    # print("k= {}".format(k) )
    
    def p():
      sum_ = 0.0
      for i in range(k):
        sum_ += binomial(n, i) * pbusy**i * (1-pbusy)**(n-i)
      return binomial(n, k-1) * pbusy**(k-1) * (1-pbusy)**(n-k+1) / sum_
      
      # sum_ = 0.0
      # for i in range(k-1, n+1):
      #   sum_ += binomial(n, i)* pbusy**i * (1-pbusy)**(n-i)
      # return sum_
      # return binomial(n, k-1)* pbusy**(k-1) * (1-pbusy)**(n-k+1)
    # return p() - 1/(n-k+1)/pbusy
    return p() - 1/(1 + (n-k)*pbusy)
  
  if k is not None:
    root = scipy.optimize.brentq(eq, 0.0001, 0.99, args = (k) )
    print("n= {}, k= {}, root= {}".format(n, k, root) )
    return root
  else:
    mew, ms = 3, 5
    for k in range(1, n+1, 20):
      if k == 1: continue
      # roots = fsolve(eq, 0.0, args=(k,), xtol=1e-06)
      roots = scipy.optimize.brentq(eq, 0.0001, 0.95, args = (k) )
      print("n= {}, k= {}, roots= {}".format(n, k, roots) )
    #   pbusy_l, eq_l = [], []
    #   for pbusy in numpy.linspace(0.01, 1, 1000):
    #     pbusy_l.append(pbusy)
    #     eq_l.append(eq(pbusy, k) )
    #   plot.plot(pbusy_l, eq_l, label=r'$k={}$'.format(k), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.legend()
    # plot.xlabel(r'pbusy', fontsize=13)
    # plot.ylabel(r'Eq', fontsize=13)
    # fig = plot.gcf()
    # # def_size = fig.get_size_inches()
    # # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
    # fig.tight_layout()
    # plot.savefig("prob_busy_complete_eq_n_{}.pdf".format(n) )
    # log(WARNING, "done; n= {}".format(n) )

# def Pr_busy(n, k):
#   return k/n * (1/(n-k+1) )**(1/k)

# *****************************  M/G/1 Approx  ***************************** #
def Pr_X_n_k_leq_x(n, k, l, x): # n:kth order stat
  def F(x):
    return 1 - math.exp(-l*x)
  return sum([binomial(n,i)*F(x)**i*(1-F(x))**(n-i) for i in range(k, n+1) ] )

def serv_tail_approx(p_empty, n, k, l, t):
  # return 1 - (1 - p_empty*math.exp(-l*t) )**(k-1)
  cdf = 0
  for e in range(k):
    cdf += binomial(k-1,e)*p_empty**e*(1-p_empty)**(k-1-e) * Pr_X_n_k_leq_x(n-k+1+e, e, l, t)
  # print("cdf= {}".format(cdf) )
  return 1 - cdf

def plot_serv_tail_approx(n, k, l):
  pe = p_empty_iteratively(n, k, l)
  x_l, y_l = [], []
  for t in numpy.linspace(0, 10, 100):
    x_l.append(t)
    y_l.append(serv_tail_approx(pe, n, k, l, t) )
  plot.plot(x_l, y_l, label=r'$\lambda= {}$'.format(l), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.legend()
  plot.title(r'$n= {}$, $k= {}$, $\lambda= {}$'.format(n, k, l) )
  plot.xlabel(r'$t$', fontsize=13)
  plot.ylabel(r'$Pr\{S > t\}$', fontsize=13)
  plot.savefig("plot_serv_tail_approx_n_{}_k_{}.png".format(n, k) )
  log(WARNING, "done; n= {}, k= {}".format(n, k) )

def serv_moment_approx(p_empty, n, k, l, i):
  # return mpmath.quad(lambda t: i*t**(i-1)*serv_tail_approx(p_empty, n, k, l, t), [0, mpmath.inf] )
  return mpmath.quad(lambda t: i*t**(i-1)*serv_tail_approx(p_empty, n, k, l, t), [0, 100000] )

def E_T_mg1_approx(n, k, l, p_empty):
  E_S = serv_moment_approx(p_empty, n, k, l, 1)
  E_S_2 = serv_moment_approx(p_empty, n, k, l, 2)
  print("n= {}, k= {}, l= {}, E_S= {}, E_S_2= {}".format(n, k, l, E_S, E_S_2) )
  E_T = E_S + l*E_S_2/2/(1-l*E_S)
  if E_T < 0: return None
  return E_T # (k-1)/n*

def p_empty_iteratively(n, k, l):
  pe = 1
  # for i in range(5):
  #   pe = mpmath.quad(lambda t: (1 - serv_tail_approx(pe, n, k, l, t) ) * l*math.exp(-l*t), [0, mpmath.inf] )
  #   # print("i= {}, pe= {}".format(i, pe) )
  
  for k_ in range(1, k+1):
    pe = mpmath.quad(lambda t: (1 - serv_tail_approx(pe, n, k_, l, t) ) * l*math.exp(-l*t), [0, mpmath.inf] )
    print("k_= {}, pe= {}".format(k_, pe) )
  return pe

if __name__ == "__main__":
  # p_empty_iteratively(n=10, k=3, l=1)
  plot_serv_tail_approx(n=10, k=9, l=0.1)
