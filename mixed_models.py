import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot

import scipy
import numpy as np
from scipy.optimize import fsolve

from patch import *
from rvs import *
from commonly_used import *

def lyapunov_stability_test():
  n, k = 100, 10
  
  V = lambda s_l: min(s_l)**(math.log(n/(n-1), 2.0) )
  for c in range(1, 10000):
    s_l = [c]*(k-1)
    
    d = -V(s_l)
    for i in range(k-1):
      s_l_ = list(s_l)
      s_l_[i] += 1
      d += V(s_l_)/n
    s_l_ = [s-1 for s in s_l]
    d += V(s_l_)*(n-k+1)/n
    print("c= {}, d= {}".format(c, d) )

def ET_mixednet_ub(n, k, l, qlambda_l=[] ):
  if len(qlambda_l):
    ET_l = []
    for i,l in enumerate(qlambda_l):
      qlambda_l_ = list(qlambda_l)
      qlambda_l_.remove(l)
      # print("l= {}, qlambda_l_= {}".format(l, qlambda_l_) )
      mu = sum(qlambda_l_[0:k] )
      ET_l.append(1/(mu-l) )
    log(WARNING, "n= {}, k= {}, qlambda_l= {}\n\t ET_l= {}".format(n, k, qlambda_l_, ET_l) )
    return ET_l
  else:
    EV = 1/l * (H(n-1) - H(n-k) )
    EV2 = 1/l**2 * (H_2(n-1) - H_2(n-k) ) + EV**2
    ET = EV + l*EV2/2/(1-l*EV)
    log(WARNING, "n= {}, k= {}, l= {}\n\t ET= {}".format(n, k, l, ET) )
    if ET < 0: return None
    return ET

def ET_mixednet_lb(n, k, l):
  EV = 1/(n-k+1)/l
  EV2 = 2/((n-k+1)*l)**2
  ET = EV + l*EV2/2/(1-l*EV)
  log(WARNING, "n= {}, k= {}, l= {}\n\t ET= {}".format(n, k, l, ET) )
  if ET < 0: return None
  return ET

def ET_mixednet_approx(n, k, l):
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
        sum_ += binom(n, i) * pbusy**i * (1-pbusy)**(n-i)
      return binom(n, k-1) * pbusy**(k-1) * (1-pbusy)**(n-k+1) / sum_
      
      # sum_ = 0.0
      # for i in range(k-1, n+1):
      #   sum_ += binom(n, i)* pbusy**i * (1-pbusy)**(n-i)
      # return sum_
      # return binom(n, k-1)* pbusy**(k-1) * (1-pbusy)**(n-k+1)
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
    #   for pbusy in np.linspace(0.01, 1, 1000):
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
def serv_tail_approx(pe, n, k, t, dist_m):
  cdf = 0
  for e in range(k):
    cdf += binom(k-1,e)*pe**e*(1-pe)**(k-1-e) * Pr_Xnk_leq_x(n-k+1+e, e, t, dist_m)
  return 1 - cdf

def approx_serv_tail_approx(pe, n, k, t, dist_m):
  return 1 - I(F(t, dist_m), (k-1)*pe, n-k+2)

def plot_serv_tail_approx(n, k, dist_m):
  pe = pempty(n, k, l)
  x_l, y_l = [], []
  for t in np.linspace(0, 10, 100):
    x_l.append(t)
    y_l.append(serv_tail_approx(pe, n, k, t, dist_m) )
  plot.plot(x_l, y_l, label=r'$\lambda= {}$'.format(l), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.legend()
  plot.title(r'$n= {}$, $k= {}$, $\lambda= {}$'.format(n, k, l) )
  plot.xlabel(r'$t$', fontsize=13)
  plot.ylabel(r'$Pr\{S > t\}$', fontsize=13)
  plot.savefig("plot_serv_tail_approx_n_{}_k_{}.png".format(n, k) )
  log(WARNING, "done; n= {}, k= {}".format(n, k) )

def serv_moment_approx(pe, n, k, m, dist_m):
  # return mpmath.quad(lambda t: m*t**(m-1)*serv_tail_approx(pe, n, k, t, dist_m), [0, mpmath.inf] ) # 100000
  return scipy.integrate.quad(lambda t: m*t**(m-1)*serv_tail_approx(pe, n, k, t, dist_m), 0, np.inf)[0]

def ET_mg1approx(n, k, dist_m):
  # pe = pempty_approx(n, k)
  pe = pempty(n, k, dist_m)
  EV = serv_moment_approx(pe, n, k, 1, dist_m)
  EV2 = serv_moment_approx(pe, n, k, 2, dist_m)
  print("k= {}, pe= {}".format(k, pe) )
  dist = dist_m['dist']
  if dist == 'Exp':
    ar = dist_m['mu']
    ET = EV + ar*EV2/2/(1-ar*EV)
  elif dist == 'Pareto':
    rv = Pareto(dist_m['loc'], dist_m['a'] )
    EX, VX = rv.mean(), rv.var()
    ar = 1/EX
    coeffvar_ar2 = VX/EX**2
    coeffvar_serv2 = (EV2 - EV**2)/EX**2
    ro = ar*EV
    
    ET = (ro/(1-ro) ) * (coeffvar_ar2 + coeffvar_serv2)/2 * EV
  print("ET= {}".format(ET) )
  if ET < 0: return None
  return ET

def ET_mg1approx_(n, k, dist_m):
  pe = pempty_approx(n, k, dist_m)
  print("k= {}, pe= {}".format(k, pe) )
  ar = dist_m['mu']
  
  EB = serv_moment_approx(pe, n, k, 1, dist_m)
  EB2 = serv_moment_approx(pe, n, k, 2, dist_m)
  ET = EB + ar*EB2/2/(1 - ar*EB)
  return ET

def pempty(n, k, dist_m):
  pe = 1
  
  mu = dist_m['mu']
  x_pdf = lambda x: mu*math.exp(-mu*x)
  for k_ in range(1, k+1):
    # pe = mpmath.quad(lambda t: (1 - serv_tail_approx(pe, n, k_, t, dist_m) ) * f(t, dist_m), [0, mpmath.inf] )
    pe = scipy.integrate.quad(lambda t: (1 - serv_tail_approx(pe, n, k_, t, dist_m) ) * x_pdf(t), 0, np.inf)[0]
    # print("k_= {}, pe= {}".format(k_, pe) )
    
  # for _ in range(10):
  #   print("... pe= {}".format(pe) )
  #   pe = scipy.integrate.quad(lambda t: (1 - serv_tail_approx(pe, n, k, t, dist_m) ) * x_pdf(t), 0, np.inf)[0]
  return pe
  # return 1 - (k-1)/n

def pempty_approx(n, k, dist_m):
  # pe = 1
  # a = (k-1)/(n-k+2)
  # if a == 0: return None
  # return (-1 + math.sqrt(1 + 4*a) )/2/a
  
  ar = dist_m['mu']
  def ro(pe):
    return ar*serv_moment_approx(pe, n, k, 1, dist_m)
  eq = lambda pe: pe - (1 - ro(pe) )
  pe = scipy.optimize.brentq(eq, 0.0001, 1)
  return pe

def plot_qoi():
  dist_m = {'dist': 'Exp', 'mu': 1}
  # dist_m = {'dist': 'Pareto', 'loc': 1, 'a': 50}
  n = 10
  print("n= {}, dist_m= {}".format(n, dist_m) )
  
  x_l, y_l, y_approx_l = [], [], []
  def plot_ar_forfixdelay(d=1):
    ET_base = ET_mg1approx(n, 2, dist_m)
    for k in range(3, n):
      def eq(ar, data):
        k = data
        return ET_mg1approx(n, k, {'dist': 'Exp', 'mu': ar} ) - ET_base
      ar = scipy.optimize.brentq(eq, 0.0001, 100, args = (k) )
      print("ar= {}".format(ar) )
  
  def plot_ED_vs_k():
    for k in range(2, n):
      x_l.append(k)
      ED = ET_mg1approx(n, k, dist_m)
      # E = ED / () if ED is not None else ED
      y_l.append(E)
      
      pe = pempty(n, k, dist_m)
      EV = serv_moment_approx(pe, n, k, 1, dist_m)
      EV2 = serv_moment_approx(pe, n, k, 2, dist_m)
      # y_l.append(EV**2)
      pe = 1
      j, i = n - (k-1)*(1-pe), (k-1)*pe
      # print("k= {}, i= {}, j= {}".format(k, i, j) )
      # y_approx_l.append(math.log(1 + j/(j-i))**2 )
      # y_approx_l.append(math.log(1 + (n+1)/(n-k+1))**2 )
      y_approx_l.append(math.log(math.sqrt((k-1)/(n-k+2) ) )**2)
      # y_approx_l.append(H_cont(j) - H_cont(j-i) )
    plot.plot(x_l, y_l, label='actual', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.plot(x_l, y_approx_l, label='approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.xlabel(r'$k$', fontsize=13)
  
  def plot_avgdelay():
    for k in range(2, n):
    # for k in np.linspace(2, n-1, 10):
    #   k = int(k)
      x_l.append(k)
      y_l.append(ET_mg1approx(n, k, dist_m) )
    plot.plot(x_l, y_l, color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.xlabel(r'$k$', fontsize=13)
    plot.ylabel(r'$E[D]$', fontsize=14)
  
  def plot_pe():
    for k in range(2, n):
      x_l.append(k)
      y_l.append(pempty(n, k, dist_m) )
      y_approx_l.append(pempty_approx(n, k) )
    plot.plot(x_l, y_l, label='Iterative', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.plot(x_l, y_approx_l, label='Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.xlabel(r'$k$', fontsize=13)
    plot.ylabel(r'$p_0$', fontsize=14)
  
  def plot_avgnumempty():
    for k in range(2, n):
      x_l.append(k)
      
      pe = pempty(n, k, dist_m)
      y_l.append((k - 1)*pe)
      pe = pempty_approx(n, k)
      y_approx_l.append((k - 1)*pe)
    plot.plot(x_l, y_l, label='Iterative', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.plot(x_l, y_approx_l, label='Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.xlabel(r'$k$', fontsize=13)
    plot.ylabel(r'$E[N_e]$', fontsize=14)
  
  # plot_ar_forfixdelay()
  plot_ED_vs_k()
  # plot_avgdelay()
  # plot_avgnumempty()
  
  plot.legend()
  plot.title(r'$n= {}$, $X \sim {}$'.format(n, dist_m) )
  plot.savefig("plot_qoi_n_{}.png".format(n) )
  log(WARNING, "done.")

def EL_n_2(n):
  return 1/2/(n-2)

def EL2_n_2(n):
  return n/2/(n-2)**2

def ET_n_2(n, ar):
  p0 = (n-2)/2/(n-1)
  ro = (1-p0)/n
  EL = 1/2/(n-2)
  return 1/(n-1)/ar * (p0 + ro + EL)

def ET2_n_2(n, ar):
  p0 = (n-2)/2/(n-1)
  ro = (1-p0)/n
  EL = 1/2/(n-2)
  EL2 = n/2/(n-2)**2
  return 1/((n-1)*ar)**2 * (2*p0 + 2*ro + EL2 + 3*EL)

def tail_exponent(n, k, dist_m):
  ar = dist_m['mu']
  pe = pempty(n, k, dist_m)
  # pe = pempty_approx(n, k)
  k_ = (k-1)*pe
  n_ = n - (k-1)*(1-pe)
  def eq(s):
    Vs = B(k_, n_-k_+1+s/ar)/B(k_, n_-k_+1)
    return ar + (s - ar)/Vs
  mu = scipy.optimize.brentq(eq, -20, -0.001)
  return mu

if __name__ == "__main__":
  # n = 50
  # dist_m = {'dist': 'Exp', 'mu': 1}
  # print("n= {}, X ~ {}".format(n, dist_m) )
  # for k in range(1, n, 5):
  #   pe = pempty(n, k, dist_m)
  #   print("k= {}, pe= {}".format(k, pe) )
  
  # plot_serv_tail_approx(n=10, k=9, {'dist': 'Exp', 'mu': 1})
  # plot_qoi()
  
  lyapunov_stability_test()
  
  