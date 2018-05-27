import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot

import numpy as np

from patch import *
from rvs import *
from commonly_used import *

def ET_randmix(n, k, ar, pd=None):
  if pd is None:
    pd = 1/n
  po = pd*(k-1)/(n-1) + (1 - pd)*k/(n-1)
  ro = (1-pd)/(n-1)/po
  """
  mu = ar*(pd + (n-1)*po)
  EBe, EBe2 = (1 - pd)*1/mu, (1 - pd)*2/mu**2
  EB, EB2 = 1/mu, 2/mu**2
  
  ET = EBe/(1 - ar*(EB - EBe) ) + ar*EB2/2/(1 - ar*EB) + ar*(EBe2 - EB2)/2/(1 - ar*(EB - EBe) )
  return ET
  """
  return ro/ar/(1 - ro)

def plot_Pr_randmix_attack_sim_vs_model():
  n = 20
  N = 100
  def sim(k, num_samples=10*1000):
    ro = (n - k)/k/(n - 1)
    
    m_l = []
    for _ in range(num_samples):
      M0_l, Mothers = [], np.zeros((N, n-1) )
      for r in range(N):
        if np.random.uniform(0, 1) <= k/n:
          to0, numto_others = 1, k-1
        else:
          to0, numto_others = 0, k
        M0_l.append(to0)
        for i in np.random.choice(n-1, numto_others, replace=False):
          if np.random.uniform(0, 1) <= ro:
            Mothers[r, i] = 1
        
        m_l.append(sum(M0_l) - np.amax(np.sum(Mothers, axis=0) ) )
    return m_l
  
  def plot_(k):
    m_l, Pr_l = [], []
    for m in range(N):
      m_l.append(m)
      Pr_l.append(Pr_randmix_attack_is_mcertain(n, k, N, m) )
    plot.plot(m_l, Pr_l, label=r'$k= {}$, upper-bound'.format(k), color=next(dark_color) )
    
    m_l = sim(k)
    m_l = np.sort(m_l)
    x_l = m_l[::-1]
    y_l = np.arange(m_l.size)/m_l.size
    plot.plot(x_l, y_l, label=r'$k= {}$, sim'.format(k), color=next(dark_color) )
  
  plot_(k=2)
  plot_(k=5)
  plot_(k=10)
  plot_(k=15)
  plot_(k=19)
  
  plot.xlabel(r'$m$', fontsize=13)
  plot.ylabel(r'$Pr\{M_0 - \max(M_1, \dots, M_{n-1}) \geq m\}$', fontsize=13)
  plot.title(r'$n= {}$, $N= {}$'.format(n, N) )
  plot.legend()
  plot.gcf().tight_layout()
  plot.savefig("plot_sim_vs_model_n{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

def Pr_randmix_attack_is_mcertain(n, k, N, m):
  # Pr{M(0) - max{M(1), ..., M(n-1) } >= m} in N observations
  M0 = Binomial(N, k/n)
  ro = (n - k)/k/(n - 1)
  Mi = Binomial(N, ro*(k/n*(k-1)/(n-1) + (1 - k/n)*k/(n-1) ) )
  return sum([M0.pdf(w)*(Mi.cdf(w-m) )**(n-1) for w in range(m, N+1) ] )

def plot_attack():
  n = 20
  N = 100
  def plot_(k):
    m_l, Pr_l = [], []
    for m in range(N):
      m_l.append(m)
      Pr_l.append(Pr_randmix_attack_is_mcertain(n, k, N, m) )
    plot.plot(m_l, Pr_l, label=r'$k= {}$'.format(k), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  # plot_(k=2)
  # plot_(k=5)
  for k in range(2, n, 5):
    plot_(k)
  # plot_(n)
  
  plot.legend()
  plot.xlabel(r'$m$', fontsize=13)
  plot.ylabel(r'$Pr\{M_0 - \max(M_1, \dots, M_{n-1}) \geq m\}$', fontsize=13)
  plot.title(r'$n= {}$, $N= {}$'.format(n, N) )
  plot.gcf().tight_layout()
  plot.savefig("plot_attack_n{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

def optimal_pd():
  n = 20
  
  def po(k, pd): return pd*(k-1)/(n-1) + (1 - pd)*k/(n-1)
  def ro(k, pd): return (1-pd)/(n-1)/po(k, pd)
  def pd(k):
    eq = lambda pd: pd - (1-pd)*k/(n-1)*ro(k, pd)
    pd = scipy.optimize.brentq(eq, 0.00001, 1)
    print("k= {}, pd= {}".format(k, pd) )
    return pd
  
  N = 100
  def sim(k, pd, num_samples=10*1000):
    ro_ = ro(k, pd)
    
    m_l = []
    for _ in range(num_samples):
      M0_l, Mothers = [], np.zeros((N, n-1) )
      for r in range(N):
        if np.random.uniform(0, 1) <= pd:
          to0, numto_others = 1, k-1
        else:
          to0, numto_others = 0, k
        M0_l.append(to0)
        for i in np.random.choice(n-1, numto_others, replace=False):
          if np.random.uniform(0, 1) <= ro_:
            Mothers[r, i] = 1
        
        m_l.append(sum(M0_l) - np.amax(np.sum(Mothers, axis=0) ) )
    return m_l
  
  def plot_(k, pd):
    m_l = sim(k, pd)
    m_l = np.sort(m_l)
    x_l = m_l[::-1]
    y_l = np.arange(m_l.size)/m_l.size
    plot.plot(x_l, y_l, label=r'$k= {}$, $p_d= {}$'.format(k, pd), color=next(dark_color) )
  
  # for k in range(2, n):
  #   pd(k)
  
  k = n-1
  plot_(k, k/n)
  plot_(k, pd(k) )
  plot.xlabel(r'$m$', fontsize=13)
  plot.ylabel(r'$Pr\{M_0 - \max(M_1, \dots, M_{n-1}) \geq m\}$', fontsize=13)
  plot.title(r'$n= {}$, $N= {}$'.format(n, N) )
  plot.legend()
  plot.gcf().tight_layout()
  plot.savefig("optimal_pd_n{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

def plot_load():
  n = 20
  # pd: Pr{releasing the q that received the arrival}
  def load(k, pd):
    po = pd*(k-1)/(n-1) + (1 - pd)*k/(n-1)
    return (1-pd)/(n-1)/po
  
  def plot_(k):
    pd_l, load_l = [], []
    for pd in np.linspace(0, 1, 20):
      pd_l.append(pd)
      load_l.append(load(k, pd) )
    plot.plot(pd_l, load_l, label=r'$k= {}$'.format(k), color=next(dark_color) )
  
  plot_(k=2)
  plot_(k=5)
  
  plot.legend()
  plot.xlabel(r'$p_d$', fontsize=13)
  plot.ylabel(r'$\rho$', fontsize=13)
  plot.title(r'$n= {}$'.format(n) )
  plot.gcf().tight_layout()
  plot.savefig("plot_load_n{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  # plot_attack()
  # plot_load()
  optimal_pd()
  # plot_Pr_randmix_attack_sim_vs_model()
