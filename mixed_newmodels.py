import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot

import scipy
import numpy as np
import matplotlib.pyplot as plot
from scipy.optimize import fsolve

from patch import *
from rvs import *
from commonly_used import *

def Pr_B_g_t(pe, n, k, t, dist_m):
  return 1 - sum([binom(k-1, e)*pe**e*(1-pe)**(k-1-e) * Pr_Xnk_leq_x(n-k+1+e, e, t, dist_m) for e in range(k) ] )

def Pr_B_g_t_(pe, n, k, t, dist_m):
  cdf = 0
  if pe >= 1:
    Pr_E_l_k = 1
  else:
    Pr_E_l_k = sum([binom(n-1, e)*pe**e*(1-pe)**(n-1-e) for e in range(k) ] )
  if Pr_E_l_k == 0: Pr_E_l_k = 0.0001
  # print("pe= {}, Pr_E_l_k= {}".format(pe, Pr_E_l_k) )
  for e in range(k):
    cdf += binom(n-1, e)*pe**e*(1-pe)**(n-1-e) / Pr_E_l_k * Pr_Xnk_leq_x(n-k+1+e, e, t, dist_m)
  return 1 - cdf

def Pr_Be_g_t(pe, n, k, t, dist_m): # exceptional first service time
  cdf = 0
  for e in range(k):
    s = sum([Pr_Xnk_leq_x(n-k+1+j, j, t, dist_m) for j in range(e+1) ] )
    cdf += binom(k-1, e)*pe**e*(1-pe)**(k-1-e) * s/(e+1)
  return 1 - cdf
  
def Pr_Be_g_t_(pe, n, k, t, dist_m):
  cdf = 0
  if pe >= 1:
    Pr_E_l_k = 1
  else:
    Pr_E_l_k = sum([binom(n-1, e)*pe**e*(1-pe)**(n-1-e) for e in range(k) ] )
  for e in range(k):
    s = sum([Pr_Xnk_leq_x(n-k+1+j, j, t, dist_m) for j in range(e+1) ] )
    cdf += binom(n-1, e)*pe**e*(1-pe)**(n-1-e) / Pr_E_l_k * s/(e+1)
  return 1 - cdf

def EB_mth(pe, n, k, m, dist_m):
  return scipy.integrate.quad(lambda t: m*t**(m-1)*Pr_B_g_t(pe, n, k, t, dist_m), 0, np.inf)[0]
def EB_mth_(pe, n, k, m, dist_m):
  return scipy.integrate.quad(lambda t: m*t**(m-1)*Pr_B_g_t_(pe, n, k, t, dist_m), 0, np.inf)[0]
def EBe_mth(pe, n, k, m, dist_m):
  return scipy.integrate.quad(lambda t: m*t**(m-1)*Pr_Be_g_t(pe, n, k, t, dist_m), 0, np.inf)[0]
def EBe_mth_(pe, n, k, m, dist_m):
  return scipy.integrate.quad(lambda t: m*t**(m-1)*Pr_Be_g_t_(pe, n, k, t, dist_m), 0, np.inf)[0]

def pempty(n, k, dist_m):
  pe = 1
  mu = dist_m['mu']
  x_pdf = lambda x: mu*math.exp(-mu*x)
  for k_ in range(1, k+1):
    pe = scipy.integrate.quad(lambda t: (1 - Pr_B_g_t(pe, n, k_, t, dist_m) ) * x_pdf(t), 0, np.inf)[0]
  return pe

def pempty_mg1efsapprox(n, k, dist_m):
  ar = dist_m['mu']
  # pe = 1
  # x_pdf = lambda x: mu*math.exp(-mu*x)
  # for k_ in range(1, k+1):
  #   pe = scipy.integrate.quad(lambda t: (1 - Pr_Be_g_t(pe, n, k_, t, dist_m) ) * x_pdf(t), 0, np.inf)[0]
  # return pe
  
  def ro(pe):
    EB = EB_mth(pe, n, k, 1, dist_m)
    EBe = EBe_mth(pe, n, k, 1, dist_m)
    return ar*EBe/(1 - ar*(EB - EBe) )
  eq = lambda pe: pe - (1 - ro(pe) )
  # pe = scipy.optimize.brentq(eq, 0.01, 0.99)
  pe = scipy.optimize.brentq(eq, 0.00001, 1)
  return pe

def pempty_mg1efsapprox_(n, k, dist_m):
  ar = dist_m['mu']
  
  def ro(pe):
    EB = EB_mth_(pe, n, k, 1, dist_m)
    EBe = EBe_mth_(pe, n, k, 1, dist_m)
    return ar*EBe/(1 - ar*(EB - EBe) )
  eq = lambda pe: pe - (1 - ro(pe) )
  # pe = scipy.optimize.brentq(eq, 0.00001, 1)
  # pe = scipy.optimize.bisect(eq, 0, 1, xtol=1e-3)
  pe = scipy.optimize.fsolve(eq, 0)
  return pe

def ET_mg1efsapprox(n, k, dist_m):
  ar = dist_m['mu']
  pe = pempty_mg1efsapprox(n, k, dist_m)
  EB, EB2 = EB_mth(pe, n, k, 1, dist_m), EB_mth(pe, n, k, 2, dist_m)
  EBe, EBe2 = EBe_mth(pe, n, k, 1, dist_m), EBe_mth(pe, n, k, 2, dist_m)
  
  ET = EBe/(1 - ar*(EB - EBe) ) + ar*EB2/2/(1 - ar*EB) + ar*(EBe2 - EB2)/2/(1 - ar*(EB - EBe) )
  return ET

def ET_mg1efsapprox_(n, k, dist_m):
  ar = dist_m['mu']
  # pe = pempty(n, k, dist_m)
  pe = pempty_mg1efsapprox_(n, k, dist_m)
  EB, EB2 = EB_mth_(pe, n, k, 1, dist_m), EB_mth(pe, n, k, 2, dist_m)
  EBe, EBe2 = EBe_mth_(pe, n, k, 1, dist_m), EBe_mth(pe, n, k, 2, dist_m)
  
  ET = EBe/(1 - ar*(EB - EBe) ) + ar*EB2/2/(1 - ar*EB) + ar*(EBe2 - EB2)/2/(1 - ar*(EB - EBe) )
  return ET

def pempty_mg1approx(n, k, dist_m):
  ar = dist_m['mu']
  
  ro = lambda pe: ar*EB_mth(pe, n, k, 1, dist_m)
  eq = lambda pe: pe - (1 - ro(pe) )
  pe = scipy.optimize.brentq(eq, 0.00001, 1)
  return pe

def ET_mg1approx_(n, k, dist_m):
  pe = pempty_mg1approx(n, k, dist_m)
  EB, EB2 = EB_mth(pe, n, k, 1, dist_m), EB_mth(pe, n, k, 2, dist_m)
  # EBe, EBe2 = EBe_mth(pe, n, k, 1, dist_m), EBe_mth(pe, n, k, 2, dist_m)
  
  ar = dist_m['mu']
  # print("k= {}, pe= {}".format(k, pe) )
  
  ET = EB + ar*EB2/2/(1-ar*EB)
  # ET = EBe/(1 - ar*(EB - EBe) ) + ar*EB2/2/(1 - ar*EB) + ar*(EBe2 - EB2)/2/(1 - ar*(EB - EBe) )
  return ET

def ET_lb(n, k, dist_m):
  ET = 0
  for b in range(1, k+1):
    tail = lambda t: 1 - Pr_Xnk_leq_x(n-b, k-b, t, dist_m)
    ET += scipy.integrate.quad(lambda t: tail(t), 0, np.inf)[0]
  return ET/(k+1)

def ET_tighterlb(n, k, dist_m):
  ar = dist_m['mu']
  ET = 0
  for b in range(1, k+1):
    tail = lambda t: 1 - Pr_Xnk_leq_x(n-b, k-b, t, dist_m)
    EB = scipy.integrate.quad(lambda t: tail(t), 0, np.inf)[0]
    
    ET += EB * (1 + ar*EB)
  return ET/(k+1)

def ET_ub(n, k, dist_m):
  ar = dist_m['mu']
  """
  EB = lambda pd: 1/(n*ar*pd)
  EBe = lambda pd: (1 - pd)/(n*ar*pd)
  ro = lambda pd: ar*EB(pd)/(1 + ar*(EB(pd) - EBe(pd) ) )
  
  def eq(pd):
    ro_ = ro(pd)
    return pd - sum([binom(n, k_)*ro_**k_*(1-ro_)**(n-k_) for k_ in range(k, n+1) ] )
    # return pd - binom(n, k)*ro_**k*(1-ro_)**(n-k)
  pd = scipy.optimize.brentq(eq, 0.00001, 1)
  # pd = scipy.optimize.brentq(eq, -0.1, 1.1)
  print("pd= {}".format(pd) )
  
  EB, EB2 = 1/(n*ar*pd), 2/(n*ar*pd)**2
  print("EB= {}, EB2= {}".format(EB, EB2) )
  q = 1 - pd
  EBe, EBe2 = q/(n*ar*pd), q*2/(n*ar*pd)**2
  print("EBe= {}, EBe2= {}".format(EBe, EBe2) )
  """
  EB = lambda pd: 1/(ar*pd)
  EBe = lambda pd: (1 - pd)/(n*ar*pd)
  ro = lambda pd: ar*EB(pd)/(1 + ar*(EB(pd) - EBe(pd) ) )
  
  ET = EBe/(1 - ar*(EB - EBe) ) + ar*EB2/2/(1 - ar*EB) + ar*(EBe2 - EB2)/2/(1 - ar*(EB - EBe) )
  return ET

if __name__ == "__main__":
  n = 20
  dist_m = {'dist': 'Exp', 'mu': 1}
  print("n= {}, X ~ {}".format(n, dist_m) )
  # for k in range(2, n):
  for k in range(2, n):
    # pe = pempty(n, k, dist_m)
    # pe_approx = pempty_mg1efsapprox(n, k, dist_m)
    # print("k= {}, pe= {}, pe_approx= {}".format(k, pe, pe_approx) )
    
    ET_mg1efs = ET_mg1efsapprox(n, k, dist_m)
    ET_mg1efs_ = ET_mg1efsapprox_(n, k, dist_m)
    ET_mg1 = ET_mg1approx_(n, k, dist_m)
    # ETlb = ET_lb(n, k, dist_m)
    # ETtighterlb = ET_tighterlb(n, k, dist_m)
    # ETub = ET_ub(n, k, dist_m)
    
    print("k= {}".format(k) )
    print("ET_mg1efs= {}, ET_mg1efs_= {}, ET_mg1= {}".format(ET_mg1efs, ET_mg1efs_, ET_mg1) )
    # print("ETlb= {}, ETtighterlb= {}".format(ETlb, ETtighterlb) )
    # print("ETub= {}".format(ETub) )
