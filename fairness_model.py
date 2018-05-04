import math, numpy, pprint

from patch import *
from rvs import *
from simplex_models import ESm_typei

r = 2
ET_UB = 17 # 20

def ff_har_ub_max(t, sdist_m):
  return 1/ESm_typei(1, t, 0, sdist_m)

def ff_ETh_lb(har, t, sdist_m):
  ES = ESm_typei(1, t, 0, sdist_m)
  ES2 = ESm_typei(2, t, 0, sdist_m)
  
  ET = ES + (har/2)*ES2/(1 - har*ES)
  return ET if ET < ET_UB else None

def ff_tmin(t):
  ns = t*r + 1
  n_sym = int(math.log(ns + 1, 2) )
  return int((ns - 1 - (n_sym-1)*2)/r)

def ff_har_ub_min(t, sdist_m):
  return 1/ESm_typei(1, ff_tmin(t), 0, sdist_m)

def ff_ETh_ub(har, t, sdist_m):
  tmin = ff_tmin(t)
  ES = ESm_typei(1, tmin, 0, sdist_m)
  ES2 = ESm_typei(2, tmin, 0, sdist_m)
  
  ET = ES + (har/2)*ES2/(1 - har*ES)
  return ET if ET >= 0 and ET < ET_UB else None

def ff_har_ub(car, t, sdist_m):
  # if t == 1:
  #   p_busyb = car/mu
  #   ES = p_busyb*1/mu + (1-p_busyb)*ES_avq(r, 1, mu)
  # elif t == 3:
  #   p_busyb = p_busyc = car/mu
  #   ES = p_busyb*p_busyc*ES_avq(r, 1, mu) + \
  #         ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*ES_avq(r, 2, mu) + \
  #         (1-p_busyb)*(1-p_busyc)*ES_avq(r, 3, mu)
  # return float(1/ES) - 0.01
  X = Exp(car)
  V = rv_from_m(sdist_m)
  V22 = X_n_k(V, 2, 2)
  Pr_V_l_X = mpmath.quad(lambda x: V.cdf(x)*X.pdf(x), [0, mpmath.inf] )
  m = int(numpy.log2(t+1) )
  def Pr_Sgs(s):
    a = V22.tail(s)*Pr_V_l_X + 1 - Pr_V_l_X
    ro = car*moment_ith(1, V)
    return V.tail(s)*V22.tail(s)**(t-m) * (a*(1 - ro) + ro)**m
  
  ES = float(mpmath.quad(lambda s: Pr_Sgs(s), [0, mpmath.inf] ) ) # 10000*10
  return 1/ES

def Pr_serv_completes_before_carr(car, t, mu):
  if t == 1:
    X = Exp(mu) # server service time
    def S_cdf(s): # request service time
      return 1 - X.tail(s)*(1 - X.cdf(s)**2)**t
    
    Y = Exp(car)
    return mpmath.quad(lambda y: S_cdf(y)*Y.pdf(y), [0, Y.u_l] )

def ff_ETh_approx_(har, car, t, sdist_m):
  X = Exp(car)
  V = rv_from_m(sdist_m)
  ro = car*moment_ith(1, V)
  m = int(numpy.log2(t+1) )
  
  V22 = X_n_k(V, 2, 2)
  # Pr_V22_l_X = mpmath.quad(lambda x: V22.cdf(x)*X.pdf(x), [X.l_l, X.u_l] )
  Pr_V_l_X = mpmath.quad(lambda x: V.cdf(x)*X.pdf(x), [0, mpmath.inf] )
  def Pr_Sgs(s):
    # a = V22.tail(s)*Pr_V22_l_X + 1 - Pr_V22_l_X
    a = V22.tail(s)*Pr_V_l_X + 1 - Pr_V_l_X
    return V.tail(s)*V22.tail(s)**(t-m) * (a*(1 - ro) + ro)**m
  
  ES = mpmath.quad(lambda s: Pr_Sgs(s), [0, mpmath.inf] ) # 10000*10
  ES2 = mpmath.quad(lambda s: 2*s*Pr_Sgs(s), [0, mpmath.inf] )
  # log(WARNING, "ES= {}, ES2= {}".format(ES, ES2) )
  ET = ES + har*ES2/2/(1 - har*ES) if har*ES < 1 else 10**3
  return ET if ET >= 0 and ET < ET_UB else None

def ff_ETh_newapprox(har, car, t, sdist_m):
  X = Exp(car)
  V = rv_from_m(sdist_m)
  ro = car*moment_ith(1, V)
  m = int(numpy.log2(t+1) )
  
  V22 = X_n_k(V, 2, 2)
  Pr_Vprime_g_v = lambda v: 1 - V.cdf(v)*scipy.integrate.quad(lambda x: math.exp(-car*x)*V.pdf(x), 0, v)[0]
  def Pr_S_g_s(s):
    return V.tail(s)*V22.tail(s)**(t-m) * (Pr_Vprime_g_v(s)*(1 - ro) + ro)**m
  
  # ES = scipy.integrate.quad(lambda s: Pr_S_g_s(s), 0, numpy.inf)[0]
  # ES2 = scipy.integrate.quad(lambda s: 2*s*Pr_S_g_s(s), 0, numpy.inf)[0]
  ES = mpmath.quad(lambda s: Pr_S_g_s(s), [0, mpmath.inf] ) # 10000*10
  ES2 = mpmath.quad(lambda s: 2*s*Pr_S_g_s(s), [0, mpmath.inf] )
  log(WARNING, "ES= {}, ES2= {}".format(ES, ES2) )
  ET = ES + har*ES2/2/(1 - har*ES) if har*ES < 1 else 10**3
  return ET if ET >= 0 and ET < ET_UB else None
  # ET = ES + har*ES2/2/(1 - har*ES)
  # return ET

def ff_ETh_approx(har, car, t, sdist_m):
  if sdist_m['dist'] == "Exp":
    V = Exp(sdist_m['mu'] )
  
  def Si_cdf(i, v): # Si ~ min{S_0, max{S_11, S_12}, ..., max{S_i1, S_i2}}
    return 1 - V.tail(v)*(1 - V.cdf(v)**2)**i
  
  def ESk_i(k, i): # kth moment of Si
    return mpmath.quad(lambda s: k*s**(k-1) * (1 - Si_cdf(i, s)), [0, mpmath.inf] )
  
  m = int(numpy.log2(t+1) )
  def EWk_i(k, i): # kth moment of Wi
    if i == 0:
      return ESk_i(k, t-m)
    
    Y = Exp(i*car)
    p = mpmath.quad(lambda v: Si_cdf(t-m+i, v)*Y.pdf(v), [0, mpmath.inf] ) # Y.u_l, Pr{S_{t-m+i} < Exp(i*car) }
    return p*ESk_i(k, t-m+i) + (1-p)*EWk_i(k, i-1)
  
  p_b = car*V.mean()
  def ESk(k):
    return sum([EWk_i(k, i)*binomial(m, i)*(1-p_b)**i * p_b**(m-i) for i in range(m+1)] )
  
  # if t == 1:
  #   mu = sdist_m['mu']
  #   p_busyb = car/mu
  #   p_fast = (1-p_busyb)*Pr_serv_completes_before_carr(car, t, mu)
  #   ES = (1-p_fast)*1/mu + p_fast*ES_avq(r, 1, mu)
  #   ES2 = (1-p_fast)*2/mu**2 + p_fast*ES2_avq(r, 1, mu)
  #   # print("ES= {}, ES2= {}".format(ES, ES2) )
  
  ES, ES2 = ESk(1), ESk(2)
  # log(WARNING, "ES= {}, ES2= {}".format(ES, ES2) )
  ET = ES + har*ES2/2/(1 - har*ES)
  return ET if ET >= 0 and ET < ET_UB else None

if __name__ == "__main__":
  t = 1
  car = 0.5
  sdist_m = {'dist': 'Exp', 'mu': 1}
  har_ub = ff_har_ub(car, t, sdist_m)
  print("har_ub= {}".format(har_ub) )
  for har in numpy.linspace(0.05, har_ub, 10):
    ETh = ff_ETh_approx(har, car, t, sdist_m)
    ETh_ = ff_ETh_approx_(har, car, t, sdist_m)
    print("har= {}, ETh= {}, ETh_= {}".format(har, ETh, ETh_) )
