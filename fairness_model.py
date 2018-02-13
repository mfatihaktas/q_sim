import math, numpy, pprint

from patch import *
from rvs import *
from simplex_models import ES_avq, ES2_avq

r = 2
ET_UB = 50

def ff_har_ub_max(t, sdist_m):
  return 1/ES_typei_(1, t, t, sdist_m)

def ff_ETh_lb(har, t, sdist_m):
  ES = ES_typei_(1, t, t, sdist_m)
  ES2 = ES_typei_(2, t, t, sdist_m)
  
  ET = ES + (har/2)*ES2/(1 - har*ES)
  if ET > ET_UB: return None
  return ET

def ff_har_ub_min(t, mu):
  n_server = t*r + 1
  n_sym = int(math.log(n_server+1, 2) )
  t_min = int((n_server - 1 - (n_sym-1)*2)/r)
  ES = ES_avq(r, t_min, mu)
  return float(1/ES) - 0.01

def ff_ETh_ub(har, t, mu):
  n_server = t*r + 1
  n_sym = int(math.log(n_server+1, 2) )
  t_min = int((n_server - 1 - (n_sym-1)*2)/r)
  
  ES = ES_avq(r, t_min, mu)
  ES2 = ES2_avq(r, t_min, mu)
  
  ET = ES + (har/2)*ES2/(1 - har*ES)
  if ET < 0 or ET > ET_UB: return None
  return ET

def ff_har_ub(car, t, mu):
  # if t == 1:
  #   p_busyb = car/mu
  #   ES = p_busyb*1/mu + (1-p_busyb)*ES_avq(r, 1, mu)
  # elif t == 3:
  #   p_busyb = p_busyc = car/mu
  #   ES = p_busyb*p_busyc*ES_avq(r, 1, mu) + \
  #         ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*ES_avq(r, 2, mu) + \
  #         (1-p_busyb)*(1-p_busyc)*ES_avq(r, 3, mu)
  # return float(1/ES) - 0.01
  X = rv_from_m(dist_m)
  

def Pr_serv_completes_before_carr(car, t, mu):
  if t == 1:
    X = Exp(mu) # server service time
    def S_cdf(s): # request service time
      return 1 - X.tail(s)*(1 - X.cdf(s)**2)**t
    
    Y = Exp(car)
    return mpmath.quad(lambda y: S_cdf(y)*Y.pdf(y), [0, Y.u_l] )

"""
def ff_ETh_approx(har, car, t, serv, serv_dist_m):
  mu = serv_dist_m['mu']
  if t == 1:
    p_busyb = car/mu
    p_fast = (1-p_busyb)*Pr_serv_completes_before_carr(car, t, mu)
    
    ES = (1-p_fast)*1/mu + p_fast*ES_avq(r, 1, mu)
    ES2 = (1-p_fast)*2/mu**2 + p_fast*ES2_avq(r, 1, mu)
  elif t == 3:
    p_busyb = p_busyc = car/mu
    
    ES = p_busyb*p_busyc*ES_avq(r, 1, mu) + \
          ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*ES_avq(r, 2, mu) + \
          (1-p_busyb)*(1-p_busyc)*ES_avq(r, 3, mu)
    ES2 = p_busyb*p_busyc*ES2_avq(r, 1, mu) + \
            ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*ES2_avq(r, 2, mu) + \
            (1-p_busyb)*(1-p_busyc)*ES2_avq(r, 3, mu)
  else:
    log(ERROR, "Unexpected t= {}".format(t) )
  
  ET = ES + (har/2)*ES2/(1 - har*ES)
  if ET < 0 or ET > ET_UB: return None
  return ET
"""
def ff_ETh_approx(har, car, t, serv, serv_dist_m):
  if serv == "Exp":
    S = Exp(serv_dist_m['mu'] )
  
  def V_cdf(i, v): # V_i ~ min{S_0, max{S_11, S_12}, ..., max{S_i1, S_i2}}
    return 1 - (1 - S.cdf(v))*(1 - S.cdf(v)**2)**i
  
  def EV_i_j(i, j): # jth moment of V_i
    return mpmath.quad(lambda v: j*v**(j-1) * (1 - V_cdf(i, v)), [0, S.u_l] )
  
  m = int(numpy.log2(t+1) )
  def E_W_i_j(i, j): # jth moment of W_i
    if i == 0:
      return EV_i_j(t-m, j)
    # Pr{V_i < Exp((m-i)mu) }
    def p():
      Y = Exp(i*car)
      return mpmath.quad(lambda v: V_cdf(t-m+i, v)*Y.pdf(v), [0, Y.u_l] )
    p = p()
    return p*EV_i_j(t-m+i, j) + (1-p)*E_W_i_j(i-1, j)
  
  p_b = car*S.mean()
  def E_V_j(j): # jth moment of V
    return sum([E_W_i_j(i, j)*binomial(m, i)*(1-p_b)**i * p_b**(m-i) for i in range(m+1)] )
  
  if t == 1:
    mu = serv_dist_m['mu']
    p_busyb = car/mu
    p_fast = (1-p_busyb)*Pr_serv_completes_before_carr(car, t, mu)
    # print("prev:")
    # print("E_V_1= {}, E_V_0= {}".format(ES_avq(r, 1, mu), 1/mu) )
    E_V = (1-p_fast)*1/mu + p_fast*ES_avq(r, 1, mu)
    E_V_2 = (1-p_fast)*2/mu**2 + p_fast*ES2_avq(r, 1, mu)
    # print("E_V= {}, E_V_2= {}".format(E_V, E_V_2) )
  
  # print("now:")
  # print("E_V_1= {}, E_V_0= {}".format(EV_i_j(1, 1), EV_i_j(0, 1)) )
  # print("E_W_1= {}, E_W_0= {}".format(E_W_i_j(1, 1), E_W_i_j(0, 1) ) )
  E_V, E_V_2 = E_V_j(1), E_V_j(2)
  # print("E_V= {}, E_V_2= {}".format(E_V, E_V_2) )
  ET = E_V + har*E_V_2/2/(1 - har*E_V)
  if ET < 0 or ET > ET_UB: return None
  return ET

