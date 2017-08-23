import math, numpy, pprint

from patch import *
from rvs import *
from simplex_models import E_S_avq, E_S_2_avq

r = 2
E_T_UB = 50

def hot_ar_ub_max_ff_simplex(t, mu):
  E_S = E_S_avq(r, t, mu)
  return float(1/E_S) - 0.01

def E_T_hot_lb_ff_simplex(hot_ar, t, mu):
  E_S = E_S_avq(r, t, mu)
  E_S_2 = E_S_2_avq(r, t, mu)
  
  E_T = E_S + (hot_ar/2)*E_S_2/(1 - hot_ar*E_S)
  if E_T > E_T_UB: return None
  return E_T

def hot_ar_ub_min_ff_simplex(t, mu):
  n_server = t*r + 1
  n_sym = int(math.log(n_server+1, 2) )
  t_min = int((n_server - 1 - (n_sym-1)*2)/r)
  E_S = E_S_avq(r, t_min, mu)
  return float(1/E_S) - 0.01

def E_T_hot_ub_ff_simplex(hot_ar, t, mu):
  n_server = t*r + 1
  n_sym = int(math.log(n_server+1, 2) )
  t_min = int((n_server - 1 - (n_sym-1)*2)/r)
  
  E_S = E_S_avq(r, t_min, mu)
  E_S_2 = E_S_2_avq(r, t_min, mu)
  
  E_T = E_S + (hot_ar/2)*E_S_2/(1 - hot_ar*E_S)
  if E_T < 0 or E_T > E_T_UB: return None
  return E_T

def hot_ar_ub_ff_simplex_approx(cold_ar, t, mu):
  if t == 1:
    p_busyb = cold_ar/mu
    E_S = p_busyb*1/mu + (1-p_busyb)*E_S_avq(r, 1, mu)
  elif t == 3:
    p_busyb = p_busyc = cold_ar/mu
    E_S = p_busyb*p_busyc*E_S_avq(r, 1, mu) + \
          ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*E_S_avq(r, 2, mu) + \
          (1-p_busyb)*(1-p_busyc)*E_S_avq(r, 3, mu)
  return float(1/E_S) - 0.01

def Pr_serv_completes_before_cold_arr(cold_ar, t, mu):
  if t == 1:
    X = Exp(mu) # server service time
    def S_cdf(s): # request service time
      return 1 - X.tail(s)*(1 - X.cdf(s)**2)**t
    
    Y = Exp(cold_ar)
    return mpmath.quad(lambda y: S_cdf(y)*Y.pdf(y), [0, Y.u_l] )

"""
def E_T_hot_approx_ff_simplex(hot_ar, cold_ar, t, serv, serv_dist_m):
  mu = serv_dist_m['mu']
  if t == 1:
    p_busyb = cold_ar/mu
    p_fast = (1-p_busyb)*Pr_serv_completes_before_cold_arr(cold_ar, t, mu)
    
    E_S = (1-p_fast)*1/mu + p_fast*E_S_avq(r, 1, mu)
    E_S_2 = (1-p_fast)*2/mu**2 + p_fast*E_S_2_avq(r, 1, mu)
  elif t == 3:
    p_busyb = p_busyc = cold_ar/mu
    
    E_S = p_busyb*p_busyc*E_S_avq(r, 1, mu) + \
          ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*E_S_avq(r, 2, mu) + \
          (1-p_busyb)*(1-p_busyc)*E_S_avq(r, 3, mu)
    E_S_2 = p_busyb*p_busyc*E_S_2_avq(r, 1, mu) + \
            ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*E_S_2_avq(r, 2, mu) + \
            (1-p_busyb)*(1-p_busyc)*E_S_2_avq(r, 3, mu)
  else:
    log(ERROR, "Unexpected t= {}".format(t) )
  
  E_T = E_S + (hot_ar/2)*E_S_2/(1 - hot_ar*E_S)
  if E_T < 0 or E_T > E_T_UB: return None
  return E_T
"""
def E_T_hot_approx_ff_simplex(hot_ar, cold_ar, t, serv, serv_dist_m):
  if serv == "Exp":
    S = Exp(serv_dist_m['mu'] )
  
  def V_cdf(i, v): # V_i ~ min{S_0, max{S_11, S_12}, ..., max{S_i1, S_i2}}
    return 1 - (1 - S.cdf(v))*(1 - S.cdf(v)**2)**i
  
  def E_V_i_j(i, j): # jth moment of V_i
    return mpmath.quad(lambda v: j*v**(j-1) * (1 - V_cdf(i, v)), [0, S.u_l] )
  
  m = int(numpy.log2(t+1) )
  def E_W_i_j(i, j): # jth moment of W_i
    if i == 0:
      return E_V_i_j(t-m, j)
    # Pr{V_i < Exp((m-i)mu) }
    def p():
      Y = Exp(i*cold_ar)
      return mpmath.quad(lambda v: V_cdf(t-m+i, v)*Y.pdf(v), [0, Y.u_l] )
    p = p()
    return p*E_V_i_j(t-m+i, j) + (1-p)*E_W_i_j(i-1, j)
  
  p_b = cold_ar*S.mean()
  def E_V_j(j): # jth moment of V
    return sum([E_W_i_j(i, j)*binomial(m, i)*(1-p_b)**i * p_b**(m-i) for i in range(m+1)] )
  
  if t == 1:
    mu = serv_dist_m['mu']
    p_busyb = cold_ar/mu
    p_fast = (1-p_busyb)*Pr_serv_completes_before_cold_arr(cold_ar, t, mu)
    # print("prev:")
    # print("E_V_1= {}, E_V_0= {}".format(E_S_avq(r, 1, mu), 1/mu) )
    E_V = (1-p_fast)*1/mu + p_fast*E_S_avq(r, 1, mu)
    E_V_2 = (1-p_fast)*2/mu**2 + p_fast*E_S_2_avq(r, 1, mu)
    # print("E_V= {}, E_V_2= {}".format(E_V, E_V_2) )
  
  # print("now:")
  # print("E_V_1= {}, E_V_0= {}".format(E_V_i_j(1, 1), E_V_i_j(0, 1)) )
  # print("E_W_1= {}, E_W_0= {}".format(E_W_i_j(1, 1), E_W_i_j(0, 1) ) )
  E_V, E_V_2 = E_V_j(1), E_V_j(2)
  # print("E_V= {}, E_V_2= {}".format(E_V, E_V_2) )
  E_T = E_V + hot_ar*E_V_2/2/(1 - hot_ar*E_V)
  if E_T < 0 or E_T > E_T_UB: return None
  return E_T

