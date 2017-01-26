import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy, sympy, simpy, getopt
from math import factorial
from numpy import linalg
from patch import *

# -------------------------------------  MDS: Hot-Data Download  ---------------------------------- #
def E_T_sys__mds_n_2(arr_rate, mu, n):
  gamma = mu
  
  E_S_c = n/(gamma+(n-1)*mu) - (n-1)/(gamma+n*mu)
  E_S_c_2 = 2*n/(gamma+(n-1)*mu)**2 - 2*(n-1)/(gamma+n*mu)**2
  E_S_p = 1/(gamma+mu)
  E_S_p_2 = 2/(gamma+mu)**2
  
  tau = mu/(gamma + (n-1)*mu) # p_i_tor/p_i_tol
  # p_01/p_10 = n*tau
  p_0 = 1/(1+n*tau/(1-tau) )
  def p_i(i):
    if i == 1:
      return n*tau*p_0
    else:
      return n*tau**i * p_0
  nu = gamma + n*mu
  f_jd = 1 + (n-1)*mu/nu/(1 + (1-tau)/n/tau)
  f_to_c = gamma/nu * p_0 + (gamma+(n-1)*mu)/nu * p_i(1)
  f_c = f_to_c/f_jd
  
  E_S = f_c*E_S_c + (1-f_c)*E_S_p
  E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
  return E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)

# ------------------------------------  MDS: Full-Data Download  ---------------------------------- #
def E_T_mds_n_r_2(arr_rate, mu, n, r):
  return E_T_mds_n_2(r/n * arr_rate, mu, r)
  
def E_T_mds_n_2(arr_rate, mu, n):
  # high-regime assumption
  f_jc = (n - 2)/(n - 1)
  f_jp = 1 - f_jc
  mds_E_S_p = 1/((n - 1)*mu)
  mds_E_S_c = (1/mu)*(harmonic_sum(n) - harmonic_sum(n-2) )
  mds_E_S = f_jp*mds_E_S_p + f_jc*mds_E_S_c
  mds_E_S_p_2 = 2/(((n - 1)**2)*(mu**2) )
  mds_E_S_c_2 = (1/(mu**2) )*(harmonic_2_sum(n) - harmonic_2_sum(n-2) ) + mds_E_S_c**2
  mds_E_S_2 = f_jp*mds_E_S_p_2 + f_jc*mds_E_S_c_2
  # mds_E_T = mds_E_S/(1 - arr_rate*mds_E_S) # w/ M/M/1 assumption
  mds_E_T = mds_E_S + (arr_rate/2)*mds_E_S_2/(1 - arr_rate*mds_E_S) # w/ M/G/1 assumption
  return mds_E_T
  
def E_T_mds_n_2_adjustable(arr_rate, mu, n):
  # gamma = min(arr_rate, mu)
  # gamma = min(arr_rate*mu/(arr_rate+mu), \
  #             mu/(harmonic_sum(n) - harmonic_sum(n-2) ) )
  gamma = mu
  ro = gamma/mu
  p_0 = 1/(1 + n*ro/(n - 1 - ro) )
  def p_i(i):
    return (ro/(n - 1) )**i * n*p_0
  
  nu = (n-1)*mu + gamma
  sum_for_pi = p_0*n*gamma + nu*(1-p_0)
  pi_0 = p_0*n*gamma/sum_for_pi
  pi_1 = p_i(1)*nu/sum_for_pi
  
  f_jd = (n-1)*mu/nu*(1-pi_0)
  f_c = pi_1*(n-1)*mu/nu
  f_jc = f_c/f_jd
  f_jp = 1 - f_jc
  
  mds_E_S_p = 1/((n - 1)*mu)
  mds_E_S_c = (1/mu)*(harmonic_sum(n) - harmonic_sum(n-2) )
  mds_E_S = f_jp*mds_E_S_p + f_jc*mds_E_S_c
  mds_E_S_p_2 = 2/(((n - 1)**2)*(mu**2) )
  mds_E_S_c_2 = (1/(mu**2) )*(harmonic_2_sum(n) - harmonic_2_sum(n-2) ) + mds_E_S_c**2
  mds_E_S_2 = f_jp*mds_E_S_p_2 + f_jc*mds_E_S_c_2
  mds_E_T = mds_E_S + (arr_rate/2)*mds_E_S_2/(1 - arr_rate*mds_E_S)
  return mds_E_T

def E_T_mds_n_k_sm(arr_rate, mu, n, k): # works for k >= 1
  harmonic_diff = harmonic_sum(n) - harmonic_sum(n-k)
  E_S = harmonic_diff/mu
  if 1/arr_rate < E_S:
    return None
  E_S_2 = (harmonic_2_sum(n)-harmonic_2_sum(n-k) ) / mu**2 + E_S**2
  return E_S + arr_rate*E_S_2/2/(1 - arr_rate*E_S)

def adj_E_T_mds_n_k_sm(arr_rate, mu, n, k):
  return E_T_mds_n_k_sm(arr_rate, mu, n, k-1) + E_T_mds_n_k_sm(arr_rate, mu, n-k+1, 1)

def adj_2_E_T_mds_n_k_sm(arr_rate, mu, n, k):
  return E_T_mds_n_k_sm(arr_rate, mu, n, k-2) + E_T_mds_n_k_sm(arr_rate, mu, n-k+2, 1) + E_T_mds_n_k_sm(arr_rate, mu, n-k+1, 1)

def E_T_mds_n_k_sm_recur(arr_rate, mu, n, k):
  if n > k:
    if k > 1:
      return E_T_mds_n_k_sm_recur(arr_rate, mu, n, k-1) + E_T_mds_n_k_sm(arr_rate, mu, n-k+1, 1)
    elif k == 1:
      return E_T_mds_n_k_sm(arr_rate, mu, n, k)
    else:
      log(ERROR, "Unexpected k= {}".format(k) )
      sys.exit(1)
  else:
    log(ERROR, "Unexpected n= {} <= k= {}".format(n, k) )
    return E_T_mds_n_k_sm(arr_rate, mu, n, k)

def mds_inner_bound_on_arr_rate(mu, n, k):
  sm_mds_n_k_E_S = 1/mu * (harmonic_sum(n) - harmonic_sum(n-k) )
  return 1/sm_mds_n_k_E_S

def mds_exact_bound_on_arr_rate(mu, n, k):
  return n*mu/k

def E_T_mds_n_k_varki_gauri_lb(arr_rate, mu, n, k):
  ro = float(arr_rate/mu)
  return 1/mu * (harmonic_sum(n) - harmonic_sum(n-k) ) + \
         1/mu * ro*(gen_harmonic_sum(n, ro) - gen_harmonic_sum(n-k, ro) )