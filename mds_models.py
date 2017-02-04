import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy, sympy, simpy, getopt
from math import factorial
from numpy import linalg
from patch import *

# -------------------------------------  MDS: Hot-Data Download  ---------------------------------- #
def E_T_sys__mds_n_2(l, mu, n):
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
  return E_S + (l/2)*E_S_2/(1 - l*E_S)

# ------------------------------------  MDS: Full-Data Download  ---------------------------------- #
def E_T_mds_n_r_2(l, mu, n, r):
  # Exact when n is divided into (r,2)'s where each arriving job is sent to one uniformly, independently
  return E_T_mds_n_2(r/n * l, mu, r)

def E_T_fj_n_2(l, mu, n, l_s):
  return 1/(l_s + l*2/n) * H(2)

def E_T_mds_n_2_no_task_cancel_lb(l, mu, n, l_s):
  # Arriving jobs are split into n but once the job is completed, remaining tasks don't get cancelled
  return 1/(l_s + l) * (H(n) - H(n-2) )

def E_T_mds_n_k(l, mu, n, k, l_s=None):
  if n == 1:
    return 0
  elif k == 1:
    return 1/(n*mu - l)
  elif k == 2:
    return E_T_mds_n_2(l, mu, n, l_s)
  
def E_T_mds_n_2(l, mu, n, l_s=None):
  # """
  if n == 2:
    # if l_s is not None:
    #   mu = mu - l_s
    # ro = l/mu
    # return (12 - ro)/(8*(mu - l) )
    
    # A LB
    # ro = l/mu
    # return (12 - ro)/(8*(mu - l) ) - 1/(mu - l) + 1/(mu - l_s - l)
    # An approximation
    if l_s is not None:
      l = l + l_s
    ro = l/mu
    return (12 - ro)/(8*(mu - l) )
    
    # if l_s is not None:
    #   ro = l/mu
    #   # return (12 - ro)/(8*(mu - l) ) + 1/(mu - l_s - l)
    #   return (12 - ro)/(8*(mu - l) ) * (l+l_s)/l
    '''
    # W/ High-traffic approximation
    d = 0.02 # mu/4
    tau = mu/(mu+d)
    ro = 2*tau + d/(mu+d)
    p_0 = 1/(1 + ro/(1-tau) )
    
    f_jd = ro*(1 - p_0)
    f_to_c = ro * ro*p_0
    f_c = f_to_c/f_jd
    print("f_c= {}".format(f_c) )
    
    E_S_c = 1/mu * H(2)
    E_S_c_2 = 5 / mu**2
    E_S_p = 1/mu
    E_S_p_2 = 1 / mu**2
    
    E_S = f_c*E_S_c + (1-f_c)*E_S_p
    E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
    
    return E_S + l/2 * E_S_2/(1 - l*E_S)
    '''
  # """
  if l_s is None:
    # high-regime assumption
    f_jc = (n - 2)/(n - 1)
    f_jp = 1 - f_jc
    E_S_p = 1/((n - 1)*mu)
    E_S_c = 1/mu * (H(n) - H(n-2) )
    if n == 2:
      E_S_p = 1/mu
      E_S_c = 1/mu * H(2)
    E_S = f_jp*E_S_p + f_jc*E_S_c
    E_S_p_2 = 2 / ((n-1)*mu)**2
    E_S_c_2 = 1/(mu**2) * (H_2(n) - H_2(n-2) ) + E_S_c**2
    E_S_2 = f_jp*E_S_p_2 + f_jc*E_S_c_2
    # mds_E_T = E_S/(1 - l*E_S) # w/ M/M/1 assumption
    mds_E_T = E_S + l/2 * E_S_2/(1 - l*E_S) # w/ M/G/1 assumption
    return mds_E_T
  else:
    # A lower-bound -- TODO: turned out to be an upper-bound, don't know WHY
    # lambda_ = l_s + 2/n * l
    # return 1/(mu - lambda_) * (H(n) - H(n-2) )
    
    # An UB: no task-cancellation after job termination
    return 1/(mu - l - l_s) * (H(n) - H(n-2) )
    
    # return E_T_mds_n_2(l, mu-l_s, n) # l/(l+l_s) * E_T_mds_n_2(l, mu, n) + l_s/(l+l_s) * (H(n)-H(n-2) )/(mu-l_s)
    # return E_T_mds_n_2(l, mu*l/(l+l_s), n)
    # return E_T_mds_n_2(l, mu, int(n*(1-l_s/mu) ) )

def E_T_mds_n_2_adjustable(l, mu, n):
  # gamma = min(l, mu)
  # gamma = min(l*mu/(l+mu), \
  #             mu/(H(n) - H(n-2) ) )
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
  
  E_S_p = 1/((n - 1)*mu)
  E_S_c = (1/mu)*(H(n) - H(n-2) )
  E_S = f_jp*E_S_p + f_jc*E_S_c
  E_S_p_2 = 2/(((n - 1)**2)*(mu**2) )
  E_S_c_2 = (1/(mu**2) )*(H_2(n) - H_2(n-2) ) + E_S_c**2
  E_S_2 = f_jp*E_S_p_2 + f_jc*E_S_c_2
  mds_E_T = E_S + (l/2)*E_S_2/(1 - l*E_S)
  return mds_E_T

def E_T_mds_n_k_sm(l, mu, n, k): # works for k >= 1
  E_S = (H(n) - H(n-k) )/mu
  if 1/l <= E_S:
    return None
  E_S_2 = (H_2(n)-H_2(n-k) ) / mu**2 + E_S**2
  return E_S + l*E_S_2/2/(1 - l*E_S)

def adj_E_T_mds_n_k_sm(l, mu, n, k):
  return E_T_mds_n_k_sm(l, mu, n, k-1) + E_T_mds_n_k_sm(l, mu, n-k+1, 1)

def adj_2_E_T_mds_n_k_sm(l, mu, n, k):
  return E_T_mds_n_k_sm(l, mu, n, k-2) + E_T_mds_n_k_sm(l, mu, n-k+2, 1) + E_T_mds_n_k_sm(l, mu, n-k+1, 1)

def E_T_mds_n_k_sm_recur(l, mu, n, k):
  if n > k:
    if k > 1:
      return E_T_mds_n_k_sm_recur(l, mu, n, k-1) + E_T_mds_n_k_sm(l, mu, n-k+1, 1)
    elif k == 1:
      return E_T_mds_n_k_sm(l, mu, n, k)
    else:
      log(ERROR, "Unexpected k= {}".format(k) )
      sys.exit(1)
  else:
    log(ERROR, "Unexpected n= {} <= k= {}".format(n, k) )
    return E_T_mds_n_k_sm(l, mu, n, k)

def mds_inner_bound_on_arr_rate(mu, n, k):
  if n == k and k == 2:
    return mu
  sm_mds_n_k_E_S = 1/mu * (H(n) - H(n-k) )
  return 1/sm_mds_n_k_E_S

def mds_exact_bound_on_arr_rate(mu, n, k):
  return n*mu/k

def E_T_mds_n_k_varki_gauri_lb(l, mu, n, k):
  ro = float(l/mu)
  return 1/mu * (H(n) - H(n-k) ) + \
         1/mu * ro*(gen_H(n, ro) - gen_H(n-k, ro) )
