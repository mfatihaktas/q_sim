import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy, sympy, simpy, getopt
from math import factorial
from numpy import linalg

from patch import *
from commonly_used import *

'''
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
def mds_inner_bound_on_arr_rate(mu, n, k):
  if n == k and k == 2:
    return mu
  sm_mds_n_k_E_S = 1/mu * (H(n) - H(n-k) )
  return 1/sm_mds_n_k_E_S

def mds_exactbound_on_ar(mu, n, k):
  return n*mu/k

def E_T_mds_n_r_2(l, mu, n, r):
  # Exact when n is divided into (r,2)'s where each arriving job is sent to one uniformly, independently
  return E_T_mds_n_2(r/n * l, mu, r)

def E_T_fj_n_2(l, mu, n, l_s):
  return 1/(l_s + l*2/n) * H(2)

def E_T_mds_n_2_no_task_cancel_lb(l, mu, n, l_s):
  # Arriving jobs are split into n but once the job is completed, remaining tasks don't get cancelled
  return 1/(l_s + l) * (H(n) - H(n-2) )

def E_S_mds_type_i(mu, n, k, i): # i in [0, k-1]
  return 1/mu * (H(n-i) - H(n-k) )
def E_S_2_mds_type_i(mu, n, k, i): # i in [0, k-1]
  return 1/mu**2 * (H_2(n-i) - H_2(n-k) ) + E_S_mds_type_i(mu, n, k, i)**2

def plot_moments_for_types():
  mu = 1
  n, k = 20, 15
  gamma = mu
  i_l = []
  E_S_mds_type_i_l, E_S_2_mds_type_i_l = [], []
  E_S_rep_t_l, E_S_2_rep_t_l = [], []
  for i in range(1, k):
    i_l.append(i)
    E_S_mds_type_i_l.append(E_S_mds_type_i(mu, n, k, i) )
    E_S_2_mds_type_i_l.append(E_S_2_mds_type_i(mu, n, k, i) )
    # E_S_rep_t_l.append(1/(t+1)/mu)
    # E_S_2_rep_t_l.append(2/((t+1)*mu)**2 )
  plot.plot(i_l, E_S_mds_type_i_l, color=next(dark_color), label=r'$E[S], simplex$', marker=next(marker), linestyle=':', mew=2)
  plot.plot(i_l, E_S_2_mds_type_i_l, color=next(dark_color), label=r'$E[S^2]$, simplex', marker=next(marker), linestyle=':', mew=2)
  # plot.plot(i_l, E_S_rep_t_l, color=next(dark_color), label=r'$E[S]$, rep', marker=next(marker), linestyle=':', mew=2)
  # plot.plot(i_l, E_S_2_rep_t_l, color=next(dark_color), label=r'$E[S^2]$, rep', marker=next(marker), linestyle=':', mew=2)
  
  plot.legend()
  plot.xlabel(r'$i$')
  plot.ylabel("")
  plot.title(r'$mu= {}, n= {}, k= {}$'.format(mu, n, k) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_moments_for_types.pdf")
  plot.gcf().clear()
  log(WARNING, "done.")

def E_T_mds_n_k(l, mu, n, k, p_i_l=[], l_s=None, naive=False, incremental=False):
  if k == 1:
    return 1/(n*mu - l)
  elif k == 2:
    return E_T_mds_n_2(l, mu, n, l_s)
  else:
    E_S_mds_type_i_l, E_S_2_mds_type_i_l = [], []
    for i in range(k):
      E_S_mds_type_i_l.append(E_S_mds_type_i(mu, n, k, i) )
      E_S_2_mds_type_i_l.append(E_S_2_mds_type_i(mu, n, k, i) )
    # p_i_l = [1/k for i in range(k) ]
    
    if len(p_i_l) == 0:
      E_X = 1/l
      E_S_min = sum(E_S_mds_type_i_l)/len(E_S_mds_type_i_l)
      ro_max = E_S_min/E_X
      
      t = k-1
      p_0 = (1-ro_max)/(1-ro_max**(t+1) )
      def p_i(i):
        if naive:
          return 1/(t+1)
        else:
          return ro_max**i * p_0
      p_i_l = [p_i(i) for i in range(t+1) ]
      # Improves estimates of ro_m's incrementally
      if incremental:
        ro_i_l = []
        for i in range(t+1):
          # print("i= {}".format(i) )
          A = sum([numpy.prod(ro_i_l[:j] ) for j in range(i+1) ] )
          
          E_Y = E_X - E_S_min
          B = (t-i)*numpy.prod(ro_i_l[:i] )
          # print("A= {}, B= {}, (E_X-E_Y*A)/B/E_Y= {}".format(A, B, (E_X-E_Y*A)/B/E_Y) )
          ro_i = min((E_X-E_Y*A)/B/E_Y, 1)/2
          ro_i_l.append(ro_i)
          p_0 = 1/(sum([numpy.prod(ro_i_l[:j] ) for j in range(i+1) ] ) + numpy.prod(ro_i_l[:i+1] )*(t-i) )
          for i_ in range(t+1):
            p_i_l[i_] = numpy.prod(ro_i_l[:i_] )*p_0
          # print("p_0= {}, ro_i_l= {}, p_i_l= {}".format(p_0, ro_i_l, p_i_l) )
    E_S = sum([E_S_mds_type_i_l[i]*p_i for i,p_i in enumerate(p_i_l) ] )
    E_S_2 = sum([E_S_2_mds_type_i_l[i]*p_i for i,p_i in enumerate(p_i_l) ] )
    E_T = E_S + l*E_S_2/2/(1-l*E_S)
    if E_T < 0 or E_T > 30: return None
    return E_T
  
def E_T_mds_n_2(l, mu, n, l_s=None):
  if n == 2:
    # An approximation
    if l_s is not None:
      l = l + l_s
    ro = l/mu
    return (12 - ro)/(8*(mu - l) )
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
    E_S_c_2 = 1/mu**2 * (H_2(n) - H_2(n-2) ) + E_S_c**2
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

'''
def ET_mds_nk_varkigauri_lb(ar, n, k, dist_m):
  if dist_m['dist'] == 'Exp':
    mu = dist_m['mu']
    ro = float(ar/mu)
    return 1/mu * (H(n) - H(n-k) ) + \
           1/mu * ro*(gen_H(n, ro) - gen_H(n-k, ro) )

# ###################################  (n, 2)  ############################## #
def mds_innerbound_on_ar(n, k, dist_m):
  dist = dist_m['dist']
  if dist == 'Exp':
    mu = dist_m['mu']
    if n == k and k == 2:
      return mu
    EV_sm = 1/mu * (H(n) - H(n-k) )
    return 1/EV_sm

def mds_exactbound_on_ar(n, k, dist_m):
  dist = dist_m['dist']
  if dist == 'Exp':
    mu = dist_m['mu']
    return n*mu/k
  else:
    EV = EXm(1, dist_m)
    return float(n/k/EV)

def ET_mds_nk_sm(ar, n, k, dist_m):
  dist = dist_m['dist']
  if dist == 'Exp':
    mu = dist_m['mu']
    EV = (H(n) - H(n-k) )/mu
    if 1/ar <= EV:
      return None
    EV2 = (H_2(n)-H_2(n-k) ) / mu**2 + EV**2
    return EV + ar*EV2/2/(1 - ar*EV)
  else:
    EV = EXm_n_k(1, n, k, dist_m)
    if 1/ar <= EV:
      return None
    EV2 = EXm_n_k(2, n, k, dist_m)
    return EV + ar*EV2/2/(1 - ar*EV)

def V_tail(pd, n, k, t, dist_m):
  cdf = 0
  for d in range(k):
    cdf += binomial(k-1, d) * pd**d * (1-pd)**(k-1-d) * Pr_X_n_k_leq_x(n-d, k-d, t, dist_m)
  return 1 - cdf

def V_moment(pd, n, k, m, dist_m):
  return mpmath.quad(lambda t: m*t**(m-1)*V_tail(pd, n, k, t, dist_m), [0, mpmath.inf] ) # [0, 100000]

def V_moment(n, k, m, dist_m, pi_l):
  def Vi_moment(i):
    return mpmath.quad(lambda t: m*t**(m-1)*(1 - Pr_X_n_k_leq_x(n-i, k-i, t, dist_m) ), [0, mpmath.inf] )
  
  # return sum([1/k*Vi_moment(i) for i in range(0, k) ] )
  return sum([pi_l[i]*Vi_moment(i) for i in range(0, k) ] )

def pdeparted(n, k, dist_m):
  pd = 1 # 0.001 # 1
  
  # for i in range(5):
  #   pd = mpmath.quad(lambda t: V_tail(pd, n, k, t, dist_m) * f(t, dist_m), [0, mpmath.inf] )
  #   print("i= {}, pd= {}".format(i, pd) )
  
  for k_ in range(1, k+1):
    pd = mpmath.quad(lambda t: V_tail(pd, n, k_, t, dist_m) * f(t, dist_m), [0, mpmath.inf] )
    print("k_= {}, pd= {}".format(k_, pd) )
  return pd

def ET_mds_nk_approx(ar, n, k, dist_m, pi_l=[] ):
  dist = dist_m['dist']
  if dist == 'Exp':
    '''
    mu = dist_m['mu']
    # pd = pdeparted(n, k, dist_m)
    # EV = V_moment(pd, n, k, 1, dist_m)
    # EV2 = V_moment(pd, n, k, 2, dist_m)
    EV = V_moment(n, k, 1, dist_m, pi_l)
    EV2 = V_moment(n, k, 2, dist_m, pi_l)
    ET = EV + ar/2 * EV2/(1 - ar*EV)
    
    # EVf = EXm_n_k(1, n, k, dist_m)
    # EVf2 = EXm_n_k(2, n, k, dist_m)
    # rof = ar*EVf*(1 - ar*EV)/(1 + ar*EVf - ar*EV)
    # ro = ar*EVf*ar*EV/(1 + ar*EVf - ar*EV)
    # ER = EV2/2/EV
    # ERf = EVf2/2/EVf
    # ET = (rof+ro)*EV + (1-rof-ro)*EVf + (rof*ERf + ro*ER)/(1 - ar*EV)
    if ET < 0: return None
    return ET
    '''
    n_ = n+(k-1)
    # dist_m_ = dict(dist_m)
    # dist_m_['mu'] *= n/(n-k+1)
    return ET_mds_nk_sm(ar, n_, k, dist_m)

def ET_mds_n2_approx(ar, n, dist_m):
  dist = dist_m['dist']
  if dist == 'Exp':
    mu = dist_m['mu']
    if n == 2:
      # An approximation
      ro = ar/mu
      return (12 - ro)/(8*(mu - ar) )
    else:
      # High traffic assumption
      f_c = (n - 2)/(n - 1)
      f_p = 1 - f_c
      EV_p = 1/((n - 1)*mu)
      EV_c = 1/mu * (H(n) - H(n-2) )
      EV = f_p*EV_p + f_c*EV_c
      
      EV2_p = 2 / ((n-1)*mu)**2
      EV2_c = 1/mu**2 * (H_2(n) - H_2(n-2) ) + EV_c**2
      EV2 = f_p*EV2_p + f_c*EV2_c
      # ET = EV/(1 - ar*EV) # w/ M/M/1 assumption
      ET = EV + ar/2 * EV2/(1 - ar*EV) # w/ M/G/1 assumption
      return ET
  else:
    f_c = (n - 2)/(n - 1)
    f_p = 1 - f_c
    EV_p = EXm_n_k(1, n-1, 1, dist_m)
    EV_c = EXm_n_k(1, n, 2, dist_m)
    EV = f_p*EV_p + f_c*EV_c
    
    EV2_p = EXm_n_k(2, n-1, 1, dist_m)
    EV2_c = EXm_n_k(2, n, 2, dist_m)
    EV2 = f_p*EV2_p + f_c*EV2_c
    # ET = EV/(1 - ar*EV) # w/ M/M/1 assumption
    ET = EV + ar/2 * EV2/(1 - ar*EV) # w/ M/G/1 assumption
    return ET

if __name__ == "__main__":
  # plot_moments_for_types()
  
  dist_m = {'dist': 'Exp', 'mu': 1}
  n = 5
  print("n= {}".format(n) )
  for k in range(1, 2):
    pd = pdeparted(n, k, dist_m)
    print("k= {}, pd= {}".format(k, pd) )
    
  