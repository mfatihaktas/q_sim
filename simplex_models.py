import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy, simpy, getopt
from math import factorial
from numpy import linalg
from patch import *

def harmonic_sum(n):
  sum_ = 0
  for i in range(1, n+1):
    sum_ += float(1/i)
  return sum_

def harmonic_2_sum(n):
  sum_ = 0
  for i in range(1, n+1):
    sum_ += float(1/(i**2) )
  return sum_

def gen_harmonic_sum(n, k):
  sum_ = 0
  for i in range(1, n+1):
    if (i - k) == 0:
      continue
    sum_ += float(1/(i*(i - k) ) )
  return sum_

def fj_2_2_E_T(arr_rate, mu):
  ro = arr_rate/mu
  return (12 - ro)/(8*(mu - arr_rate) )

# -----------------------------------------------  MDS  ------------------------------------------ #
def mds_n_2_E_T(arr_rate, mu, n): # for only k = 2
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
  
def adjustable_mds_n_2_E_T(arr_rate, mu, n): # for only k = 2
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

def split_merge_mds_n_k_E_T(arr_rate, mu, n, k): # works for k >= 1
  harmonic_diff = harmonic_sum(n) - harmonic_sum(n-k)
  if 1/arr_rate < harmonic_diff/mu:
    return None
  ro = arr_rate/mu
  return harmonic_diff/mu + \
    arr_rate*(harmonic_2_sum(n)-harmonic_2_sum(n-k) + harmonic_diff**2)/(2*mu**2 * (1 - ro*harmonic_diff) )

def adj_split_merge_mds_n_k_E_T(arr_rate, mu, n, k):
  return split_merge_mds_n_k_E_T(arr_rate, mu, n, k-1) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)

def adj_2_split_merge_mds_n_k_E_T(arr_rate, mu, n, k):
  return split_merge_mds_n_k_E_T(arr_rate, mu, n, k-2) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+2, 1) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)

def recur_split_merge_mds_n_k_E_T(arr_rate, mu, n, k):
  if n > k:
    if k > 1:
      return recur_split_merge_mds_n_k_E_T(arr_rate, mu, n, k-1) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)
    elif k == 1:
      return split_merge_mds_n_k_E_T(arr_rate, mu, n, k)
    else:
      log(ERROR, "Unexpected k= {}".format(k) )
      sys.exit(1)
  else:
    log(ERROR, "Unexpected n= {} <= k= {}".format(n, k) )
    return split_merge_mds_n_k_E_T(arr_rate, mu, n, k)

def mds_inner_bound_on_arr_rate(n, k, mu):
  sm_mds_n_k_E_S = 1/mu * (harmonic_sum(n) - harmonic_sum(n-k) )
  return 1/sm_mds_n_k_E_S
  
def simplex_inner_bound_on_arr_rate(r, t, mu):
  def beta(x, y):
    return math.gamma(x)*math.gamma(y)/math.gamma(x+y)
  sm_simplex_E_S = 1/(mu*r)*beta(t+1, 1/r)
  return 1/sm_simplex_E_S
  
# -----------------------------------------  Simplex 1 repair  ----------------------------------- #
def simplex_w_one_repair__sys_time(arr_rate, mu, c=None):
  if c == None:
    E_S_c = 0.665/mu
    E_S_c_2 = 7/(9*(mu**2) )
    E_S_p = 0.5/mu
    E_S_p_2 = 0.5/(mu**2)
    E_S = 0.6*E_S_c + 0.4*E_S_p
    E_S_2 = 0.6*E_S_c_2 + 0.4*E_S_p_2
    # E_T = E_S/(1 - arr_rate*E_S) # w/ M/M/1 assumption
    E_T = E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S) # w/ M/G/1 assumption
    return E_T
  else:
    C = 3*mu
    f_c = 1/(1 + c/2/(1+c) )
    mu = C/2/(c+1)
    m = 0.165
    E_S = 1/mu * (m*f_c + 0.5)
    n = 5/18
    E_S_2 = 1/(mu**2) * (n*f_c + 0.5)
    return E_S + arr_rate*E_S_2/2/(1-arr_rate*E_S)

def simplex_w_one_repair__sys_time_trial(c, t, arr_rate, mu):
  # mu_1 = mu - arr_rate
  # mu_2 = 2*mu - arr_rate
  gamma = mu
  # return 1/gamma + 1/mu * (mu_2/mu_1 - mu_1/mu_2 + mu_1/(gamma+mu_2) - mu_2/(gamma+mu_1) )
  return 1/(c*mu + gamma - arr_rate)

# -----------------------------------------  Simplex 2 repair  ----------------------------------- #
def simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu):
  gamma = mu
  nu = 4*mu + gamma
  P = []
  state_list = []
  state_index_map = {}
  # mc_truncation_state_id = 4
  for y in range(mc_truncation_state_id+1):
    for x in range(mc_truncation_state_id+1):
      s = str(x)+"_"+str(y)
      state_list.append(s)
      state_index_map[s] = len(state_list)-1
  log(WARNING, "state_list=\n {}".format(pprint.pformat(state_list) ) )
  log(WARNING, "state_index_map=\n {}".format(pprint.pformat(state_index_map) ) )
  for i in range(len(state_list) ):
    P.append(len(state_list)*[0] )
  for y in range(mc_truncation_state_id+1):
    for x in range(mc_truncation_state_id+1):
      s = str(x)+"_"+str(y)
      if x == 0 and y == 0:
        P[state_index_map[s] ][state_index_map[s] ] = gamma/nu
        P[state_index_map[s] ][state_index_map["0_1"] ] = 2*mu/nu
        P[state_index_map[s] ][state_index_map["1_0"] ] = 2*mu/nu
      elif x == 0:
        P[state_index_map[s] ][state_index_map[str(x)+"_"+str(y-1)] ] = (mu+gamma)/nu
        P[state_index_map[s] ][state_index_map[str(x+1)+"_"+str(y)] ] = 2*mu/nu
        if (y == mc_truncation_state_id):
          P[state_index_map[s] ][state_index_map[s] ] = mu/nu
        else:
          P[state_index_map[s] ][state_index_map[str(x)+"_"+str(y+1)] ] = mu/nu
      elif y == 0:
        P[state_index_map[s] ][state_index_map[str(x-1)+"_"+str(y)] ] = (mu+gamma)/nu
        P[state_index_map[s] ][state_index_map[str(x)+"_"+str(y+1)] ] = 2*mu/nu
        if (x == mc_truncation_state_id):
          P[state_index_map[s] ][state_index_map[s] ] = mu/nu
        else:
          P[state_index_map[s] ][state_index_map[str(x+1)+"_"+str(y)] ] = mu/nu
      else:
        P[state_index_map[s] ][state_index_map[str(x-1)+"_"+str(y-1)] ] += (2*mu+gamma)/nu
        if (x == mc_truncation_state_id):
          P[state_index_map[s] ][state_index_map[s] ] += mu/nu
        else:
          P[state_index_map[s] ][state_index_map[str(x+1)+"_"+str(y)] ] += mu/nu
        if (y == mc_truncation_state_id):
          P[state_index_map[s] ][state_index_map[s] ] += mu/nu
        else:
          P[state_index_map[s] ][state_index_map[str(x)+"_"+str(y+1)] ] += mu/nu
  
  log(WARNING, "P= \n{}".format(numpy.matrix(P) ) )
  pi_v = linalg.matrix_power(P, 100)[0]
  log(WARNING, "pi_v= {}".format(pprint.pformat(pi_v) ) )
  state_prob_map = {state_list[i]:pi for i,pi in enumerate(pi_v) }
  log(WARNING, "state_prob_map= {}".format(pprint.pformat(state_prob_map) ) )
  """
  for i in range(len(state_list) ):
    P[i] = len(state_list)*[0]
  P[state_index_map["0_0"] ][state_index_map["0_0"] ] = gamma/nu
  P[state_index_map["0_0"] ][state_index_map["0_1"] ] = 2*mu/nu
  P[state_index_map["0_0"] ][state_index_map["1_0"] ] = 2*mu/nu
  P[state_index_map["0_1"] ][state_index_map["0_0"] ] = (mu+gamma)/nu
  P[state_index_map["0_1"] ][state_index_map["0_2"] ] = mu/nu
  P[state_index_map["0_1"] ][state_index_map["1_1"] ] = 2*mu/nu
  P[state_index_map["0_2"] ][state_index_map["0_1"] ] = (mu+gamma)/nu
  P[state_index_map["0_2"] ][state_index_map["1_2"] ] = 2*mu/nu
  P[state_index_map["0_2"] ][state_index_map["0_2"] ] = mu/nu
  # P[state_index_map["0_3"] ][state_index_map["0_2"] ] = (mu+gamma)/nu
  P[state_index_map["1_0"] ][state_index_map["0_0"] ] = (mu+gamma)/nu
  P[state_index_map["1_0"] ][state_index_map["2_0"] ] = mu/nu
  P[state_index_map["1_0"] ][state_index_map["1_1"] ] = 2*mu/nu
  P[state_index_map["1_1"] ][state_index_map["0_0"] ] = (2*mu+gamma)/nu
  P[state_index_map["1_1"] ][state_index_map["2_1"] ] = mu/nu
  P[state_index_map["1_1"] ][state_index_map["1_2"] ] = mu/nu
  P[state_index_map["1_2"] ][state_index_map["0_1"] ] = (2*mu+gamma)/nu
  P[state_index_map["1_2"] ][state_index_map["2_2"] ] = mu/nu
  P[state_index_map["1_2"] ][state_index_map["1_2"] ] = mu/nu
  P[state_index_map["2_0"] ][state_index_map["1_0"] ] = (mu+gamma)/nu
  P[state_index_map["2_0"] ][state_index_map["2_1"] ] = 2*mu/nu
  P[state_index_map["2_0"] ][state_index_map["2_0"] ] = mu/nu
  P[state_index_map["2_1"] ][state_index_map["1_0"] ] = (2*mu+gamma)/nu
  P[state_index_map["2_1"] ][state_index_map["2_2"] ] = mu/nu
  P[state_index_map["2_1"] ][state_index_map["2_1"] ] = mu/nu
  P[state_index_map["2_2"] ][state_index_map["1_1"] ] = (2*mu+gamma)/nu
  P[state_index_map["2_2"] ][state_index_map["2_2"] ] = 2*mu/nu
  log(WARNING, "LATER_P= \n{}".format(numpy.matrix(P) ) )
  # for p in numpy.arange(10, 100, 10):
  #   P_p = linalg.matrix_power(P, p)
  #   log(WARNING, "p= {}, P_p=\n{}".format(p, pprint.pformat(P_p) ) )
  pi_v = linalg.matrix_power(P, 100)[0]
  log(WARNING, "LATER pi_v= {}".format(pprint.pformat(pi_v) ) )
  state_prob_map = {state_list[i]:pi for i,pi in enumerate(pi_v) }
  log(WARNING, "LATER state_prob_map= {}".format(pprint.pformat(state_prob_map) ) )
  """
  return state_prob_map

def binomial(x, y):
  try:
    binom = factorial(x) // factorial(y) // factorial(x - y)
  except ValueError:
    binom = 0
  return binom

def avq_low_traff_serv_time_first_moment(r, t, mu):
  def beta(x, y):
    return math.gamma(x)*math.gamma(y)/math.gamma(x+y)
  return 1/(mu*r)*beta(t+1, 1/r)

def avq_low_traff_serv_time_second_moment(r, t, mu):
  second_moment = 0
  for j in range(t + 1):
    inner_term = 0
    for l in range(r*j + 1):
      inner_term += (-1)**l * binomial(r*j, l)*(2/(mu**2 * (l+1)**2) )
    second_moment += binomial(t, j) * (-1)**j * inner_term
  return second_moment

def simplex_w_two_repair__sys_time(arr_rate, mu):
  E_S_s_c = avq_low_traff_serv_time_first_moment(2, 2, mu)
  E_S_s_c_2 = avq_low_traff_serv_time_second_moment(2, 2, mu)
  E_S_p_1st_kind = 5/(12*mu)
  E_S_p_1st_kind_2 = 23/(72*mu**2)
  E_S_p_2nd_kind = 1/(3*mu)
  E_S_p_2nd_kind_2 = 2/(9*mu**2)
  
  mc_truncation_state_id = 5
  state_prob_map = simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu)
  log(WARNING, "state_prob_map=\n {}".format(pprint.pformat(state_prob_map) ) )
  
  gamma = mu
  nu = 4*mu + gamma
  state_sum = sum([state_prob_map["0_"+str(i)] + state_prob_map[str(i)+"_0"] for i in range(1, mc_truncation_state_id+1, 1) ] )
  f_jd = state_prob_map["0_0"]*gamma/nu \
       + (mu+gamma)/nu * state_sum \
       + (2*mu+gamma)/nu * (1 - state_prob_map["0_0"] - state_sum)
  f_jc = gamma/nu * state_prob_map["0_0"] \
       + (mu+gamma)/nu * (state_prob_map["0_1"] + state_prob_map["1_0"] ) \
       + (2*mu+gamma)/nu * state_prob_map["1_1"]
  f_jp1 = (2*mu+gamma)/nu * sum([state_prob_map[str(i)+"_1"] for i in range(2, mc_truncation_state_id+1, 1) ] ) \
        + (2*mu+gamma)/nu * sum([state_prob_map["1_"+str(i)] for i in range(2, mc_truncation_state_id+1, 1) ] )
  f_c = f_jc/f_jd
  f_p1 = f_jp1/f_jd
  f_p2 = 1 - f_c - f_p1
  log(WARNING, "f_c= {}, f_p1= {}, f_p2= {}".format(f_c, f_p1, f_p2) )
  
  E_S = f_c*E_S_s_c + f_p1*E_S_p_1st_kind + f_p2*E_S_p_2nd_kind
  E_S_2 = f_c*E_S_s_c_2 + f_p1*E_S_p_1st_kind_2 + f_p2*E_S_p_2nd_kind_2
  E_T = E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  return E_T

def simplex_split_merge_sys_time(t, arr_rate, mu):
  E_S = avq_low_traff_serv_time_first_moment(2, t, mu)
  E_S_2 = avq_low_traff_serv_time_second_moment(2, t, mu)
  return E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)

def simplex_w_one_repair__parametric_sys_time():
  C = 3
  def parametric_sys_time(arr_rate, r):
    f_c = 1/(1 + r/2/(1+r) )
    mu = C/2/(r+1)
    m = 0.165
    E_S = 1/mu * (m*f_c + 0.5)
    n = 5/18
    E_S_2 = 1/(mu**2) * (n*f_c + 0.5)
    return E_S + arr_rate*E_S_2/2/(1-arr_rate*E_S)
  color = iter(cm.rainbow(numpy.linspace(0, 1, 20) ) )
  for arr_rate in numpy.arange(0.05, 1.0, 0.1):
    c_list = numpy.arange(0, 1, 0.05)
    E_T_list = [parametric_sys_time(arr_rate, r) for r in c_list]
    plot.plot(c_list, E_T_list, 'o', color=next(color), label=r'$\lambda$={0:0.2f}'.format(arr_rate) )
    plot.legend()
  plot.xlabel("r")
  plot.ylabel("E[T]")
  plot.title(r'C= {}'.format(C) )
  plot.savefig("simplex_w_one_repair__parametric_sys_time.png")

if __name__ == "__main__":
  # simplex_w_two_repair__sys_time(0.9, 1.0)
  simplex_w_one_repair__parametric_sys_time()
  # mu = 1.0
  # mc_truncation_state_id__state_prob_map_map = {}
  # for mc_truncation_state_id in range(2, 10, 1):
  #   mc_truncation_state_id__state_prob_map_map[mc_truncation_state_id] = \
  #     simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu)
  # log(WARNING, "mc_truncation_state_id__state_prob_map_map= {}".format(pprint.pformat(mc_truncation_state_id__state_prob_map_map) ) )
  