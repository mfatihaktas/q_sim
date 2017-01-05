import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy, sympy, simpy, getopt
from math import factorial
from numpy import linalg
from patch import *
from sim_components import *

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

def E_T_fj_2(arr_rate, mu):
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

def sm_mds_n_k_E_T(arr_rate, mu, n, k): # works for k >= 1
  harmonic_diff = harmonic_sum(n) - harmonic_sum(n-k)
  if 1/arr_rate < harmonic_diff/mu:
    return None
  ro = arr_rate/mu
  return harmonic_diff/mu + \
    arr_rate*(harmonic_2_sum(n)-harmonic_2_sum(n-k) + harmonic_diff**2)/(2*mu**2 * (1 - ro*harmonic_diff) )

def adj_sm_mds_n_k_E_T(arr_rate, mu, n, k):
  return sm_mds_n_k_E_T(arr_rate, mu, n, k-1) + sm_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)

def adj_2_sm_mds_n_k_E_T(arr_rate, mu, n, k):
  return sm_mds_n_k_E_T(arr_rate, mu, n, k-2) + sm_mds_n_k_E_T(arr_rate, mu, n-k+2, 1) + sm_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)

def recur_sm_mds_n_k_E_T(arr_rate, mu, n, k):
  if n > k:
    if k > 1:
      return recur_sm_mds_n_k_E_T(arr_rate, mu, n, k-1) + sm_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)
    elif k == 1:
      return sm_mds_n_k_E_T(arr_rate, mu, n, k)
    else:
      log(ERROR, "Unexpected k= {}".format(k) )
      sys.exit(1)
  else:
    log(ERROR, "Unexpected n= {} <= k= {}".format(n, k) )
    return sm_mds_n_k_E_T(arr_rate, mu, n, k)

def mds_inner_bound_on_arr_rate(n, k, mu):
  sm_mds_n_k_E_S = 1/mu * (harmonic_sum(n) - harmonic_sum(n-k) )
  return 1/sm_mds_n_k_E_S
  
def simplex_inner_bound_on_arr_rate(r, t, mu, w_sys=True):
  E_S = 1
  if w_sys:
    def beta(x, y):
      return math.gamma(x)*math.gamma(y)/math.gamma(x+y)
    E_S = 1/(mu*r)*beta(t+1, 1/r)
  else:
    E_S = 1/mu * sum([binomial(t,i) * 2**i*(-1)**(t-i)/(2*t-i) for i in range(t+1) ] )
  
  return 1/E_S

# -----------------------------------  Simplex w/ split-to-one  ------------------------------- #
def arr_rate_ub_simplex_split_to_one(t, mu):
  gamma = mu
  split_prob_l = [1/(t+1) for g in range(t+1) ]
  
  return min(gamma/split_prob_l[0], mu/max(split_prob_l[1:t+1] ) )

def E_T_simplex_split_to_one(t, arr_rate, mu, p_r=None):
  gamma = mu
  split_prob_l = [0]*(t+1)
  for g in range(t+1):
    if p_r is None:
      split_prob_l[g] = 1/(t+1)
    else:
      if g == 0:
        split_prob_l[g] = 1-t*p_r
      else:
        split_prob_l[g] = p_r
  
  E_T = 0
  for g in range(t+1):
    arr_rate_ = arr_rate*split_prob_l[g]
    if g == 0:
      E_T += split_prob_l[g] * 1/(gamma-arr_rate_)
    else:
      E_T += split_prob_l[g] * E_T_fj_2(arr_rate_, mu)
  return E_T

def plot_E_T_simplex_split_to_one(t, mu):
  arr_rate_ub = simplex_inner_bound_on_arr_rate(2, t, mu, True)
  
  color = iter(cm.rainbow(numpy.linspace(0, 1, 4) ) )
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/4):
    p_r_l, E_T_l = [], []
    for p_r in numpy.arange(0.05, 1, 1/20):
      p_r_l.append(p_r)
      E_T_l.append(max(0, E_T_simplex_split_to_one(t, arr_rate, mu, p_r) ) )
    plot.plot(p_r_l, E_T_l, 'o', color=next(color), label=r't= {}, $\lambda$= {}, $\mu$= {}'.format(t, arr_rate, mu) )
    plot.legend()
  plot.xlabel(r'$p_r$')
  plot.ylabel(r'$E[T]$')
  # plot.title()
  plot.savefig("plot_E_T_simplex_split_to_one.png")
# -----------------------------------------  Simplex(t=1)  ----------------------------------- #
def simplex_w_one_repair__E_T(arr_rate, mu, c=None):
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
    Cap = 3*mu
    f_c = 1/(1 + 2/c/(c+2) )
    # mu = Cap/(c+2)
    E_S_c = 2*(c+2)/Cap/(c+1) - 1/Cap
    E_S_c_2 = (2*(c+2)/Cap/(c+1) )**2 - 2/(Cap**2)
    E_S_p = (c+2)/Cap/(c+1)
    E_S_p_2 = 2*((c+2)/Cap/(c+1) )**2
    
    E_S = f_c*E_S_c + (1-f_c)*E_S_p
    E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
    return E_S + arr_rate*E_S_2/2/(1-arr_rate*E_S)

def simplex_w_one_repair__E_T_matrix_analytic(t, arr_rate, mu):
  l = arr_rate
  a = b = g = mu
  d = a + b + g + l
  F_0 = numpy.matrix([[-l , 0  , l , 0  ], \
                      [a+g, b-d, 0 , 0  ], \
                      [g  , b  , -d, a  ], \
                      [b+g, 0  , 0 , a-d] ] )
  H_0 = numpy.matrix([[0, 0, 0, 0, 0], \
                      [0, l, 0, 0, 0], \
                      [0, 0, l, 0, 0], \
                      [0, 0, 0, l, 0] ] )
  L_0 = numpy.matrix([[0, a+g, 0  , 0  ], \
                      [0, 0  , a+g, 0  ], \
                      [0, 0  , g  , 0  ], \
                      [0, 0  , b+g, 0  ], \
                      [0, 0  , 0  , b+g] ] )
  F = numpy.matrix([[b-d, 0 , 0 , 0 , 0  ], \
                    [b  , -d, 0 , 0 , 0  ], \
                    [0  , b , -d, a , 0  ], \
                    [0  , 0 , 0 , -d, a  ], \
                    [0  , 0 , 0 , 0 , a-d] ] )
  L = numpy.matrix([[0, a+g, 0  , 0  , 0], \
                    [0, 0  , a+g, 0  , 0], \
                    [0, 0  , g  , 0  , 0], \
                    [0, 0  , b+g, 0  , 0], \
                    [0, 0  , 0  , b+g, 0] ] )
  H = numpy.matrix([[l, 0, 0, 0, 0], \
                    [0, l, 0, 0, 0], \
                    [0, 0, l, 0, 0], \
                    [0, 0, 0, l, 0], \
                    [0, 0, 0, 0, l] ] )
  def R(e):
    R = numpy.zeros((5, 5) )
    done = False
    counter = 0
    while not done:
      _R = R
      R = -(R**2 * L + H)*F**-1
      
      d = 0
      for r in range(5):
        for c in range(5):
          d_ = abs(R[r, c] - _R[r, c] )
          if d_ > d:
            d = d_
      if d < e:
        break
      counter = counter + 1
    # print("R:: counter= {}".format(counter) )
    return R
  
  R = R(0.000001)
  # print("R=\n {}".format(R) )
  FI_u = numpy.concatenate((F_0, H_0), axis=1)
  FI_l = numpy.concatenate((L_0, R*L+F), axis=1)
  FI = numpy.concatenate((FI_u, FI_l), axis=0)
  # print("FI=\n {}".format(FI) )
  c_v = numpy.concatenate((numpy.ones((4,1)), (numpy.identity(5)-R)**-1 * numpy.ones((5,1)) ), axis=0)
  # print("c_v=\n {}".format(c_v) )
  FI[:, 0] = c_v
  # print("FI=\n {}".format(FI) )
  PI_v = numpy.matrix([[1, 0, 0, 0, 0, 0, 0, 0, 0] ]) * FI**-1
  # print("PI_v=\n {}".format(PI_v) )
  PI_0_v = PI_v[:, 0:4]
  PI_1_v = PI_v[:, 4:9]
  print("PI_0_v= {}, PI_1_v= {}".format(PI_0_v, PI_1_v) )
  # for i in range(2, 10):
  #   PI_i_v = PI_1_v*R**(i-1)
  #   print("i= {}, PI_i_v= {}".format(i, PI_i_v) )
  # N_prob_map = {}
  # N_prob_map[0] = PI_0_v[0, 0]
  # N_prob_map[1] = numpy.sum(PI_0_v[:, 1:4] )
  # N_prob_map[2] = numpy.sum(PI_1_v)
  # for n in range(3, 100):
  #   PI_i_v = PI_1_v*R**(n-1)
  #   N_prob_map[n] = numpy.sum(PI_i_v)
  # # print("N_prob_map= {}".format(pprint.pformat(N_prob_map) ) )
  # return sum([n*prob for n,prob in N_prob_map.items() ] )/arr_rate
  
  E_N = PI_0_v*numpy.ones((4,1)) - PI_0_v[0, 0] + \
        PI_1_v*((numpy.identity(5)-R)**-2 + (numpy.identity(5)-R)**-1)*numpy.ones((5,1))
  print("E_N= {}".format(E_N) )
  return E_N[0, 0]/arr_rate

# -----------------------------------------  Simplex 2 repair  ----------------------------------- #
def simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu):
  gamma = mu
  nu = gamma + 4*mu
  P = []
  state_l = []
  state_index_map = {}
  # mc_truncation_state_id = 4
  for y in range(mc_truncation_state_id+1):
    for x in range(mc_truncation_state_id+1):
      s = str(x)+"_"+str(y)
      state_l.append(s)
      state_index_map[s] = len(state_l)-1
  # log(WARNING, "state_l=\n {}".format(pprint.pformat(state_l) ) )
  # log(WARNING, "state_index_map=\n {}".format(pprint.pformat(state_index_map) ) )
  for i in range(len(state_l) ):
    P.append(len(state_l)*[0] )
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
  
  # log(WARNING, "P= \n{}".format(numpy.matrix(P) ) )
  pi_v = linalg.matrix_power(P, 100)[0]
  # log(WARNING, "pi_v= {}".format(pprint.pformat(pi_v) ) )
  state_prob_map = {state_l[i]:pi for i,pi in enumerate(pi_v) }
  # log(WARNING, "state_prob_map= {}".format(pprint.pformat(state_prob_map) ) )
  """
  for i in range(len(state_l) ):
    P[i] = len(state_l)*[0]
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
  state_prob_map = {state_l[i]:pi for i,pi in enumerate(pi_v) }
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

def simplex_w_two_repair__E_T(arr_rate, mu, M):
  E_S_c = avq_low_traff_serv_time_first_moment(2, 2, mu)
  E_S_c_2 = avq_low_traff_serv_time_second_moment(2, 2, mu)
  E_S_p_1st = 5/(12*mu)
  E_S_p_1st_2 = 23/(72*mu**2)
  E_S_p_2nd = 1/(3*mu)
  E_S_p_2nd_2 = 2/(9*mu**2)
  
  mc_truncation_state_id = M # 5
  state_prob_map = simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu)
  log(WARNING, "state_prob_map=\n {}".format(pprint.pformat(state_prob_map) ) )
  
  gamma = mu
  nu = gamma + 4*mu
  state_sum = sum([state_prob_map["0_"+str(i)] + state_prob_map[str(i)+"_0"] for i in range(1, mc_truncation_state_id+1, 1) ] )
  f_jd = state_prob_map["0_0"]*gamma/nu \
       + (mu+gamma)/nu * state_sum \
       + (2*mu+gamma)/nu * (1 - state_prob_map["0_0"] - state_sum)
  f_jc = gamma/nu * state_prob_map["0_0"] \
       + (mu+gamma)/nu * (state_prob_map["0_1"] + state_prob_map["1_0"] ) \
       + (2*mu+gamma)/nu * state_prob_map["1_1"]
  f_jp1 = (mu+gamma)/nu * sum([state_prob_map[str(i)+"_0"]+state_prob_map["0_"+str(i)] for i in range(2, mc_truncation_state_id+1, 1) ] ) \
        + (2*mu+gamma)/nu * sum([state_prob_map[str(i)+"_1"]+state_prob_map["1_"+str(i)] for i in range(2, mc_truncation_state_id+1, 1) ] )
  f_c = f_jc/f_jd
  f_p1 = f_jp1/f_jd
  f_p2 = 1 - f_c - f_p1
  log(WARNING, "f_c= {}, f_p1= {}, f_p2= {}".format(f_c, f_p1, f_p2) )
  
  E_S = f_c*E_S_c + f_p1*E_S_p_1st + f_p2*E_S_p_2nd
  E_S_2 = f_c*E_S_c_2 + f_p1*E_S_p_1st_2 + f_p2*E_S_p_2nd_2
  E_T = E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  return E_T

def simplex_wo_sys_w_two_repair__E_T(arr_rate, mu):
  return 1/(4*mu-arr_rate) + 1/(3*mu-arr_rate) + 2/3/(2*mu-arr_rate)

def simplex_sm_E_T(t, arr_rate, mu, c=None):
  if c == None:
    E_S = avq_low_traff_serv_time_first_moment(2, t, mu)
    E_S_2 = avq_low_traff_serv_time_second_moment(2, t, mu)
    return E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  else:
    if t == 1:
      Cap = 3*mu
      E_S_c = 2*(c+2)/Cap/(c+1) - 1/Cap
      E_S_c_2 = (2*(c+2)/Cap/(c+1) )**2 - 2/(Cap**2)
      return E_S_c + (arr_rate/2)*E_S_c_2/(1 - arr_rate*E_S_c)
    else:
      log(ERROR, "NOT ready!; t= {}".format(t) )
      return 1

def simplex_wo_sys_sm_E_T(t, arr_rate, mu, c=None):
  E_S = 1/mu * sum([binomial(t,i) * 2**i*(-1)**(t-i)/(2*t-i) for i in range(t+1) ] )
  E_S_2 = 2/mu**2 * sum([binomial(t,i) * 2**i*(-1)**(t-i)/(2*t-i) for i in range(t+1) ] )
  return E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  
def simplex_w_one_repair__parametric_E_T():
  Cap = 3
  def parametric_E_T(arr_rate, c):
    f_c = 1/(1 + 2/c/(c+2) )
    # mu = Cap/(c+2)
    E_S_c = 2*(c+2)/Cap/(c+1) - 1/Cap
    E_S_c_2 = (2*(c+2)/Cap/(c+1) )**2 - 2/(Cap**2)
    E_S_p = (c+2)/Cap/(c+1)
    E_S_p_2 = 2*((c+2)/Cap/(c+1) )**2
    
    E_S = f_c*E_S_c + (1-f_c)*E_S_p
    E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
    return E_S + arr_rate*E_S_2/2/(1-arr_rate*E_S)
  color = iter(cm.rainbow(numpy.linspace(0, 1, 20) ) )
  for arr_rate in numpy.arange(0.05, 1.0, 0.1):
    c_l = numpy.arange(0.05, 5, 0.05)
    E_T_l = [parametric_E_T(arr_rate, c) for c in c_l]
    plot.plot(c_l, E_T_l, '.', color=next(color), label=r'$\lambda$={0:0.2f}'.format(arr_rate) )
    plot.legend()
  plot.xlabel("c")
  plot.ylabel("E[T]")
  plot.title(r'Total service rate= {}'.format(Cap) )
  plot.savefig("simplex_w_one_repair__parametric_E_T.png")

def simplex_steady_state_prob_hist():
  k, r, t = 2, 2, 1
  num_q = int(1 + t*r)
  qid_l = ["{}".format(i) for i in range(1, num_q + 1) ]
  qmu_l = [1, 1, 1]
  def get_state_prob_map(arr_rate):
    log(WARNING, "arr_rate= {}, k= {}, r= {}, t= {}, qmu_l= {}".format(arr_rate, k, r, t, pprint.pformat(qmu_l) ) )
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    a_q = AVQ("a_q", env, k, r, t, qid_l, qserv_rate_l=qmu_l)
    aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 0.1)
    a_q.join_q.out_m = aq_monitor
    pg.out = a_q
    env.run(until=50000)
    
    # print("aq_monitor.polled_state__counter_map= {}".format(pprint.pformat(aq_monitor.polled_state__counter_map) ) )
    total_counter = sum([c for rs, c in aq_monitor.polled_state__counter_map.items() ] )
    state_prob_map = {rs:float(c)/total_counter for rs, c in aq_monitor.polled_state__counter_map.items() }
    # print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
    return state_prob_map # ['0,(0,0)']
  # for arr_rate in numpy.arange(0.05, 1.2, 0.1):
  color = iter(cm.rainbow(numpy.linspace(0, 1, 20) ) )
  plot.figure(figsize=(20,10) )
  for arr_rate in numpy.arange(0.05, 1.3, 0.1):
  # for arr_rate in numpy.arange(0.05, 0.1, 0.1):
    state_prob_map = get_state_prob_map(arr_rate)
    
    def state(kp, i, j):
      return "{},({},{})".format(kp, i, j)
    i__tau_l_map = {}
    for i in range(10):
      if i not in i__tau_l_map:
        i__tau_l_map[i] = []
      for kp in range(i, 10):
        s_u, s_l = state(kp, i, 0), state(kp+1, i, 0)
        if s_u in state_prob_map and s_l in state_prob_map:
          i__tau_l_map[i].append(state_prob_map[s_l]/state_prob_map[s_u] )
        # if state(k+1, 0, i) in state_prob_map:
        #   i__tau_l_map[i].append(state_prob_map[state(k+1, 0, i) ] /state_prob_map[state(k, 0, i) ] )
    log(WARNING, "i__tau_l_map=\n {}".format(pprint.pformat(i__tau_l_map) ) )
    #
    wing_cutoff_i = 2
    wing_cutoff_sum = 0
    for s, p in state_prob_map.items():
      split_l = s.split(",")
      if int(split_l[1].split("(")[1] ) > wing_cutoff_i or int(split_l[2].split(")")[0] ) > wing_cutoff_i:
        wing_cutoff_sum += p
      
    s_l, p_l = [], []
    for s, p in state_prob_map.items():
      if p > 0.01:
        s_l.append(s)
        p_l.append(p)
    plot.bar(range(len(p_l) ), p_l, color=next(color) )
    plot.xticks([i+0.5 for i in range(len(s_l) ) ], s_l, size='small')
    plot.xlabel("State")
    plot.ylabel("Steady-state probability")
    plot.title(r't= {}, $\lambda$= {}, [$\alpha$, $\beta$, $\gamma$]= {}, sum_on_plot= {}, wing_cutoff_sum= {}'. \
      format(t, "{0:.2f}".format(arr_rate), pprint.pformat(qmu_l), "{0:.2f}".format(sum(p_l)), "{0:.2f}".format(wing_cutoff_sum) ) )
    plot.savefig("simplex_steady_state_prob_hist_ar_{0:.2f}.png".format(arr_rate) )
    plot.clf()

# ##############################  UB; Simplex(t=2:[sys, MDS(r,2)])  ############################# #
def E_S_c_avq_sys__mds_r_2(gamma, mu, r):
  return r/(gamma+(r-1)*mu) - (r-1)/(gamma+r*mu)

def E_S_c_2_avq_sys__mds_r_2(gamma, mu, r):
  return 2*r/(gamma+(r-1)*mu)**2 - 2*(r-1)/(gamma+r*mu)**2

def E_T_avq_sys__mds_r_2(arr_rate, gamma, mu, r):
  # This turns out to be a worse upper bound than split-merge
  """
  f_c = mu*(gamma+(r-2)*mu)/(r*(r-1)*mu**2 + gamma*(gamma+2*(r-1)*mu) )
  
  E_S_c = E_S_c_avq_sys__mds_r_2(gamma, mu, r)
  E_S_c_2 = E_S_c_2_avq_sys__mds_r_2(gamma, mu, r)
  E_S_p = E_S_c_avq_sys__mds_r_2(gamma, mu, r-1)
  E_S_p_2 = E_S_c_2_avq_sys__mds_r_2(gamma, mu, r-1)
  log(WARNING, "E_S_c= {}, E_S_c_2= {}, E_S_p= {}, E_S_p_2= {}".format(E_S_c, E_S_c_2, E_S_p, E_S_p_2) )
  
  E_S = f_c*E_S_c + (1-f_c)*E_S_p
  E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
  E_T = E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  """
  # After one of the r servers goes ahead by one, the next server for job termination in the repair group
  # is held fixed
  ro = mu/(mu+gamma)
  p_0 = (1-ro)/(1+(r-1)*ro)
  def p_i(i):
    return r*ro**i * p_0
  def pi_i(i):
    sum_p_i_nu_i = p_0*mu*(r-2) + 2*mu+gamma
    if i == 0:
      return p_0*(gamma+r*mu)/sum_p_i_nu_i
    else:
      return p_i(i)*(2*mu+gamma)/sum_p_i_nu_i
  f_jd = pi_i(0)*gamma/(gamma+r*mu) + (1-pi_i(0) )*(mu+gamma)/(2*mu+gamma)
  f_jc = pi_i(0)*gamma/(gamma+r*mu) + pi_i(1)*(mu+gamma)/(2*mu+gamma)
  f_c = f_jc/f_jd
  
  beta = r*mu/(gamma+r*mu)
  E_A = 1/(gamma+r*mu)
  E_A_2 = 2/(gamma+r*mu)**2
  E_B = 1/(gamma+mu)
  E_B_2 = 1/(gamma+mu)**2
  
  E_S_c = E_A + beta*E_B
  E_S_c_2 = E_A_2 + beta**2*E_B_2 + 2*beta*E_A*E_B
  E_S_p = 1/(gamma+mu)
  E_S_p_2 = 2/(gamma+mu)**2
  
  E_S = f_c*E_S_c + (1-f_c)*E_S_p
  E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
  E_T = E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  return E_T

# ########################################  LB; Simplex(t)  ###################################### #
def E_T_simplex_varki_gauri_lb(t, arr_rate, gamma, mu):
  def P_i_next(i):
    return 2*(t-i)*mu/(gamma+(2*t-i)*mu)
  def mu_i(i):
    return gamma + (2*t-i)*mu
  def lambda_i(i):
    if i == 0:
      return arr_rate
    else:
      prod = 1
      for k in range(i):
        prod = prod * P_i_next(k)
      return arr_rate*prod
  E_T = 0
  # for i in range(2):
  for i in range(t+1):
    E_T = E_T + P_i_next(i)/(mu_i(i) - lambda_i(i) )
  return E_T
  
def E_T_simplex_lb(t, arr_rate, gamma, mu, ub=False, naive=False):
  # def E_S_c(t, gamma, mu):
  #   return sum([binomial(t,i)*2**i * (-1)**(k-i) / (gamma+(2*k-i)*mu) for i in range(t+1) ] )
  # def E_S_c_2(t, gamma, mu):
  #   return sum([binomial(t,i)*2**i * (-1)**(k-i) * 2/(gamma+(2*k-i)*mu)**2 for i in range(t+1) ] )
  def E_S_p_m(m, t, gamma, mu):
    return sum([binomial(t-m,i)*2**i * (-1)**(t-m-i) / (gamma+(2*t-m-i)*mu) for i in range(t-m+1) ] )
  def E_S_p_m_2(m, t, gamma, mu):
    return sum([binomial(t-m,i)*2**i * (-1)**(t-m-i) * 2/(gamma+(2*t-m-i)*mu)**2 for i in range(t-m+1) ] )
    
  # E_S_c = E_S_c(t, gamma, mu)
  # E_S_c_2 = E_S_c_2(t, gamma, mu)
  E_S_p_m_l, E_S_p_m_2_l = [], []
  for m in range(t+1):
    E_S_p_m_l.append(E_S_p_m(m, t, gamma, mu) )
    E_S_p_m_2_l.append(E_S_p_m_2(m, t, gamma, mu) )
  # if t == 1:
  #   E_S_c = 2/(gamma+mu) - 1/(gamma+2*mu)
  #   E_S_c_2 = 4/(gamma+mu)**2 - 2/(gamma+2*mu)**2
  #   E_S_p = 1/(gamma+mu)
  #   E_S_p_2 = 2/(gamma+mu)**2
    
  #   f_c = 0.5
  #   E_S = f_c*E_S_c + (1-f_c)*E_S_p
  #   E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
  #   return E_S + arr_rate*E_S_2/2/(1-arr_rate*E_S)
  # else:
  #   log(ERROR, "not implemented!")
  #   return 1
  """
  nu = 2*t*mu+gamma
  P_01 = 0 # 2*t*mu/nu * mu/nu * (gamma+mu)/nu
  for i in range(t):
    prod = 1
    for j in range(1, i+1):
      prod *= 2*(t-j)*mu/nu
    P_01 += prod * (i+1)*mu/nu * (gamma+(i+1)*mu)/nu
  P_01 *= 2*t*mu/nu
  # print("_P_01= {}, P_01= {}".format(2*t*mu/nu * mu/nu * (gamma+mu)/nu, P_01) )
  
  P_10 = 0 # (gamma+mu)/nu + (1-(2*mu+gamma)/nu)*(2*mu+gamma)/nu
  for i in range(t):
    prod = 1
    for j in range(1, i+1):
      prod *= (1-(gamma+2*i*mu)/nu)
    P_10 += prod * (gamma+(i+1)*mu)/nu
  ro_max = P_01/P_10
  """
  E_X = 1/arr_rate
  E_S_min = sum(E_S_p_m_l)/len(E_S_p_m_l)
  # E_J_max = E_X/(E_X-E_S_min)
  # ro_max = 1 - 1/E_J_max
  # ro_max = compute_ro(t, arr_rate, gamma, mu)
  # ro_max = min(E_S_min/(E_X-E_S_min), 0.99)
  # ro_max = min(math.pow(E_S_min/(E_X-E_S_min), 1/t), 0.99)
  ro_max = E_S_min/E_X
  # print("ro_max= {}".format(ro_max) )
  
  p_0 = (1-ro_max)/(1-ro_max**(t+1) )
  if naive:
    p_0 = 1/(1+t*ro_max)
  def p_m(m):
    if naive:
      return 1/(t+1)
    else:
      return ro_max**m * p_0
    # return ro_max**(m != 0) * p_0
    # if m == t:
    #  m = m - 1
    # return 1/2**(m+1)
  # print("p_0= {}".format(p_0) )
  p_m_l = (t+1)*[0]
  for m in range(t+1):
    p_m_l[m] = p_m(m)
  # """
  # Trying to improve lower bound using incremental steps by adjusting ro_m
  if ub:
    ro_m_l = []
    for m in range(t+1):
      print("m= {}".format(m) )
      A = 0
      for i in range(m+1):
        A += numpy.prod(ro_m_l[:i] )
      E_Y = E_X - E_S_min
      B = (t-m)*numpy.prod(ro_m_l[:m] )
      print("A= {}, B= {}, (E_X-E_Y*A)/B/E_Y= {}".format(A, B, (E_X-E_Y*A)/B/E_Y) )
      ro_m = min((E_X-E_Y*A)/B/E_Y, 1)/2
      ro_m_l.append(ro_m)
      p_0 = 1/(sum([numpy.prod(ro_m_l[:i] ) for i in range(m+1) ] ) + numpy.prod(ro_m_l[:m+1] )*(t-m) )
      for m_ in range(t+1):
        p_m_l[m_] = numpy.prod(ro_m_l[:m_] )*p_0
      print("p_0= {}, ro_m_l= {}, p_m_l= {}".format(p_0, ro_m_l, p_m_l) )
  # """
  # E_S = sum(E_S_p_m_l)/len(E_S_p_m_l)
  # E_S_2 = sum(E_S_p_m_2_l)/len(E_S_p_m_2_l)
  E_S = sum([E_S_p_m_l[m]*p_m for m,p_m in enumerate(p_m_l) ] )
  E_S_2 = sum([E_S_p_m_2_l[m]*p_m for m,p_m in enumerate(p_m_l) ] )
  # E_S = 0.5*E_S_p_m_l[0] + 0.5*sum(E_S_p_m_l[1:] )/t
  # E_S_2 = 0.5*E_S_p_m_2_l[0] + 0.5*sum(E_S_p_m_2_l[1:] )/t
  
  return E_S + arr_rate*E_S_2/2/(1-arr_rate*E_S)

def compute_ro(t, arr_rate, gamma, mu):
  def E_S_p_m(m, t, gamma, mu):
    return sum([binomial(t-m,i)*2**i * (-1)**(t-m-i) / (gamma+(2*t-m-i)*mu) for i in range(t-m+1) ] )
  def E_S_p_m_2(m, t, gamma, mu):
    return sum([binomial(t-m,i)*2**i * (-1)**(t-m-i) * 2/(gamma+(2*t-m-i)*mu)**2 for i in range(t-m+1) ] )
  
  E_S_p_m_l, E_S_p_m_2_l = [], []
  for m in range(t+1):
    E_S_p_m_l.append(E_S_p_m(m, t, gamma, mu) )
    E_S_p_m_2_l.append(E_S_p_m_2(m, t, gamma, mu) )
  
  def p_m(m, ro):
    p_0 = (1-ro)/(1-ro**(t+1) )
    return ro**m * p_0
  def E_S(ro):
    return sum([E_S_p_m_l[m]*p_m(m,ro) for m in range(t+1) ] )
  
  E_X = 1/arr_rate
  E_S_min = E_S_p_m_l[0]
  ro = 0.99
  for i in range(10):
    ro = E_S(ro)/E_X
    print("i= {}, ro= {}".format(i, ro) )
  return ro

if __name__ == "__main__":
  # simplex_w_two_repair__E_T(0.9, 1.0)
  # simplex_w_one_repair__parametric_E_T()
  # simplex_steady_state_prob_hist()
  # E_T_simplex_lb(t=3, arr_rate=0.9, gamma=1, mu=1)
  # compute_ro(t=3, arr_rate=0.9, gamma=1, mu=1)
  
  # mu = 1.0
  # mc_truncation_state_id__state_prob_map_map = {}
  # for mc_truncation_state_id in range(1, 10, 1):
  #   mc_truncation_state_id__state_prob_map_map[mc_truncation_state_id] = simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu)
  # log(WARNING, "mc_truncation_state_id__state_prob_map_map= {}".format(pprint.pformat(mc_truncation_state_id__state_prob_map_map) ) )
  
  # simplex_w_one_repair__E_T_matrix_analytic(t=3, arr_rate=0.55, mu=1)
  plot_E_T_simplex_split_to_one(t=3, mu=1)
  