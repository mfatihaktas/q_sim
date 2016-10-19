import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt
from math import factorial
from numpy import linalg
from patch import *

def simplex_w_two_repair__state_prob_map(mu):
  P = []
  state_list = ["0_0", "0_1", "0_2", "1_0", "1_1", "1_2", "2_0", "2_1", "2_2"]
  for i in range(len(state_list) ):
    P.append(len(state_list)*[0] )
  state_index_map = {s:i for i,s in enumerate(state_list) }
  
  gamma = mu
  nu = 4*mu + gamma
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
  
  log(WARNING, "P= \n{}".format(numpy.matrix(P) ) )
  # for p in numpy.arange(10, 100, 10):
  #   P_p = linalg.matrix_power(P, p)
  #   log(WARNING, "p= {}, P_p=\n{}".format(p, pprint.pformat(P_p) ) )
  pi_v = linalg.matrix_power(P, 100)[0]
  # log(WARNING, "pi_v= {}".format(pprint.pformat(pi_v) ) )
  return {state_list[i]:pi for i,pi in enumerate(pi_v) }

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

def simplex_w_two_repair__serv_time(arr_rate, mu):
  E_S_s_c = avq_low_traff_serv_time_first_moment(2, 1, mu)
  E_S_s_c_2 = avq_low_traff_serv_time_second_moment(2, 1, mu)
  E_S_p_1st_kind = 5/(12*mu)
  E_S_p_1st_kind_2 = 23/(72*mu**2)
  E_S_p_2nd_kind = 1/(3*mu)
  E_S_p_2nd_kind_2 = 2/(9*mu**2)
  
  state_prob_map = simplex_w_two_repair__state_prob_map(mu)
  log(WARNING, "state_prob_map=\n {}".format(pprint.pformat(state_prob_map) ) )
  
  gamma = mu
  nu = 4*mu + gamma
  mc_truncation_state_id = int(math.sqrt(len(state_prob_map) ) ) - 1
  state_sum = sum([state_prob_map["0_"+str(i)] + state_prob_map[str(i)+"_0"] for i in range(1, mc_truncation_state_id+1, 1) ] )
  f_jd = state_prob_map["0_0"]*gamma/nu \
       + (mu+gamma)/nu * state_sum \
       + (2*mu+gamma)/nu * (1 - state_prob_map["0_0"] - state_sum)
  f_jc = state_prob_map["0_0"]*gamma/nu \
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
  
if __name__ == "__main__":
  simplex_w_two_repair__serv_time(1.0)
  