import math, numpy, pprint

from patch import *
from simplex_models import avq_low_traff_serv_time_first_moment, avq_low_traff_serv_time_second_moment

def E_T_pop_lb_ff_simplex(arr_rate, t, mu):
  r = 2
  E_S = avq_low_traff_serv_time_first_moment(r, t, mu)
  E_S_2 = avq_low_traff_serv_time_second_moment(r, t, mu)
  
  E_T = E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  if E_T > 100: return None
  return E_T

def E_T_pop_ub_ff_simplex(arr_rate, t, mu):
  r = 2
  n_server = t*r + 1
  n_sym = int(math.log(n_server+1, 2) )
  
  t_min = (n_server - 1 - (n_sym-1)*2)/r
  
  E_S = avq_low_traff_serv_time_first_moment(r, t_min, mu)
  E_S_2 = avq_low_traff_serv_time_second_moment(r, t_min, mu)
  
  E_T = E_S + (arr_rate/2)*E_S_2/(1 - arr_rate*E_S)
  if E_T > 100: return None
  return E_T
