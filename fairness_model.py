import math, numpy, pprint

from patch import *
from simplex_models import avq_low_traff_serv_time_first_moment, avq_low_traff_serv_time_second_moment

r = 2
E_T_UB = 50

def pop_ar_ub_max_ff_simplex(t, mu):
  E_S = avq_low_traff_serv_time_first_moment(r, t, mu)
  return float(1/E_S) - 0.01

def E_T_pop_lb_ff_simplex(pop_ar, t, mu):
  E_S = avq_low_traff_serv_time_first_moment(r, t, mu)
  E_S_2 = avq_low_traff_serv_time_second_moment(r, t, mu)
  
  E_T = E_S + (pop_ar/2)*E_S_2/(1 - pop_ar*E_S)
  if E_T > E_T_UB: return None
  return E_T

def pop_ar_ub_min_ff_simplex(t, mu):
  n_server = t*r + 1
  n_sym = int(math.log(n_server+1, 2) )
  t_min = int((n_server - 1 - (n_sym-1)*2)/r)
  E_S = avq_low_traff_serv_time_first_moment(r, t_min, mu)
  return float(1/E_S) - 0.01

def E_T_pop_ub_ff_simplex(pop_ar, t, mu):
  n_server = t*r + 1
  n_sym = int(math.log(n_server+1, 2) )
  t_min = int((n_server - 1 - (n_sym-1)*2)/r)
  
  E_S = avq_low_traff_serv_time_first_moment(r, t_min, mu)
  E_S_2 = avq_low_traff_serv_time_second_moment(r, t_min, mu)
  
  E_T = E_S + (pop_ar/2)*E_S_2/(1 - pop_ar*E_S)
  if E_T < 0 or E_T > E_T_UB: return None
  return E_T

def pop_ar_ub_approx_ff_simplex(unpop_ar, t, mu):
  if t == 1:
    p_busyb = unpop_ar/mu
    E_S = p_busyb*1/mu + (1-p_busyb)*avq_low_traff_serv_time_first_moment(r, 1, mu)
  elif t == 3:
    p_busyb = p_busyc = unpop_ar/mu
    E_S = p_busyb*p_busyc*avq_low_traff_serv_time_first_moment(r, 1, mu) + \
          ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*avq_low_traff_serv_time_first_moment(r, 2, mu) + \
          (1-p_busyb)*(1-p_busyc)*avq_low_traff_serv_time_first_moment(r, 3, mu)
  return float(1/E_S) - 0.01

def E_T_pop_approx_ff_simplex(pop_ar, unpop_ar, t, mu):
  if t == 1:
    p_busyb = unpop_ar/mu
    
    E_S = p_busyb*1/mu + (1-p_busyb)*avq_low_traff_serv_time_first_moment(r, 1, mu)
    E_S_2 = p_busyb*2/mu**2 + (1-p_busyb)*avq_low_traff_serv_time_second_moment(r, 1, mu)
  elif t == 3:
    p_busyb = p_busyc = unpop_ar/mu
    
    E_S = p_busyb*p_busyc*avq_low_traff_serv_time_first_moment(r, 1, mu) + \
          ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*avq_low_traff_serv_time_first_moment(r, 2, mu) + \
          (1-p_busyb)*(1-p_busyc)*avq_low_traff_serv_time_first_moment(r, 3, mu)
    E_S_2 = p_busyb*p_busyc*avq_low_traff_serv_time_second_moment(r, 1, mu) + \
            ((1-p_busyb)*p_busyc + p_busyb*(1-p_busyc) )*avq_low_traff_serv_time_second_moment(r, 2, mu) + \
            (1-p_busyb)*(1-p_busyc)*avq_low_traff_serv_time_second_moment(r, 3, mu)
  else:
    log(ERROR, "Unexpected t= {}".format(t) )
  
  E_T = E_S + (pop_ar/2)*E_S_2/(1 - pop_ar*E_S)
  if E_T > E_T_UB: return None
  return E_T
