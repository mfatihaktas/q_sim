import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy, sympy, csv
from math import factorial
from numpy import linalg
from patch import *

def simplex_t_1_mc():
  l, g, a, b, z = sympy.var('l g a b z')
  jq_limit = 2
  num_job_limit = 10
  s_i_map = {}
  state_list = []
  def state(n_j, n_a, n_b):
    return str(n_j) + ",(" + str(n_a) + "," + str(n_b) + ")"
  for n_j in range(num_job_limit + 1):
    for n_a in range(jq_limit + 1):
      for n_b in range(jq_limit + 1):
        if not (n_a == 0 or n_b == 0) or n_b > n_j or n_a > n_j:
          continue
        s = state(n_j, n_a, n_b)
        state_list.append(s)
        s_i_map[s] = len(state_list) - 1
  print("asd= {}".format("asd") )
  log(WARNING, "s_i_map= {}".format(pprint.pformat(s_i_map) ) )
  Q = sympy.zeros(len(state_list) )
  # Q = Matrix(len(state_list), len(state_list), z)
  log(WARNING, "Q= {}".format(pprint.pformat(Q) ) )
  for n_j in range(num_job_limit + 1):
    for n_a in range(jq_limit + 1):
      for n_b in range(jq_limit + 1):
        if not (n_a == 0 or n_b == 0) or n_b > n_j or n_a > n_j:
          continue
        s = state(n_j, n_a, n_b)
        if n_j == 0:
          Q[s_i_map[s], s_i_map[s]] = -l
          Q[s_i_map[s], s_i_map["1,(0,0)"]] = g
          Q[s_i_map[s], s_i_map["1,(0,1)"]] = a + g
          Q[s_i_map[s], s_i_map["1,(1,0)"]] = b + g
        elif n_j != 0 and (n_a == 0 and n_b == 0):
          Q[s_i_map[s], s_i_map[s]] = -(l + g + a + b)
          s_u = state(n_j-1, n_a, n_b)
          Q[s_i_map[s], s_i_map[s_u]] = l
          if n_j != num_job_limit:
            s_d = state(n_j+1, n_a, n_b)
            Q[s_i_map[s], s_i_map[s_d]] = g
            s_dl = state(n_j+1, n_a, n_b+1)
            Q[s_i_map[s], s_i_map[s_dl]] = a + g
            s_dr = state(n_j+1, n_a+1, n_b)
            Q[s_i_map[s], s_i_map[s_dr]] = b + g
        elif n_a == 0:
          if n_b == n_j: # all the way on the left
            Q[s_i_map[s], s_i_map[s]] = -(l + g + a)
          else:
            Q[s_i_map[s], s_i_map[s]] = -(l + g + a + b)
            s_u = state(n_j-1, n_a, n_b)
            Q[s_i_map[s], s_i_map[s_u]] = l
          s_r = state(n_j, n_a, n_b-1)
          Q[s_i_map[s], s_i_map[s_r]] = b
          if n_b != jq_limit and n_j != num_job_limit:
            s_dl = state(n_j+1, n_a, n_b+1)
            Q[s_i_map[s], s_i_map[s_dl]] = a + g
        elif n_b == 0:
          if n_a == n_j: # all the way on the right
            Q[s_i_map[s], s_i_map[s]] = -(l + g + b)
          else:
            Q[s_i_map[s], s_i_map[s]] = -(l + g + a + b)
            s_u = state(n_j-1, n_a, n_b)
            Q[s_i_map[s], s_i_map[s_u]] = l
          s_l = state(n_j, n_a-1, n_b)
          Q[s_i_map[s], s_i_map[s_l]] = a
          if n_a != jq_limit and n_j != num_job_limit:
            s_dr = state(n_j+1, n_a+1, n_b)
            Q[s_i_map[s], s_i_map[s_dr]] = b + g
  
  Q = Q.subs(g, 0)
  m = sympy.var('m')
  Q = Q.subs(a, m)
  Q = Q.subs(b, m)
  log(WARNING, "Q= {}".format(pprint.pformat(Q) ) )
  file = open('Q.csv', 'w')
  
  writer = csv.writer(file)
  writer.writerow(["", *state_list] )
  for i in range(len(state_list) ):
    writer.writerow([state_list[i], *Q.row(i) ] )
  file.close()

def simplex_t_q__steady_state_dist():
  """
  # mu, gamma= 1
  arr_rate__zero_state_prob_map= {
    0.10000000000000001: 0.932628,
    0.20000000000000001: 0.867572,
    0.30000000000000004: 0.802282,
    0.40000000000000002: 0.741966,
    0.5: 0.673994,
    0.59999999999999998: 0.614784,
    0.70000000000000007: 0.55389,
    0.80000000000000004: 0.490902,
    0.90000000000000002: 0.426464,
    1.0: 0.36783,
    1.1000000000000001: 0.302002}
  # mu= 1, gamma= 0.5
  arr_rate__zero_state_prob_map= {
  0.10000000000000001: 0.904448,
  0.20000000000000001: 0.817524,
  0.30000000000000004: 0.724832,
  0.40000000000000002: 0.6379,
  0.5: 0.555546,
  0.59999999999999998: 0.479598,
  0.70000000000000007: 0.397546,
  0.80000000000000004: 0.324622,
  0.90000000000000002: 0.25593,
  1.0: 0.18136,
  1.1000000000000001: 0.123022,
  1.2000000000000002: 0.062016}
  """
  mu, gamma = 1, 1 # 0.5
  lambda_ = 0.7 # 0.9 # 0.9 # 1.1 # 0.5
  log(WARNING, "lambda_= {}, mu= {}, gamma= {}".format(lambda_, mu, gamma) )
  ro = lambda_/mu
  # tau = lambda_/(2*mu+gamma)
  # tau = (1.02*lambda_)/(2*mu+gamma) # for lambda_= 0.1
  # tau = (1.05*lambda_)/(2*mu+gamma) # for lambda_= 0.2
  # tau = (1.07*lambda_)/(2*mu+gamma) # for lambda_= 0.3
  # tau = (1.09*lambda_)/(2*mu+gamma) # for lambda_= 0.4
  # tau = (1.13*lambda_)/(2*mu+gamma) # for lambda_= 0.5
  # tau = (1.15*lambda_)/(2*mu+gamma) # for lambda_= 0.6
  # tau = (1.18*lambda_)/(2*mu+gamma) # for lambda_= 0.7
  # tau = (1.22*lambda_)/(2*mu+gamma) # for lambda_= 0.8
  # tau = (1.26*lambda_)/(2*mu+gamma) # for lambda_= 0.9
  # tau = (1.3*lambda_)/(2*mu+gamma) # for lambda_= 1
  # tau = (1.36*lambda_)/(2*mu+gamma) # for lambda_= 1.1
  _tau = lambda_/(2*mu+gamma)
  # _tau = lambda_/(0.9*(2*mu+gamma) + 0.1*(mu+gamma) )
  # _tau = _tau*(1+_tau/(2*mu+gamma) )
  # _tau = _tau*(1+_tau)
  # tau =  _tau*(1+lambda_*_tau)
  # x = 0.98 # for lambda_=0.1
  # x = 0.95 # for lambda_=0.2
  # x = 0.93 # for lambda_=0.3
  # x = 0.90 # for lambda_=0.4
  # x = 0.86 # for lambda_=0.5
  # x = 0.84 # for lambda_=0.6
  # x = 0.81 # for lambda_=0.7
  # tau = lambda_/(x*(2*mu+gamma) + (1-x)*lambda_) # _tau
  tau = _tau*(1+0.3*lambda_)
  # d = _tau/2 # 0.15
  # tau = lambda_/(2*mu+gamma)
  s = lambda_+2*mu+gamma
  # tau_ = _tau/(1-_tau*(mu+gamma)*mu/s*(1/s+1/(s-mu) ) )
  k = mu*(mu+gamma)/s*(1/s+1/(s-mu))
  tau_ = _tau/(1-_tau*k/(1-k*lambda_/(lambda_+mu+gamma) ) )
  __tau = tau
  # tau = tau_
  # tau_0 = _tau/(1-_tau*2*mu*(mu+gamma)/s**2)
  k_ = 2*mu*(mu+gamma)/s**2
  tau_0 = tau # lambda_/(lambda_+2*mu+gamma) # _tau # _tau/(1-_tau*k_/(1-k_*lambda_/(lambda_+2*mu+gamma) ) )
  # final_tau = (-gamma*tau_0**2 + (lambda_+2*mu+2*gamma)*tau_0)/lambda_ - 1
  # tau = final_tau
  # print("final_tau= {}".format(final_tau) )
  print("__tau= {}, tau= {}, tau_0= {}".format(__tau, tau, tau_0) )
  
  r = 1/(1+gamma/mu)
  """
  def dist(k, i):
    pi_0 = 1/(1/(1-tau) + 2*tau*r/(1-ro)/(1-ro*r) )
    if i == 0:
      return pi_0*tau**k
    else:
      return pi_0*tau*ro**(k-1) * r**i
  def short_dist(k, i):
    A = (1-tau)*(1-tau*r)/(1+tau*r)
    return tau**k * A * r**i
  # B = (1-tau)*(1-tau*r)/(1+tau*r)
  # print("B= {}".format(B) )
  """
  s = lambda_ + mu + gamma
  b = s/(mu+gamma)
  a = -tau*mu/(mu+gamma)
  delta = b**2 + 4*a
  print("delta= {}".format(delta) )
  r_0 = (-b - math.sqrt(delta) )/2/a
  r_1 = (-b + math.sqrt(delta) )/2/a
  print("r_0= {}, r_1= {}".format(r_0, r_1) )
  # Y = (r_0*pi_0 + (pi_1 - pi_0*b)*r_0*r_1)/(r_0-r_1)
  # X = pi_0 - Y
  pi_1_c = (lambda_-gamma*tau_0)/2/(mu+gamma)
  # pi_0 = (1-tau)/(1 + 2*(r_0+(pi_1_c-b)*r_0*r_1 + r_1-1)/(r_0-1)/(r_1-1) )
  pi_0 = 1/(1/(1-tau_0) + 2/(1-tau)*(r_0+(pi_1_c-b)*r_0*r_1+r_1-1)/(r_0-1)/(r_1-1) )
  print("pi_0= {}".format(pi_0) )
  Y = (r_0+(pi_1_c-b)*r_0*r_1)*pi_0/(r_0-r_1)
  X = pi_0 - Y
  def final_dist(k, i):
    return tau**(k-i) * (X/(r_0**i) + Y/(r_1**i) )
  for k in range(0, 3+1):
    for i in range(0, k+1):
      # print("short_dist:: p= {}".format(short_dist(k, i) ) )
      print("k= {}, i= {}; p= {}".format(k, i, final_dist(k, i) ) )

def plot_simplex_start_type_moments(t, gamma, mu):
  def binomial(x, y):
    try:
      binom = factorial(x) // factorial(y) // factorial(x - y)
    except ValueError:
      binom = 0
    return binom
  
  def E_S_p_m(m, t, gamma, mu):
    return sum([binomial(t-m,i)*2**i * (-1)**(t-m-i) / (gamma+(2*t-m-i)*mu) for i in range(t-m+1) ] )
  def E_S_p_m_2(m, t, gamma, mu):
    return sum([binomial(t-m,i)*2**i * (-1)**(t-m-i) * 2/(gamma+(2*t-m-i)*mu)**2 for i in range(t-m+1) ] )
  
  E_S_p_m_l, E_S_p_m_2_l = [], []
  m_l = []
  for m in range(t+1):
    m_l.append(m)
    E_S_p_m_l.append(E_S_p_m(m, t, gamma, mu) )
    E_S_p_m_2_l.append(E_S_p_m_2(m, t, gamma, mu) )
  
  plot.plot(m_l, E_S_p_m_l, 'o-', label="E_S_p_m")
  plot.plot(m_l, E_S_p_m_2_l, 'o-', label="E_S_p_m_2")
  plot.legend()
  plot.xlabel("m")
  plot.ylabel("Time")
  plot.savefig("plot_simplex_start_type_moments.png")

def job_start_fractions(t):
  P = numpy.zeros(shape=(t+1, t+1) )
  for i in range(t+1):
    to_right_sum = 1/(2+0.2)
    if i == t:
      to_right_sum = 0
    j_g_i_l = [to_right_sum/(t-i) for j in range(i+1, t+1) ]
    j_leq_i_l = [(1-to_right_sum)/(i+1) for j in range(0, i+1) ]
    # j_g_i_l = [to_right_sum/(t-i+1) for j in range(i, t+1) ]
    # j_leq_i_l = [(1-to_right_sum)/i for j in range(0, i) ]
    P[i, :] = j_leq_i_l + j_g_i_l
  print("P= \n{}".format(P) )
  P_100 = linalg.matrix_power(P, 100)
  print("P_100= \n{}".format(P_100) )
  pi_v = P_100[0]
  # log(WARNING, "pi_v= {}".format(pprint.pformat(pi_v) ) )
  state_prob_map = {i:pi for i,pi in enumerate(pi_v) }
  print("state_prob_map= \n{}".format(pprint.pformat(state_prob_map) ) )

if __name__ == "__main__":
  # simplex_t_1_mc()
  # simplex_t_q__steady_state_dist()
  # job_start_fractions(5)
  plot_simplex_start_type_moments(10, 1, 1)
  
