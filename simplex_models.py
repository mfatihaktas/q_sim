import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, math, numpy
import itertools

from math import factorial
from numpy import linalg
from rvs import *
from patch import *

def simplex_sym_l__sym__rgroup_l_m(t):
  sym_l = []
  sym__rgroup_l_m = {}
  
  if t == 1: sym_l = ['a', 'b']
  elif t == 3: sym_l = ['a', 'b', 'c']
  elif t == 7: sym_l = ['a', 'b', 'c', 'd']
  for sym in sym_l:
    rgroup_l = []
    if t == 1:
      if sym == 'a':
        rgroup_l.append([0] )
        rgroup_l.append([1, 2] )
      elif sym == 'b':
        rgroup_l.append([1] )
        rgroup_l.append([0, 2] )
    elif t == 3:
      if sym == 'a':
        rgroup_l.append([0] )
        rgroup_l.append([1, 2] )
        rgroup_l.append([3, 4] )
        rgroup_l.append([5, 6] )
      elif sym == 'b':
        rgroup_l.append([1] )
        rgroup_l.append([0, 2] )
        rgroup_l.append([3, 5] )
        rgroup_l.append([4, 6] )
      elif sym == 'c':
        rgroup_l.append([3] )
        rgroup_l.append([0, 4] )
        rgroup_l.append([1, 5] )
        rgroup_l.append([2, 6] )
    elif t == 7:
      if sym == 'a':
        rgroup_l.append([0] )
        rgroup_l.append([1, 2] )
        rgroup_l.append([3, 4] )
        rgroup_l.append([5, 6] )
        rgroup_l.append([7, 8] )
        rgroup_l.append([9, 10] )
        rgroup_l.append([11, 12] )
        rgroup_l.append([13, 14] )
      elif sym == 'b':
        rgroup_l.append([1] )
        rgroup_l.append([0, 2] )
        rgroup_l.append([3, 5] )
        rgroup_l.append([4, 6] )
        rgroup_l.append([7, 9] )
        rgroup_l.append([8, 10] )
        rgroup_l.append([11, 13] )
        rgroup_l.append([12, 14] )
      elif sym == 'c':
        rgroup_l.append([3] )
        rgroup_l.append([0, 4] )
        rgroup_l.append([1, 5] )
        rgroup_l.append([2, 6] )
        rgroup_l.append([7, 11] )
        rgroup_l.append([8, 12] )
        rgroup_l.append([9, 13] )
        rgroup_l.append([10, 14] )
      elif sym == 'd':
        rgroup_l.append([7] )
        rgroup_l.append([0, 8] )
        rgroup_l.append([1, 9] )
        rgroup_l.append([2, 10] )
        rgroup_l.append([3, 11] )
        rgroup_l.append([4, 12] )
        rgroup_l.append([5, 13] )
        rgroup_l.append([6, 14] )
    sym__rgroup_l_m[sym] = rgroup_l
  return (sym_l, sym__rgroup_l_m)

def tompecs_sym_l__sym__rgroup_l_m(scheme):
  sym_l = ['a', 'b', 'c', 'd', 'e', 'f']
  sym__rgroup_l_m = {}
  
  if scheme == 'lrc':
    for sym in sym_l:
      rgroup_l = []
      if sym == 'a':
        rgroup_l.append([0] )
        rgroup_l.append([1, 2, 6] )
        for c in itertools.combinations([1, 2, 3, 4, 5, 8, 9], 6):
          rgroup_l.append(c)
      elif sym == 'b':
        rgroup_l.append([1] )
        rgroup_l.append([0, 2, 6] )
        for c in itertools.combinations([0, 2, 3, 4, 5, 8, 9], 6):
          rgroup_l.append(c)
      elif sym == 'c':
        rgroup_l.append([2] )
        rgroup_l.append([0, 1, 6] )
        for c in itertools.combinations([0, 1, 3, 4, 5, 8, 9], 6):
          rgroup_l.append(c)
      elif sym == 'd':
        rgroup_l.append([3] )
        rgroup_l.append([4, 5, 7] )
        for c in itertools.combinations([0, 1, 2, 4, 5, 8, 9], 6):
          rgroup_l.append(c)
      elif sym == 'e':
        rgroup_l.append([4] )
        rgroup_l.append([3, 5, 7] )
        for c in itertools.combinations([0, 1, 2, 3, 5, 8, 9], 6):
          rgroup_l.append(c)
      elif sym == 'f':
        rgroup_l.append([5] )
        rgroup_l.append([3, 4, 7] )
        for c in itertools.combinations([0, 1, 2, 3, 4, 8, 9], 6):
          rgroup_l.append(c)
      sym__rgroup_l_m[sym] = rgroup_l
  elif scheme == 'mds':
    for sym in sym_l:
      rgroup_l = []
      if sym == 'a':
        rgroup_l.append([0] )
        for c in itertools.combinations([1, 2, 3, 4, 5, 6, 7, 8], 6):
          rgroup_l.append(c)
      elif sym == 'b':
        rgroup_l.append([1] )
        for c in itertools.combinations([0, 2, 3, 4, 5, 6, 7, 8], 6):
          rgroup_l.append(c)
      elif sym == 'c':
        rgroup_l.append([2] )
        for c in itertools.combinations([0, 1, 3, 4, 5, 6, 7, 8], 6):
          rgroup_l.append(c)
      elif sym == 'd':
        rgroup_l.append([3] )
        for c in itertools.combinations([0, 1, 2, 4, 5, 6, 7, 8], 6):
          rgroup_l.append(c)
      elif sym == 'e':
        rgroup_l.append([4] )
        for c in itertools.combinations([0, 1, 2, 3, 5, 6, 7, 8], 6):
          rgroup_l.append(c)
      elif sym == 'f':
        rgroup_l.append([5] )
        for c in itertools.combinations([0, 1, 2, 3, 4, 6, 7, 8], 6):
          rgroup_l.append(c)
      sym__rgroup_l_m[sym] = rgroup_l
  elif scheme == 'rep':
    rgroup_l = []
    if sym == 'a':
      rgroup_l.append([0] )
      rgroup_l.append([1] )
      rgroup_l.append([2] )
    elif sym == 'b':
      rgroup_l.append([3] )
      rgroup_l.append([4] )
      rgroup_l.append([5] )
    elif sym == 'c':
      rgroup_l.append([6] )
      rgroup_l.append([7] )
      rgroup_l.append([8] )
    elif sym == 'd':
      rgroup_l.append([9] )
      rgroup_l.append([10] )
      rgroup_l.append([11] )
    elif sym == 'e':
      rgroup_l.append([12] )
      rgroup_l.append([13] )
      rgroup_l.append([14] )
    elif sym == 'f':
      rgroup_l.append([15] )
      rgroup_l.append([16] )
      rgroup_l.append([17] )
    sym__rgroup_l_m[sym] = rgroup_l
  else:
    log(ERROR, "Unrecognized scheme= {}".format(scheme) )
    return
  
  return sym_l, sym__rgroup_l_m

MAX_ET = 100 # 35 # 70

def E_T_rep_n_1(ar, mu, n):
  E_S = 1/n/mu
  E_S_2 = 2/(n*mu)**2/mu**2
  
  E_T = E_S + ar*E_S_2/2/(1-ar*E_S)
  if E_T < 0 or E_T > MAX_ET: return None
  # if E_T < 0: return None
  return E_T

def E_T_rep_n_1_split_to_one(ar, mu, n):
  E_T = 1/(mu - ar/n)
  if E_T < 0 or E_T > MAX_ET: return None
  return E_T

def E_T_fj_2(ar, mu):
  ro = ar/mu
  E_T = (12 - ro)/(8*(mu - ar) )
  if E_T < 0: return None
  return E_T

# ****************************  Simplex possible serv time moments  **************************** #
def ES(r, t, mu):
  if t == 0:
    return 1/mu
  else:
    return 1/(mu*r)*B(t+1, 1/r)

def ES2(r, t, mu):
  if t == 0:
    return 2/mu**2
  else:
    s = 0
    for j in range(t + 1):
      _s = 0
      for l in range(r*j + 1):
        _s += (-1)**l * binom(r*j, l)*(2/(mu**2 * (l+1)**2) )
      s += binom(t, j) * (-1)**j * _s
    return s

def ESm_typei(m, t, i, sdist_m):
  V = rv_from_m(sdist_m)
  V_2_2 = X_n_k(V, 2, 2)
  def Pr_Sgs(s):
    return V.tail(s)**(i+1) * V_2_2.tail(s)**(t-i)
  
  return float(mpmath.quad(lambda s: m*s**(m-1) * Pr_Sgs(s), [0, 10000*10] ) )

def ES_typei(t, i, sdist_m):
  sdist = sdist_m['dist']
  if sdist == "Exp":
    mu = sdist_m['mu']
    gamma = mu
    return sum([binom(t-i,j)*2**j * (-1)**(t-i-j) / (gamma+(2*t-i-j)*mu) for j in range(t-i+1) ] )
  elif sdist == "Pareto":
    loc, a = sdist_m['loc'], sdist_m['a']
    if a*(t+1) < 1:
      return None
    return sum([binom(t-i,j)*2**j * (-1)**(t-i-j) * loc*(1 + 1/(a*(2*t+1-i-j)-1) ) for j in range(t-i+1) ] )
  elif sdist == "Bern":
    U, L, p = sdist_m['U'], sdist_m['L'], sdist_m['p']
    ps = p**(i+1) * (1 - (1-p)**2)**(t-i)
    return U + (L-U)*ps
  elif sdist == "BernPareto":
    U, L, p, loc, a = sdist_m['U'], sdist_m['L'], sdist_m['p'], sdist_m['loc'], sdist_m['a']
    if a <= 1:
      return None
    ps = p**(i+1) * (1 - (1-p)**2)**(t-i)
    return (U + (L-U)*ps) * loc*(1 + 1/(a-1) )

def ES2_typei(t, i, sdist_m):
  sdist = sdist_m['dist']
  if sdist == "Exp":
    mu = sdist_m['mu']
    gamma = mu
    return sum([binom(t-i,j)*2**j * (-1)**(t-i-j) * 2/(gamma+(2*t-i-j)*mu)**2 for j in range(t-i+1) ] )
  elif sdist == "Pareto":
    loc, a = sdist_m['loc'], sdist_m['a']
    if a*(t+1)/2 < 1:
      return None
    return sum([binom(t-i,j)*2**j * (-1)**(t-i-j) * loc**2*(1 + 1/(a*(2*t+1-i-j)/2-1) ) for j in range(t-i+1) ] )
  elif sdist == "Bern":
    U, L, p = sdist_m['U'], sdist_m['L'], sdist_m['p']
    ps = p**(i+1) * (1 - (1-p)**2)**(t-i)
    return U**2 + (L**2 - U**2)*ps
  elif sdist == "BernPareto":
    U, L, p, loc, a = sdist_m['U'], sdist_m['L'], sdist_m['p'], sdist_m['loc'], sdist_m['a']
    if a <= 2:
      return None
    ps = p**(i+1) * (1 - (1-p)**2)**(t-i)
    return (U**2 + (L**2 - U**2)*ps) * loc**2*(1 + 1/(a/2-1) )
  else:
    return None

def reptoall_innerbound_on_ar(t, sdist_m):
  # if serv == "Exp":
  #   if w_sys:
  #     # E_S = 1/(mu*r)*B(t+1, 1/r)
  #     E_S = ES_typei(mu, t, i=0) #   else:
  #     E_S = 1/mu * sum([binom(t,i) * 2**i*(-1)**(t-i)/(2*t-i) for i in range(t+1) ] )
  E_S = ES_typei(t, 0, sdist_m)
  return float(1/E_S)

# -----------------------------------  Simplex w/ split-to-one  ------------------------------- #
def arub_simplex_selectone(t, serv, sdist_m):
  if serv == "Exp":
    rv = Exp(sdist_m['mu'] )
  elif serv == "Pareto":
    rv = Pareto(sdist_m['loc'], sdist_m['a'] )
  elif serv == "Bern*Pareto":
    rv = BernPareto(sdist_m['U'], sdist_m['L'], sdist_m['p'], sdist_m['loc'], sdist_m['a'] )
  return (t+1)*1/rv.mean()

def ET_selectone(t, ar, mu, p_r=None):
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
    ar_ = ar*split_prob_l[g]
    if g == 0:
      E_T_mm1 = 1/(gamma-ar_)
      if E_T_mm1 < 0: return None
      E_T += split_prob_l[g] * E_T_mm1
    else:
      E_T_fj = E_T_fj_2(ar_, mu)
      if E_T_fj == None:
        return None
      E_T += split_prob_l[g] * E_T_fj
  if E_T < 0 or E_T > 20: return None
  return E_T

def plot_ET_selectone(t, mu):
  ar_ub = ar_ub_simplex_split_to_one(t, mu)
  
  # for ar in numpy.linspace(0.05, ar_ub, 5):
  for ar in numpy.linspace(0.05, ar_ub, 5):
    p_r_l, E_T_l = [], []
    for p_r in numpy.linspace(0.05, 1/t, 20):
      p_r_l.append(p_r)
      E_T_l.append(ET_selectone(t, ar, mu, p_r=p_r) )
    l_str = "{0:.2f}".format(ar)
    plot.plot(p_r_l, E_T_l, label=r'$t= {}, \mu= {}, \lambda= {}$'.format(t, mu, l_str), color=next(dark_color), marker=next(marker) )
    plot.legend()
  plot.xlabel(r'$p_r$')
  plot.ylabel(r'$E[T]$')
  plot.savefig("plot_ET_selectone.png")

def ET_selectone_min(t, ar, mu):
  # print("ar= {}".format(ar) )
  ET_pre = None
  for p_r in numpy.linspace(0.05, 1/t, 100):
    ET = ET_selectone(t, ar, mu, p_r)
    if ET_pre == None:
      ET_pre = ET
    else:
      if ET == None or ET >= ET_pre:
        return ET_pre
      else:
        ET_pre = ET

# -----------------------------------------  Simplex(t=1)  ----------------------------------- #
def ET_reptoall_t1(ar, mu, c=None):
  if c == None:
    E_S_c = 0.665/mu
    E_S_c_2 = 7/(9*(mu**2) )
    E_S_p = 0.5/mu
    E_S_p_2 = 0.5/(mu**2)
    E_S = 0.6*E_S_c + 0.4*E_S_p
    E_S_2 = 0.6*E_S_c_2 + 0.4*E_S_p_2
    E_T = E_S + (ar/2)*E_S_2/(1 - ar*E_S) # w/ M/G/1 assumption
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
    return E_S + ar*E_S_2/2/(1-ar*E_S)

def ET_reptoall_t1_matrixanalytic(t, ar, mu):
  l = ar
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
  # print("PI_0_v= {}, PI_1_v= {}".format(PI_0_v, PI_1_v) )
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
  # return sum([n*prob for n,prob in N_prob_map.items() ] )/ar
  
  E_N = PI_0_v*numpy.ones((4,1)) - PI_0_v[0, 0] + \
        PI_1_v*((numpy.identity(5)-R)**-2 + (numpy.identity(5)-R)**-1)*numpy.ones((5,1))
  # print("E_N= {}".format(E_N) )
  E_T = E_N[0, 0]/ar
  if E_T > 50: return None
  return E_T

# -----------------------------------------  Simplex(t=2)  ----------------------------------- #
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

def E_T_simplex_t_2(ar, mu, M):
  E_S_c = ES(2, 2, mu)
  E_S_c_2 = ES2(2, 2, mu)
  E_S_p_1st = 5/(12*mu)
  E_S_p_1st_2 = 23/(72*mu**2)
  E_S_p_2nd = 1/(3*mu)
  E_S_p_2nd_2 = 2/(9*mu**2)
  
  mc_truncation_state_id = M # 5
  state_prob_map = simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu)
  # log(WARNING, "state_prob_map=\n {}".format(pprint.pformat(state_prob_map) ) )
  
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
  E_T = E_S + (ar/2)*E_S_2/(1 - ar*E_S)
  return E_T

def E_T_simplex_t_2_wo_sys(ar, mu):
  return 1/(4*mu-ar) + 1/(3*mu-ar) + 2/3/(2*mu-ar)

def E_T_simplex_t_1_parametric():
  Cap = 3
  def parametric_E_T(ar, c):
    f_c = 1/(1 + 2/c/(c+2) )
    # mu = Cap/(c+2)
    E_S_c = 2*(c+2)/Cap/(c+1) - 1/Cap
    E_S_c_2 = (2*(c+2)/Cap/(c+1) )**2 - 2/(Cap**2)
    E_S_p = (c+2)/Cap/(c+1)
    E_S_p_2 = 2*((c+2)/Cap/(c+1) )**2
    
    E_S = f_c*E_S_c + (1-f_c)*E_S_p
    E_S_2 = f_c*E_S_c_2 + (1-f_c)*E_S_p_2
    return E_S + ar*E_S_2/2/(1-ar*E_S)
  color = iter(cm.rainbow(numpy.linspace(0, 1, 20) ) )
  for ar in numpy.arange(0.05, 1.0, 0.1):
    c_l = numpy.arange(0.05, 5, 0.05)
    E_T_l = [parametric_E_T(ar, c) for c in c_l]
    plot.plot(c_l, E_T_l, '.', color=next(color), label=r'$\lambda$={0:0.2f}'.format(ar) )
    plot.legend()
  plot.xlabel("c")
  plot.ylabel("E[T]")
  plot.title(r'Total service rate= {}'.format(Cap) )
  plot.savefig("E_T_simplex_t_1_parametric.png")

# ##############################  UB; Simplex(t=2:[sys, MDS(r,2)])  ############################# #
def E_S_c_avq_sys__mds_r_2(gamma, mu, r):
  return r/(gamma+(r-1)*mu) - (r-1)/(gamma+r*mu)

def E_S_c_2_avq_sys__mds_r_2(gamma, mu, r):
  return 2*r/(gamma+(r-1)*mu)**2 - 2*(r-1)/(gamma+r*mu)**2

def E_T_avq_sys__mds_r_2(ar, gamma, mu, r):
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
  E_T = E_S + (ar/2)*E_S_2/(1 - ar*E_S)
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
  E_T = E_S + (ar/2)*E_S_2/(1 - ar*E_S)
  return E_T

# ########################################  Simplex(t)  ###################################### #
def ET_simplex_lb(t, ar, sdist_m):
  E_S = ES_typei(t, t, sdist_m)
  E_S_2 = ES2_typei(t, t, sdist_m)
  E_T = E_S + (ar/2)*E_S_2/(1 - ar*E_S)
  if E_T < 0 or E_T > MAX_ET: return None
  return E_T

def ET_simplex_sm(t, ar, sdist_m):
  E_S = ES_typei(t, 0, sdist_m)
  E_S_2 = ES2_typei(t, 0, sdist_m)
  E_T = E_S + (ar/2)*E_S_2/(1 - ar*E_S)
  if E_T < 0 or E_T > MAX_ET: return None
  return E_T

def ET_simplex_splitmerge_wo_sys(t, ar, mu, c=None):
  E_S = 1/mu * sum([binom(t,i) * 2**i*(-1)**(t-i)/(2*t-i) for i in range(t+1) ] )
  E_S_2 = 2/mu**2 * sum([binom(t,i) * 2**i*(-1)**(t-i)/(2*t-i) for i in range(t+1) ] )
  return E_S + (ar/2)*E_S_2/(1 - ar*E_S)

def plot_diminishing_return_as_t_incs():
  serv = "Exp"
  sdist_m = {'mu': 1}
  mu = 1
  m_l = []
  ES_simplex_typei_l, ES2_simplex_typei_l = [], []
  E_S_rep_t_l, E_S_2_rep_t_l = [], []
  for t in range(1, 20):
    m_l.append(t)
    ES_simplex_typei_l.append(ES_typei(t, t, serv, sdist_m) )
    ES2_simplex_typei_l.append(ES2_typei(t, t, serv, sdist_m) )
    E_S_rep_t_l.append(1/(t+1)/mu)
    E_S_2_rep_t_l.append(2/((t+1)*mu)**2 )
  plot.plot(m_l, ES_simplex_typei_l, color=next(dark_color), label=r'$E[S], simplex$', marker=next(marker), linestyle=':', mew=2)
  plot.plot(m_l, ES2_simplex_typei_l, color=next(dark_color), label=r'$E[S^2]$, simplex', marker=next(marker), linestyle=':', mew=2)
  plot.plot(m_l, E_S_rep_t_l, color=next(dark_color), label=r'$E[S]$, rep', marker=next(marker), linestyle=':', mew=2)
  plot.plot(m_l, E_S_2_rep_t_l, color=next(dark_color), label=r'$E[S^2]$, rep', marker=next(marker), linestyle=':', mew=2)
  
  plot.legend()
  plot.xlabel(r'$t$')
  plot.ylabel("")
  plot.title(r'type-$t$ (Fastest) service start, $mu= {}$'.format(mu) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_diminishing_return_as_t_incs.pdf")
  plot.gcf().clear()
  log(WARNING, "done.")

def ET_simplex_approx(t, ar, sdist_m, p_i_l=[], naive=False, incremental=False, ar_ub=False):
  ES_simplex_typei_l, ES2_simplex_typei_l = [], []
  for i in range(t+1):
    ES_simplex_typei_l.append(ES_typei(t, i, sdist_m) )
    ES2_simplex_typei_l.append(ES2_typei(t, i, sdist_m) )
  
  if len(p_i_l) == 0:
    E_X = 1/ar
    E_S_min = sum(ES_simplex_typei_l)/len(ES_simplex_typei_l)
    ro_max = E_S_min/E_X
    
    p_0 = (1-ro_max)/(1-ro_max**(t+1) )
    # if naive:
    #   p_0 = 1/(1+t*ro_max)
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
        print(">> i= {}".format(i) )
        A = sum([numpy.prod(ro_i_l[:j] ) for j in range(i+1) ] )
        
        E_Y = E_X - E_S_min
        B = (t-i)*numpy.prod(ro_i_l[:i] )
        print("A= {}, B= {}, (E_X-E_Y*A)/B/E_Y= {}".format(A, B, (E_X-E_Y*A)/B/E_Y) )
        ro_i = min((E_X-E_Y*A)/B/E_Y, 1)/2 # !!! /2 is necessary while comparing Straightforward vs Better vs Fine-grained
        ro_i_l.append(ro_i)
        p_0 = 1/(sum([numpy.prod(ro_i_l[:j] ) for j in range(i+1) ] ) + numpy.prod(ro_i_l[:i+1] )*(t-i) )
        for i_ in range(t+1):
          p_i_l[i_] = numpy.prod(ro_i_l[:i_] )*p_0
        print("p_0= {}, ro_i_l= {}, p_i_l= {}".format(p_0, ro_i_l, p_i_l) )
  E_S = sum([ES_simplex_typei_l[i]*p_i for i,p_i in enumerate(p_i_l) ] )
  if ar_ub:
    return 1/E_S
  E_S_2 = sum([ES2_simplex_typei_l[i]*p_i for i,p_i in enumerate(p_i_l) ] )
  E_T = E_S + ar*E_S_2/2/(1-ar*E_S)
  if E_T < 0 or E_T > MAX_ET: return None, None
  return E_T, p_i_l

def E_T_simplex_varki_gauri_lb(t, ar, gamma, mu):
  def P_i_next(i):
    return 2*(t-i)*mu/(gamma+(2*t-i)*mu)
  def mu_i(i):
    return gamma + (2*t-i)*mu
  def lambda_i(i):
    if i == 0:
      return ar
    else:
      prod = 1
      for k in range(i):
        prod = prod * P_i_next(k)
      return ar*prod
  E_T = 0
  # for i in range(2):
  for i in range(t+1):
    E_T = E_T + P_i_next(i)/(mu_i(i) - lambda_i(i) )
  return E_T

if __name__ == "__main__":
  # mu = 1.0
  # mc_truncation_state_id__state_prob_map_iap = {}
  # for mc_truncation_state_id in range(1, 10, 1):
  #   mc_truncation_state_id__state_prob_map_iap[mc_truncation_state_id] = simplex_w_two_repair__state_prob_map(mc_truncation_state_id, mu)
  # log(WARNING, "mc_truncation_state_id__state_prob_map_iap= {}".format(pprint.pformat(mc_truncation_state_id__state_prob_map_iap) ) )
  
  # ET_reptoall_t1_matrixanalytic(t=3, ar=0.55, mu=1)
  # plot_ET_selectone(t=1, mu=1)
  
  sdist_m = {'dist': 'Exp', 'mu': 1}
  t = 3
  for i in range(t+1):
    ESi = ES_typei(t, i, sdist_m)
    ESi_ = ESm_typei(1, t, i, sdist_m)
    print("ESi= {}, ESi_= {}".format(ESi, ESi_) )
