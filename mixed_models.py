from patch import *

def E_T_mixed_net_ub(n, k, l, qlambda_l=[] ):
  if len(qlambda_l):
    E_T_l = []
    for i,l in enumerate(qlambda_l):
      qlambda_l_ = list(qlambda_l)
      qlambda_l_.remove(l)
      # print("l= {}, qlambda_l_= {}".format(l, qlambda_l_) )
      mu = sum(qlambda_l_[0:k] )
      E_T_l.append(1/(mu-l) )
    log(WARNING, "n= {}, k= {}, qlambda_l= {}\n\t E_T_l= {}".format(n, k, qlambda_l_, E_T_l) )
    return E_T_l
  else:
    E_S = 1/l * (H(n-1) - H(n-k) )
    E_S_2 = 1/l**2 * (H_2(n-1) - H_2(n-k) ) + E_S**2
    E_T = E_S + l*E_S_2/2/(1-l*E_S)
    log(WARNING, "n= {}, k= {}, l= {}\n\t E_T= {}".format(n, k, l, E_T) )
    if E_T < 0: return None
    return E_T

def E_T_mixed_net_lb(n, k, l):
  E_S = 1/(n-k+1)/l
  E_S_2 = 2/((n-k+1)*l)**2
  E_T = E_S + l*E_S_2/2/(1-l*E_S)
  log(WARNING, "n= {}, k= {}, l= {}\n\t E_T= {}".format(n, k, l, E_T) )
  if E_T < 0: return None
  return E_T