import math

from rvs import *
from patch import *

def f(x, dist_m):
  dist = dist_m['dist']
  if dist == 'Exp':
    mu = dist_m['mu']
    rv = Exp(mu)
  elif dist == 'SExp':
    D, mu = dist_m['D'], dist_m['mu']
    rv = Exp(mu, D)
  elif dist == 'Pareto':
    loc, a = dist_m['l'], dist_m['a']
    rv = Pareto(loc, a)
  return rv.pdf(x)

def F(x, dist_m):
  dist = dist_m['dist']
  if dist == 'Exp':
    mu = dist_m['mu']
    rv = Exp(mu)
  elif dist == 'SExp':
    D, mu = dist_m['D'], dist_m['mu']
    rv = Exp(mu, D)
  elif dist == 'Pareto':
    loc, a = dist_m['l'], dist_m['a']
    rv = Pareto(loc, a)
  return rv.cdf(x)

def Pr_X_n_k_leq_x(n, k, x, dist_m):
  p = F(x, dist_m)
  return sum([binomial(n,i)*p**i*(1-p)**(n-i) for i in range(k, n+1) ] )

def EXm_n_k(m, n, k, dist_m):
  return mpmath.quad(lambda t: m*t**(m-1)*(1 - Pr_X_n_k_leq_x(n, k, t, dist_m) ), [0, mpmath.inf] )

def EXm(m, dist_m):
  return mpmath.quad(lambda t: m*t**(m-1)*(1 - F(t, dist_m) ), [0, mpmath.inf] )

def scale_dist(dist_m, k):
  dist = dist_m['dist']
  dist_m_ = {'dist': dist}
  if dist == 'SExp':
    dist_m_['D'] = dist_m['D']/k
    dist_m_['mu'] = dist_m['mu']
  elif dist == 'Pareto':
    dist_m_['l'] = dist_m['l']/k
    dist_m_['a'] = dist_m['a']
  return dist_m_