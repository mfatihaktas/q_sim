import sys, pprint, random, math, numpy, simpy, getopt, itertools, sympy

from patch import *

def try_approx():
  mu = 1
  d = 2
  k = 20
  actual = sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] )
  
  p = math.exp(-mu*d)
  # approx = -1/mu*(H(k) - (k*p + 1 - (1-p)**k)/2 - mu*d)
  approx = -1/mu*(H(k) - mu*d)
  print("actual= {}, approx= {}, err= {}".format(actual, approx, abs(actual-approx) ) )

# Send k initially, send n-k more after d
def E_T_k_n(mu, d, k, n):
  p = math.exp(-mu*d)
  q = 1-math.exp(-mu*d)
  
  E_T = d*q**k - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] )
  sum_ = 0 # d*(1 - q**k)
  for r in range(k):
     sum_ += (d + (H(n-r)-H(n-k) )/mu) * binomial(k, r)*(1-p)**r * p**(k-r)
  E_T += sum_
  
  return E_T

# Is equal to E_T_k_n
def E_T_k_n_alt(mu, d, k, n):
  p = math.exp(-mu*d)
  q = 1-math.exp(-mu*d)
  E_H_n_r = 0
  for r in range(k+1):
     E_H_n_r += H(n-r) * binomial(k, r) * q**r * (1-q)**(k-r)
  
  # return H(k)/mu * q**k + (d*(1 - q**k) + E_H_n_r/mu - H(n-k)/mu)*(1 - q**k)
  # return d*q**k + \
  # return d*q**k + -1/mu*(H(k) - (k*p + 1 - (1-p)**k)/2 - mu*d) + \
  # return d*q**k + d - 1/mu*H(k) + \
  # return d*q**k + q*H(k) + \
  
  # E_H_n_r = H(n-math.ceil(k*q) ) # H(n-math.floor(k*q) ) # H(n-math.ceil(k*q) ) # H(n-k) # (H(n) + H(n-k) )/2 # + 1/n*(1-q)**k
  
  # return d - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
  # return d - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
  # return d + -H(k) + \
  # return d + 1/mu*(H(k) - mu*d - sympy.mpmath.quad(lambda x: (1 - (1-x)**k)/x, [0, p] ) ) + \ # great
  return d + sympy.mpmath.quad(lambda x: x**k * 1/(1-x), [0, q] ) + \
         (E_H_n_r/mu - H(n-k)/mu)

# E[H(k-r)] where r ~ Binomial(k, p)
def E_H_k_r(k, p):
  # p = random.random()
  E = 0
  for r in range(k+1):
    E += H(k-r) * binomial(k, r) * p**r * (1-p)**(k-r)
  return E

def E_T_k_n_lb(mu, d, k, n):
  p = math.exp(-mu*d)
  q = 1-math.exp(-mu*d)
  E_T = d*q**k - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] )
  
  # sum_ += math.log(n-k*p) - k*p*(1-p)/(n-k*p)**2 \
  #       - (math.log(k-k*p) - k*p*(1-p)/(k-k*p)**2) \
  #       + 0.5*(1/(n-k*p) + k*p*(1-p)/(n-k*p)**3 \
  #               - (1/(k-k*p) + k*p*(1-p)/(k-k*p)**3) ) \
  #       - H(n-k)*(1-p)**k
  E_H_n_r = H(n) # math.log(n-k*p) - k*p*(1-p)/(n-k*p)**2 + 0.5772156649
  # E_H_n_r = math.log(2*n+1-2*k*p) - 2*k*p*(1-p)/(2*n+1-2*k*p)**2 - math.log(2) + 0.5772156649
  sum_ = d*(1 - q**k) + E_H_n_r/mu - H(n-k)/mu
  E_T += sum_
  
  return E_T

if __name__ == "__main__":
  try_approx()