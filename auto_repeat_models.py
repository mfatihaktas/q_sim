import sys, pprint, random, math, numpy, simpy, getopt, itertools, sympy

from patch import *

# Pr{T >= t}
def prob_T_k_n_geq_t(mu, d, k, n, t):
  q = 1 - math.exp(-mu*d)
  
  # Pr{T >= t | T < d}*Pr{T < d}
  def lhs():
    if t > d:
      return 0
    q_ = 1 - math.exp(-mu*t)
    # return q**k - q_**k
    # return 1 - q_**k / q**k
    return q**k - q_**k
  # Pr{T >= t | T > d}*Pr{T > d}
  def rhs():
    # if t <= d:
    #   return 0
    def prob_X_n_r__k_r_leq_tau(r):
      tau = max(0, t - d)
      q_ = 1 - math.exp(-mu*tau)
      # print("q_= {}".format(q_) )
      sum_ = 0
      for j in range(k-r, n-r+1):
        sum_ += binomial(n-r,j) * q_**j * (1-q_)**(n-r-j)
      return sum_
    sum_ = 0
    for r in range(k):
      # print("prob_X_n_r__k_r_leq_tau(r= {})= {}".format(r, prob_X_n_r__k_r_leq_tau(r) ) )
      sum_ += prob_X_n_r__k_r_leq_tau(r) * binomial(k,r) * q**r * (1-q)**(k-r)
    return (1 - q**k - sum_) # / q**k
  # print("lhs= {}, rhs= {}".format(lhs(), rhs() ) )
  return lhs() + rhs()

def prob_T_k_n_geq_t_alt(mu, d, k, n, t):
  q = 1 - math.exp(-mu*d)
  
  # Pr{T >= t | T < d}*Pr{T < d}
  def lhs():
    if t > d:
      return 0
    q_ = 1 - math.exp(-mu*t)
    # return q**k - q_**k
    # return 1 - q_**k / q**k
    return q**k - q_**k
  # Pr{T >= t | T > d}*Pr{T > d}
  def rhs():
    # if t <= d:
    #   return 0
    def prob_X_n_r__k_r_leq_tau(r):
      tau = max(0, t - d)
      q_ = 1 - math.exp(-mu*tau)
      
      return sympy.mpmath.quad(lambda x: x**(k-r-1) * (1-x)**(n-k), [0, q_] ) / \
             sympy.mpmath.quad(lambda x: x**(k-r-1) * (1-x)**(n-k), [0, 1] )
    
    # sum_ = 0
    # for r in range(k):
    #   # print("prob_X_n_r__k_r_leq_tau(r= {})= {}".format(r, prob_X_n_r__k_r_leq_tau(r) ) )
    #   sum_ += prob_X_n_r__k_r_leq_tau(r) * binomial(k,r) * q**r * (1-q)**(k-r)
    
    # sum_ = prob_X_n_r__k_r_leq_tau((math.ceil(k*q) + math.floor(k*q) )/2) - prob_X_n_r__k_r_leq_tau(k)*q**k # math.ceil(k*q)
    sum_ = prob_X_n_r__k_r_leq_tau(k*q) - prob_X_n_r__k_r_leq_tau(k)*q**k # math.ceil(k*q)
    return (1 - q**k - sum_) # / q**k
  # print("lhs= {}, rhs= {}".format(lhs(), rhs() ) )
  return lhs() + rhs()

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
  q = 1 - math.exp(-mu*d)
  
  E_H_n_r = 0
  for r in range(k):
    E_H_n_r += ((H(n-r) )/mu) * binomial(k, r)*(1-p)**r * p**(k-r)
  
  print("d= {}, k= {}, n= {}\n\t E_H_n_r - H(n-k)= {}".format(d, k, n, E_H_n_r - H(n-k) ) )
  
  return d - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
         1/mu*(E_H_n_r - H(n-k) )

# Is equal to E_T_k_n
def E_T_k_n_alt(mu, d, k, n):
  p = math.exp(-mu*d)
  q = 1-math.exp(-mu*d)
  E_H_n_r = 0
  for r in range(k+1):
    E_H_n_r += H(n-r) * binomial(k, r) * q**r * (1-q)**(k-r)
  
  # E_H_n_r = H(n-math.ceil(k*q) ) # H(n-math.floor(k*q) ) # H(n-math.ceil(k*q) ) # H(n-k) # (H(n) + H(n-k) )/2 # + 1/n*(1-q)**k
  # E_H_n_r = math.log(n-k*q) + 0.5772156649 + 1/2/(n-k*q)
  
  # return d - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
  # return d - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
  # return d + -H(k) + \
  # return d + 1/mu*(H(k) - mu*d - sympy.mpmath.quad(lambda x: (1 - (1-x)**k)/x, [0, p] ) ) + \ # great
  
  # return d - 1/mu*sympy.mpmath.quad(lambda x: x**k * 1/(1-x), [0, q] ) + \
  #       (E_H_n_r/mu - H(n-k)/mu)
  if n == k:
    return d - 1/mu*sympy.mpmath.quad(lambda x: x**k * 1/(1-x), [0, q] )
  
  return d - 1/mu*sympy.mpmath.quad(lambda x: x**k * 1/(1-x), [0, q] ) + \
         1/mu*(E_H_n_r - H(n-k) )
        # 1/mu*(math.log((n-k*q)/(n-k) ) )
        # 1/mu*(math.log((n-k*q)/(n-k) ) + 1/2/(n-k*q) - 1/2/(n-k) )

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