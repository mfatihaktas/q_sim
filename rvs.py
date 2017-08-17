import math, random, numpy, csv
from scipy.stats import *
from patch import *

class RV(): # Random Variable
  def __init__(self, l_l, u_l):
    self.l_l = l_l
    self.u_l = u_l

class Exp(RV):
  def __init__(self, mu, D=0):
    RV.__init__(self, l_l=0, u_l=float("inf") )
    self.D = D
    self.mu = mu
    if D < 0 or mu < 0:
      log(ERROR, "Unexpected arg; {}".format(self) )
      return 1
  
  def __str__(self):
    return "Exp(D= {0:.2f}, mu= {0:.2f})".format(self.D, self.mu)
  
  def tail(self, x):
    if x <= self.l_l or x < self.D:
      return 1
    return math.exp(-self.mu*(x - self.D) )
  
  def cdf(self, x):
    if x <= self.l_l or x < self.D:
      return 0
    return 1 - math.exp(-self.mu*(x - self.D) )
  
  def pdf(self, x):
    if x <= self.l_l or x < self.D:
      return 0
    return self.mu*math.exp(-self.mu*(x - self.D) )
  
  def mean(self):
    return self.D + 1/self.mu
  
  def var(self):
    return 1/self.mu**2
  
  def gen_sample(self):
    return self.D + random.expovariate(self.mu)

class Pareto(RV):
  def __init__(self, loc, a):
    RV.__init__(self, l_l=loc, u_l=float("inf") )
    self.loc = loc
    self.a = a
  
  def __str__(self):
    return "Pareto(loc= {}, a= {})".format(self.loc, self.a)
  
  def tail(self, x):
    if x <= self.l_l:
      return 1
    return (self.loc/x)**self.a
  
  def cdf(self, x):
    if x <= self.l_l:
      return 0
    return 1 - (self.loc/x)**self.a
  
  def pdf(self, x):
    if x <= self.l_l:
      return 0
    return self.a*self.loc**self.a / x**(self.a+1)
  
  def dpdf_dx(self, x):
    if x <= self.l_l:
      return 0
    return sympy.mpmath.diff(lambda y: self.a*self.loc**self.a / y**(self.a+1), x)
  
  def mean(self):
    if self.a <= 1:
      log(WARNING, "Mean is Infinity; a= {} <= 1".format(self.a) )
      return float("inf")
    else:
      return self.loc*self.a/(self.a-1)
  
  def var(self):
    if self.a <= 2:
      log(WARNING, "Variance is Infinity; a= {} <= 2".format(self.a) )
      return float("inf")
    else:
      return self.a*self.loc**2 / (self.a-1)**2/(self.a-2)
  
  def gen_sample(self):
    return ((numpy.random.pareto(self.a, 1) + 1)*self.loc)[0]
    # return pareto.ppf(numpy.random.uniform(0, 1), b=self.a, scale=self.loc)

class Google(RV):
  def __init__(self, k):
    RV.__init__(self, l_l=0, u_l=float("inf") )
    
    self.k = k
    self.sample_l = []
    with open("filtered_task_lifetimes_for_jobs_w_num_task_{}.dat".format(k), mode="rt") as f:
      reader = csv.reader(f)
      for line in reader:
        self.sample_l.append(float(line[0] ) )
    self.sample_l.sort()
    self.num_sample = len(self.sample_l)
  
  def __str__(self):
    return "Google(k= ".format(self.k)
  
  def mean(self):
    return sum(self.sample_l)/self.num_sample
  
  def gen_sample(self):
    return self.sample_l[math.floor(self.num_sample*random.random() ) ]

class Dolly(RV):
  # Kristen et al. A Better Model for Job Redundancy: Decoupling Server Slowdown and Job Size
  def __init__(self):
    RV.__init__(self, l_l=1, u_l=12)
  
  def __str__(self):
    return "Dolly"
  
  def gen_sample(self):
    u = random.uniform(0, 1)
    if u <= 0.23: return 1 + u/100
    u -= 0.23
    if u <= 0.14: return 2 + u/100
    u -= 0.14
    if u <= 0.09: return 3 + u/100
    u -= 0.09
    if u <= 0.03: return 4 + u/100
    u -= 0.03
    if u <= 0.08: return 5 + u/100
    u -= 0.08
    if u <= 0.1: return 6 + u/100
    u -= 0.1
    if u <= 0.04: return 7 + u/100
    u -= 0.04
    if u <= 0.14: return 8 + u/100
    u -= 0.14
    if u <= 0.12: return 9 + u/100
    u -= 0.12
    if u <= 0.021: return 10 + u/100
    u -= 0.021
    if u <= 0.007: return 11 + u/100
    u -= 0.007
    if u <= 0.002: return 12 + u/100
    return 12 + u/100 # for safety

class Bern(RV):
  def __init__(self, U, L, p):
    RV.__init__(self, l_l=U, u_l=L)
    
    self.p = p
  
  def __str__(self):
    return "Bern"
  
  def mean(self):
    return self.p*self.u_l + (1 - self.p)*self.l_l
  
  def gen_sample(self):
    u = random.uniform(0, 1)
    return self.u_l + u/100 if u <= self.p else self.l_l + u/100

# class BernPareto(RV):
#   def __init__(self, U, L, p, loc, a):
#     RV.__init__(self, l_l=U*loc, u_l=float("Inf") )
    
#     self.bern = Bern(U, L, p)
#     self.pareto = Pareto(loc, a)
  
#   def __str__(self):
#     return "Bern*Pareto"
  
#   def mean(self):
#     return self.bern.mean()*self.pareto.mean()
  
#   def gen_sample(self):
#     return self.bern.gen_sample()*self.pareto.gen_sample()
