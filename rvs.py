import math, random, numpy, csv
from scipy.stats import *
from patch import *

class RV(): # Random Variable
  def __init__(self, l_l, u_l):
    self.l_l = l_l
    self.u_l = u_l

class Pareto(RV):
  def __init__(self, a, loc=1):
    RV.__init__(self, l_l=loc, u_l=float("inf") )
    self.a = a
    self.loc = loc
  
  def __str__(self):
    return "Pareto(a= {}, loc= {})".format(self.a, self.loc)
  
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
      return self.a*self.loc/(self.a-1)
  
  def var(self):
    if self.a <= 2:
      log(WARNING, "Variance is Infinity; a= {} <= 2".format(self.a) )
      return float("inf")
    else:
      return self.a*self.loc**2 / (self.a-1)**2/(self.a-2)
  
  def gen_sample(self):
    return ((numpy.random.pareto(self.a, 1) + 1)*self.loc)[0]
    # return pareto.ppf(numpy.random.uniform(0, 1), b=self.a, scale=self.loc)

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

class Google(RV):
  def __init__(self, k):
    RV.__init__(self, l_l=0, u_l=float("inf") )
    
    self.k = k
    self.sample_l = []
    with open("task_lifetimes_for_jobs_w_num_task_{}.dat".format(k), mode="rt") as f:
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
