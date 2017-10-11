import math, random, numpy, csv
from scipy.stats import *
from patch import *

class RV(): # Random Variable
  def __init__(self, l_l, u_l):
    self.l_l = l_l
    self.u_l = u_l

class Exp(RV):
  def __init__(self, mu, D=0):
    RV.__init__(self, l_l=D, u_l=float("inf") )
    self.D = D
    self.mu = mu
  
  def __str__(self):
    return "Exp(D={}, mu={})".format(self.D, self.mu)
  
  def tail(self, x):
    if x <= self.l_l:
      return 1
    return math.exp(-self.mu*(x - self.D) )
  
  def cdf(self, x):
    if x <= self.l_l:
      return 0
    return 1 - math.exp(-self.mu*(x - self.D) )
  
  def pdf(self, x):
    if x <= self.l_l:
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

class TPareto(): # Truncated
  def __init__(self, l, u, a):
    RV.__init__(self, l_l=l, u_l=u)
    self.l = l
    self.u = u
    self.a = a
  
  def __str__(self):
    return "Pareto(l= {}, u= {}, a= {})".format(self.l, self.u, self.a)
  
  def cdf(self, x):
    if x < self.l: return 0
    elif x >= self.u: return 1
    else:
      return (1 - (self.l/x)**self.a)/(1 - (self.l/self.u)**self.a)
  
  def tail(self, x):
    return 1 - self.cdf(x)
  
  def mean(self):
    return self.moment(1)
  
  def moment(self, k):
    if k == self.a:
      return math.log(self.u_l/self.l)
    else:
      return self.a*self.l**k/(self.a-k) * \
             (1 - (self.l/self.u)**(self.a-k))/(1 - (self.l/self.u)**self.a)
  
  def gen_sample(self):
    u = random.uniform(0, 1)
    return self.l*(1 - u*(1-(self.l/self.u)**self.a) )**(-1/self.a)

class Google(RV):
  def __init__(self, k):
    RV.__init__(self, l_l=0, u_l=float("inf") )
    
    self.k = k
    self.sample_l = []
    # with open("filtered_task_lifetimes_for_jobs_w_num_task_{}.dat".format(k), mode="rt") as f:
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
  def __init__(self, L, U, p):
    RV.__init__(self, l_l=L, u_l=U)
    
    self.p = p
  
  def __str__(self):
    return "Bern(l= {}, u= {}, p= {})".format(self.l_l, self.u_l, self.p)
  
  def mean(self):
    return (1 - self.p)*self.l_l + self.p*self.u_l
  
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

class DUniform(): # Discrete
  def __init__(self, lb, ub):
    RV.__init__(self, l_l=lb, u_l=ub)
  
  def __str__(self):
    return "DUniform[{}, {}]".format(self.l_l, self.u_l)
  
  def mean(self):
    return (self.u_l + self.l_l)/2
  
  def pdf(self, x):
    return float(1/(self.u_l - self.l_l + 1) )
  
  def gen_sample(self):
    return random.randint(self.l_l, self.u_l)

class BoundedZipf():
  def __init__(self, lb, ub, a=1):
    RV.__init__(self, l_l=lb, u_l=ub)
    self.a = a
    
    self.v = numpy.arange(self.l_l, self.u_l+1) # values
    w_l = [float(v)**(-a) for v in self.v] # self.v**(-a) # weights
    self.p = [w/sum(w_l) for w in w_l]
    self.dist = scipy.stats.rv_discrete(name='bounded_zipf', values=(self.v, self.p) )
  
  def __str__(self):
    return "BoundedZipf([{}, {}], a= {})".format(self.l_l, self.u_l, self.a)
  
  def pdf(self, x):
    return self.dist.pmf(x)
  
  def cdf(self, x):
    # if x < self.l_l: return 0
    # elif x >= self.u_l: return 1
    # else:
    #   return sum(self.p[:(x-self.l_l+1) ] )
    return self.dist.cdf(x)
  
  def inv_cdf(self, p):
    return self.dist.ppf(p)
  
  def tail(self, x):
    return 1 - self.cfd(x)
  
  def mean(self):
    # return sum([v*self.p(i) for i,v in enumerate(self.v) ] )
    return self.dist.mean()
  
  def gen_sample(self):
    return self.dist.rvs(size=1)

class NegBinomial():
  def __init__(self, num_succ, p):
    RV.__init__(self, l_l=num_succ, u_l=float("Inf") )
    self.p = p
    
    self.dist = scipy.stats.nbinom(num_succ, p)
  
  def __str__(self):
    return "NegBinom[num_succ= {}, p= {}]".format(self.l_l, self.p)
  
  def cdf(self, x):
    return self.dist.cdf(x - self.l_l)
  
  def tail(self, x):
    return 1 - self.cdf(x)
  
  def gen_sample(self):
    return self.dist.rvs(size=1)

class Gamma():
  def __init__(self, num_exp, rate):
    RV.__init__(self, l_l=0, u_l=float("Inf") )
    
    self.shape, self.scale = num_exp, 1/rate
    # self.dist = numpy.random.gamma(shape, scale, size=1)
    self.dist = scipy.stats.gamma(self.shape, self.scale)
  
  def __str__(self):
    return "Gamma[shape= {}, scale= {}]".format(self.shape, self.scale)
  
  def cdf(self, x):
    return self.dist.cdf(x)
  
  def tail(self, x):
    return 1 - self.cdf(x)
  
  def gen_sample(self):
    return self.dist.rvs(size=1)
