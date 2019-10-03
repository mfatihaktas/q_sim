import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

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
    # return "Exp(D={}, mu={})".format(self.D, self.mu)
    return r'Exp(\mu={})'.format(self.mu)
  
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
  
  def sample(self):
    return self.D + random.expovariate(self.mu)

class Pareto(RV):
  def __init__(self, loc, a):
    RV.__init__(self, l_l=loc, u_l=float("inf") )
    self.loc = loc
    self.a = a
  
  def __str__(self):
    # return "Pareto(loc= {}, a= {})".format(self.loc, self.a)
    return r'Pareto(s= {}, \alpha= {})'.format(self.loc, self.a)
  
  def to_latex(self):
    return r"${}(\min= {}, \alpha= {})$".format(r'\mathrm{Pareto}', round(self.loc, 2), round(self.a, 2) )
  
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
  
  def sample(self):
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
  
  def to_latex(self):
    return r"${}(\min= {}, \max= {}, \alpha= {})$".format(r'\mathrm{TPareto}', round(self.l, 2), round(self.u, 2), round(self.a, 2) )
    
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
  
  def sample(self):
    u = random.uniform(0, 1)
    return self.l*(1 - u*(1-(self.l/self.u)**self.a) )**(-1/self.a)

def plot_gensample_check():
  l, u, a = 1, 10**5, 2
  rv = TPareto(l, u, a)
  
  x_l = []
  for i in range(10**5):
    x_l.append(rv.sample() )
  x_l = numpy.sort(x_l)
  x_l = x_l[::-1]
  # i_ = None
  # for i in range(len(x_l)-1, 0, -1):
  #   if x_l[i] > 1.01: i_ = i; break
  # x_l = x_l[:i_]
  y_l = numpy.arange(x_l.size)/x_l.size
  plot.plot(x_l, y_l, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  y_l = []
  for x in x_l:
    y_l.append(rv.tail(x) )
  plot.plot(x_l, y_l, label=r'$Pareto(l= %.2f, u= %.2f, \alpha= %.2f)$' % (l, u, a), color=next(dark_color), linestyle='-')
  plot.legend()
  plot.xscale('log')
  plot.yscale('log')
  plot.xlabel(r'$x$', fontsize=13)
  plot.ylabel(r'$p(X > x)$', fontsize=13)
  plot.title(r'$X \sim$ {}'.format(rv) )
  plot.savefig("plot_gensample_check.png")
  plot.gcf().clear()

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
  
  def sample(self):
    return self.sample_l[math.floor(self.num_sample*random.random() ) ]

class SimRV(RV):
  def __init__(self, sample_l):
    RV.__init__(self, l_l=min(sample_l), u_l=max(sample_l) )
    
    self.sample_l = sample_l
    self.num_sample = len(self.sample_l)
  
  def __str__(self):
    return "SimRV"
  
  def mean(self):
    return sum(self.sample_l)/self.num_sample
  
  def sample(self):
    return self.sample_l[math.floor(self.num_sample*random.random() ) ]

class ExplicitRV(RV):
  def __init__(self, v_l, p_l):
    RV.__init__(self, l_l=min(v_l), u_l=max(v_l) )
    
    self.v_l = v_l
    self.p_l = p_l
    self.dist = scipy.stats.rv_discrete(name='dolly', values=(v_l, p_l) )
  
  def __repr__(self):
    return "ExplicitRV[\n\tv_l= {}, \n\tp_l= {}]".format(self.v_l, self.p_l)
  
  def sample(self):
    return self.dist.rvs()

class Dolly(RV):
  ## Kristen et al. A Better Model for Job Redundancy: Decoupling Server Slowdown and Job Size
  def __init__(self):
    RV.__init__(self, l_l=1, u_l=12)
    
    self.v = numpy.arange(1, 13)
    self.p = [0.23, 0.14, 0.09, 0.03, 0.08, 0.1, 0.04, 0.14, 0.12, 0.021, 0.007, 0.002]
    self.dist = scipy.stats.rv_discrete(name='dolly', values=(self.v, self.p) )
  
  def __str__(self):
    return "Dolly[{}, {}]".format(self.l_l, self.u_l)
  
  def pdf(self, x):
    return self.dist.pmf(x) if (x >= self.l_l and x <= self.u_l) else 0
  
  def cdf(self, x):
    if x < self.l_l:
      return 0
    elif x > self.u_l:
      return 1
    return float(self.dist.cdf(x) )
  
  def sample(self):
    u = random.uniform(0, 1)
    # if u <= 0.23: return 1 + u/100
    # u -= 0.23
    # if u <= 0.14: return 2 + u/100
    # u -= 0.14
    # if u <= 0.09: return 3 + u/100
    # u -= 0.09
    # if u <= 0.03: return 4 + u/100
    # u -= 0.03
    # if u <= 0.08: return 5 + u/100
    # u -= 0.08
    # if u <= 0.1: return 6 + u/100
    # u -= 0.1
    # if u <= 0.04: return 7 + u/100
    # u -= 0.04
    # if u <= 0.14: return 8 + u/100
    # u -= 0.14
    # if u <= 0.12: return 9 + u/100
    # u -= 0.12
    # if u <= 0.021: return 10 + u/100
    # u -= 0.021
    # if u <= 0.007: return 11 + u/100
    # u -= 0.007
    # if u <= 0.002: return 12 + u/100
    # return 12 + u/100 # for safety
    return self.dist.rvs() + u/100

class Bern(RV):
  def __init__(self, L, U, p):
    RV.__init__(self, l_l=L, u_l=U)
    
    self.p = p
  
  def __str__(self):
    return "Bern(l= {}, u= {}, p= {})".format(self.l_l, self.u_l, self.p)
  
  def mean(self):
    return (1 - self.p)*self.l_l + self.p*self.u_l
  
  def sample(self):
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
  
#   def sample(self):
#     return self.bern.sample()*self.pareto.sample()

class DUniform(): # Discrete
  def __init__(self, lb, ub):
    RV.__init__(self, l_l=lb, u_l=ub)
  
  def __str__(self):
    return "DUniform[{}, {}]".format(self.l_l, self.u_l)
  
  def mean(self):
    return (self.u_l + self.l_l)/2
  
  def pdf(self, x):
    return float(1/(self.u_l - self.l_l + 1) )
  
  def sample(self):
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
  
  def sample(self):
    return self.dist.rvs(size=1)

class Binomial():
  def __init__(self, n, p):
    RV.__init__(self, l_l=0, u_l=n)
    self.n = n
    self.p = p
    
    self.dist = scipy.stats.binom(n, p)
  
  def __str__(self):
    return "Binom[n= {}, p= {}]".format(self.n, self.p)
  
  def pdf(self, x):
    return self.dist.pmf(x)
    
  def cdf(self, x):
    return self.dist.cdf(x)
  
  def tail(self, x):
    return 1 - self.cdf(x)
  
  def sample(self):
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
  
  def sample(self):
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
  
  def sample(self):
    return self.dist.rvs(size=1)

class X_n_k():
  def __init__(self, X, n, k):
    RV.__init__(self, l_l=X.l_l, u_l=X.u_l)
    self.X, self.n, self.k = X, n, k
  
  def __str__(self):
    return "{}_{{}:{}}".format(self.X, self.n, self.k)
  
  def pdf(self, x):
    return self.n*self.X.pdf(x) * binom(self.n-1, self.k-1) * self.X.cdf(x)**(self.k-1) * self.X.tail(x)**(self.n-self.k)
  
  def cdf(self, x):
    # return cdf_n_k(self.X, self.n, self.k, x)
    cdf = 0
    for i in range(self.k, self.n+1):
      cdf += binom(self.n, i) * self.X.cdf(x)**i * self.X.tail(x)**(self.n-i)
    return cdf
  
  def tail(self, x):
    return 1 - self.cdf(x)
  
  def sample(self):
    return gen_orderstat_sample(self.X, self.n, self.k)

def moment_ith(i, X):
  return float(mpmath.quad(lambda x: i*x**(i-1) * (1 - X.cdf(x) ), [0, mpmath.inf] ) ) # 10000*10

def rv_from_m(dist_m):
  d = dist_m['dist']
  if d == 'Exp':
    return Exp(dist_m['mu'] )
  elif d == "SExp":
    return Exp(dist_m['mu'], dist_m['D'] )
  elif d == "Pareto":
    return Pareto(dist_m['l'], dist_m['a'] )
  elif d == "TPareto":
    return TPareto(dist_m['l'], dist_m['u'], dist_m['a'] )
  elif d == "Bern":
    return Bern(dist_m['U'], dist_m['L'], dist_m['p'] )
  elif d == 'Dolly':
    return Dolly()

def binom(n, k):
  # if n == k:
  #   return 1
  # elif k == 1:
  #   return n
  # elif k == 0:
  #   return 1
  # elif k > n:
  #   return 0
  # else:
  #   return math.factorial(n)/math.factorial(k)/math.factorial(n-k)
  return scipy.special.binom(n, k)

if __name__ == "__main__":
  plot_gensample_check()
