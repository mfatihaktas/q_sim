import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, random, math, numpy, getopt, itertools, mpmath, textwrap

from rvs import *
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

def prob_T_k_n_geq_t_approx(mu, d, k, n, t):
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
      
      return mpmath.quad(lambda x: x**(k-r-1) * (1-x)**(n-k), [0, q_] ) / \
             mpmath.quad(lambda x: x**(k-r-1) * (1-x)**(n-k), [0, 1] )
    
    # sum_ = 0
    # for r in range(k):
    #   # print("prob_X_n_r__k_r_leq_tau(r= {})= {}".format(r, prob_X_n_r__k_r_leq_tau(r) ) )
    #   sum_ += prob_X_n_r__k_r_leq_tau(r) * binomial(k,r) * q**r * (1-q)**(k-r)
    
    # sum_ = prob_X_n_r__k_r_leq_tau((math.ceil(k*q) + math.floor(k*q) )/2) - prob_X_n_r__k_r_leq_tau(k)*q**k # math.ceil(k*q)
    sum_ = prob_X_n_r__k_r_leq_tau(k*q) - prob_X_n_r__k_r_leq_tau(k)*q**k # math.ceil(k*q)
    return (1 - q**k - sum_) # / q**k
  # print("lhs= {}, rhs= {}".format(lhs(), rhs() ) )
  return lhs() + rhs()

# ##################  Send n initially any k is enough, X_i ~ Exp(mu), each packet drops ~ Exp(gamma)  ################ #
def Pr_succ_n_k_w_drop(mu, gamma, n, k):
  p = mu/(gamma+mu)
  def prob_F_f(f):
    return binomial(f+k-1, f) * p**k * (1-p)**f
  prob = 0
  for f in range(n-k+1):
    prob += prob_F_f(f)
  return prob

def E_T_n_k_w_drop_given_succ(mu, gamma, n, k):
  p = mu/(gamma+mu)
  def prob_F_f(f):
    return binomial(f+k-1, f) * p**k * (1-p)**f
  
  E_T = 0
  for f in range(n-k+1):
    E_T += 1/mu*(H(n) - H(n-k-f-1) ) * prob_F_f(f)
  return E_T/Pr_succ_n_k_w_drop(mu, gamma, n, k)

def E_T_n_k_w_drop_given_succ_approx(mu, gamma, n, k):
  p = mu/(gamma+mu)
  E_T = 1/mu*H(n)*Pr_succ_n_k_w_drop(mu, gamma, n, k)
  # if n*p-k > 0:
  #   E_T -= 1/mu*(math.log(n*p-k) + 0.5772156649 + 1/2/(n*p-k) )
  #     # - (math.ceil(n*p)-k > 0)*1/mu*H(math.ceil(n*p)-k)
  
  rhs = 0
  def prob_F_f(f):
    return binomial(f+k-1, f) * p**k * (1-p)**f
  # for f in range(n-k+1):
  #   rhs += H(n-k-f-1)*prob_F_f(f)
  
  for i in range(1, n-k):
    prob = 0
    for f in range(n-k-i):
      prob += prob_F_f(f)
    rhs += 1/i * prob
    # rhs += 1/i * sum([binomial(n,j) * p**j * (1-p)**(n-j) for j in range(k+i+1, n+1) ] )
  
  for i in range(1, n-k):
    rhs += 1/i * sum([binomial(n,j) * p**j * (1-p)**(n-j) for j in range(k+i+1, n+1) ] )
  
  # for j in range(k+2, n+1):
  #   rhs += H(j-k-1) * binomial(n,j) * p**j * (1-p)**(n-j)
  
  E_T -= 1/mu*rhs
  
  return E_T/Pr_succ_n_k_w_drop(mu, gamma, n, k)

# ##################  X_i ~ G, (l=k, k, n=k+1, \Delta)  ################ #
def E_T_G_1red(task_t_rv, d, k):
  E_X_k_k = mpmath.quad(lambda x: x * k*task_t_rv.cdf(x)**(k-1)*task_t_rv.pdf(x), [0, mpmath.inf] )
  return E_X_k_k - \
        mpmath.quad(lambda x: k*task_t_rv.cdf(x-d)*(1-task_t_rv.cdf(x))*task_t_rv.cdf(x)**(k-1), [d, mpmath.inf] )
  
  # def Pr_T_g_t(t):
  #   t_ = max(d, t)
  #   return (t <= d)*(task_t_rv.cdf(d)**k - task_t_rv.cdf(t)**k) \
  #           + 1 - task_t_rv.cdf(t_)**(k-1) * (k*task_t_rv.cdf(t-d)*(1 - task_t_rv.cdf(t_) ) + task_t_rv.cdf(t_) )
  # return mpmath.quad(Pr_T_g_t, [0, mpmath.inf] )

def E_T_G_1red_approx(task_t_rv, d, k):
  # mu = task_t_rv.mean()
  # return task_t_rv.cdf(d+mu) + task_t_rv.var()/2 * task_t_rv.dpdf_dx(d+mu)
  E_X_k_k = mpmath.quad(lambda x: x * k*task_t_rv.cdf(x)**(k-1)*task_t_rv.pdf(x), [0, mpmath.inf] )
  # return E_X_k_k - \
  #       mpmath.quad(lambda x: k*task_t_rv.cdf(x-d)*(1-task_t_rv.cdf(x))*task_t_rv.cdf(x)**(k-1), [d, mpmath.inf] )

def E_C_G_1red(task_t_rv, d, k, w_cancel=True):
  mu = task_t_rv.mean()
  if w_cancel:
    def Pr_T_g_t(t):
      t_ = max(d, t)
      return (t <= d)*(task_t_rv.cdf(d)**k - task_t_rv.cdf(t)**k) \
              + 1 - task_t_rv.cdf(t_)**(k-1) * (k*task_t_rv.cdf(t-d)*(1 - task_t_rv.cdf(t_) ) + task_t_rv.cdf(t_) )
    
    def E_X_k__k_minus_1(k_):
      return k_*mu - mpmath.quad(lambda x: x * k_*task_t_rv.cdf(x)**(k_-1)*task_t_rv.pdf(x), [0, mpmath.inf] )
    
    # E_X_k_k = mpmath.quad(lambda x: x * k*task_t_rv.cdf(x)**(k-1)*task_t_rv.pdf(x), [0, mpmath.inf] )
    E_T = E_T_G_1red(task_t_rv, d, k)
    # return k*mu + 2*E_T_G_1red(task_t_rv, d, k) - E_X_k_k - mpmath.quad(Pr_T_g_t, [0, d] )
    
    Pr_X_k_k__l__d_plus_X = mpmath.quad(lambda x: task_t_rv.cdf(d+x)**k*task_t_rv.pdf(x), [0, mpmath.inf] )
    Pr_X_k_k_m_1__l__d_plus_X = mpmath.quad(lambda x: (k*task_t_rv.cdf(d+x)**(k-1)*(1-task_t_rv.cdf(d+x)) + task_t_rv.cdf(d+x)**k)*task_t_rv.pdf(x), [0, mpmath.inf] )
    
    Pr_T_g_d = 1 - task_t_rv.cdf(d)**k
    sum_ = E_X_k__k_minus_1(k) + E_T \
           + mpmath.quad(Pr_T_g_t, [d, mpmath.inf] ) # - max(E_T-d-mu, 0)*(1-Pr_X_k_k__l__d_plus_X)
    return min(sum_, E_X_k__k_minus_1(k+1) + E_T_G_1red(task_t_rv, 0, k) )
  else:
    return k*mu + mu*(1 - task_t_rv.cdf(d)**k)

def plot_Pr_T_g_t_G_1red():
  def Pr_T_g_t(k, d, t):
    t_ = max(d, t)
    return (t <= d)*(task_t_rv.cdf(d)**k - task_t_rv.cdf(t)**k) \
            + 1 - task_t_rv.cdf(t_)**(k-1) * (k*task_t_rv.cdf(t-d)*(1 - task_t_rv.cdf(t_) ) + task_t_rv.cdf(t_) )
  task_t_rv = Exp(mu=2/3)
  k = 7
  D = task_t_rv.mean()
  
  d__x_l_m, d__y_l_m = {}, {}
  for d in numpy.arange(D, 5*D, D):
    d__x_l_m[d] = []
    d__y_l_m[d] = []
    for t in numpy.arange(0, 5*task_t_rv.mean(), 0.1):
      d__x_l_m[d].append(t)
      d__y_l_m[d].append(Pr_T_g_t(k, d, t) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  color = iter(cm.rainbow(numpy.linspace(0, 1, 6) ) )
  for d in numpy.arange(D, 5*D, D):
    m = next(marker)
    c = next(color)
    plot.plot(d__x_l_m[d], d__y_l_m[d], label=r'd:{}'.format(d), color=c, marker=m, linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$t$ (s)')
  plot.ylabel(r'$Pr\{T_+ > t\}$')
  plot.title(r'$k$= {}, $\Delta$= {}, $X_i \sim$ {}'.format(k, d, task_t_rv) )
  plot.savefig("plot_Pr_T_g_t_G_1red_{}.png".format(task_t_rv) )
  plot.gcf().clear()
  log(WARNING, "done; k= {}, d= {}, task_t_rv= {}".format(k, d, task_t_rv) )
  
# ####################  X_i ~ Exp(mu), (l, k, n, \Delta)  ####################### #
def E_C_exp_k_l_n(mu, d, k, l, n, w_cancel=False):
  if w_cancel:
    return k/mu
  q = 1 - math.exp(-mu*d)
  a = I(q,k,l-k+1)
  return l/mu*a + n/mu*(1-a)

def E_T_exp_k_l_n(mu, d, k, l, n):
  q = 1 - math.exp(-mu*d)
  
  # E_T = d - mpmath.quad(lambda x: sum([binomial(l,i)*(1-math.exp(-mu*x) )**i * math.exp(-mu*x)**(l-i) for i in range(k, l+1) ] ), [0, d] )
  # for r in range(k):
  #   E_T += (H(n-r) - H(n-k) ) * binomial(l,r) * q**r * (1-q)**(l-r)
  # return E_T
  E_T = H(l)-H(l-k) + I(1-q, l-k+1, k)*(H(l-k)-H(n-k) )
  for r in range(k):
    E_T += (H(n-r) - H(l-r) ) * binomial(l,r) * q**r * (1-q)**(l-r)
  return E_T

def E_T_exp_k_l_n_approx(mu, d, k, l, n):
  q = 1 - math.exp(-mu*d)
  
  # E_T = 0
  # for r in range(k+1):
  #   E_T += (H(n) - H(n-r) ) * binomial(k,r) * q**r * (1-q)**(k-r)
  
  # E_T = d - mpmath.quad(lambda x: sum([binomial(l,i)*(1-math.exp(-mu*x) )**i * math.exp(-mu*x)**(l-i) for i in range(k, l+1) ] ), [0, d] )
  
  # E_T = d - 1/mu*sum([binomial(l,i)*B(i+1, l-i, u_l=q) for i in range(k, l+1) ] )
  # E_T = d + 1/mu*math.log(1-q)*I(q,k,l-k+1) - 1/mu/B(k,l-k+1)*mpmath.quad(lambda x: math.log(1-x)*x**(k-1) * (1-x)**(l-k), [0, q] )
  # E_T = d + 1/mu*math.log(1-q)*I(q,k,l-k+1) - 1/mu/B(k,l-k+1)*sum([-1/i * B(k+i,l-k+1,u_l=q) for i in range(1, 100) ] )
  # E_T = d + 1/mu*math.log(1-q)*I(q,k,l-k+1) - 1/mu/B(k,l-k+1)*sum([-1/i * B(k+i,l-k+1,u_l=q) for i in range(1, 2) ] )
  
  # for r in range(k):
  #   E_T += (H(n-r) - H(n-k) ) * binomial(l,r) * q**r * (1-q)**(l-r)
  # sum_ = 0
  # for r in range(k):
  #   sum_ += H(n-r) * binomial(l,r) * q**r * (1-q)**(l-r)
  # for r in range(k):
  #   sum_ += H(n-r) * binomial(k,r) * q**r * (1-q)**(k-r) * (1-q)**(l-k)
  # for r in range(k):
  #   sum_ -= H(n-k) * binomial(l,r) * q**r * (1-q)**(l-r)
  # return E_T + sum_
  
  # for r in range(k):
  #   E_T += (H(n-r) - H(n-k) ) * binomial(k,r) * q**r * (1-q)**(k-r)
  
  # return E_T + (H(n) - H(n-k) )*I(1-q,l-k,k+1) + \
  #         -(l*q*1/n*I(1-q,l-k,k) )
  # #       -(l*q*(1/n+1/2/n**2)*I(1-q,l-k,k) + binomial(l,2)*(q/n)**2 * I(1-q, l-k, k-1) )
  # #       # -(l*q*1/n*I(1-q,l-k,k) )
  
  # return E_T + math.log((n-l*q)/(n-k) ) # works well for large l, k
  
  # Did not turn out to be very good
  # for r in range(k):
  #   sum_ += math.log((n-r)/(l-r) ) * binomial(l,r) * q**r * (1-q)**(l-r)
  # return math.log(l/(l-k) ) + sum_ + I(1-q, l-k+1, k)*math.log((l-k)/(n-k) )
  
  # return (math.ceil(l*q) == l)*(H(l) - H(l-k) ) + (math.ceil(l*q) != l)*(d + 1/mu*math.log(math.exp(-mu*d)*l/(n-k) + (n-l)/(n-k) ) )
  # return q**k * (H(l) - H(l-k) ) + (1-q**k)*(d + 1/mu*math.log(math.exp(-mu*d)*l/(n-k) + (n-l)/(n-k) ) )
  # return E_T + I(1-q, l-k+1, k)*(1/mu*math.log(math.exp(-mu*d)*l/(n-k) + (n-l)/(n-k) ) )
  
  E_T = H(l)-H(l-k) + I(1-q, l-k+1, k)*(H(l-k)-H(n-k) )
  # for r in range(k):
  #   E_T += math.log((n-r)/(n-k) ) * binomial(l,r) * q**r * (1-q)**(l-r)
  r_ = k*q # math.ceil(k*q)
  # E_T += H(n-r_) - H(l-r_)
  E_T += math.log(n-r_) - math.log(l-r_)
  return E_T

# ####################  X_i ~ Exp(mu), (l=k, k, n, \Delta)  ####################### #
def E_T_exp_k_n(mu, d, k, n):
  if d == 0:
    return 1/mu*(H(n) - H(n-k) )
  q = 1 - math.exp(-mu*d)
  
  E_H_n_r = 0
  for r in range(k+1):
    E_H_n_r += H(n-r) * binomial(k,r) * q**r * (1-q)**(k-r)
  
  return d - mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
         1/mu*(E_H_n_r - H(n-k) )
        # 1/mu*sum_

def E_T_exp_k_n_approx(mu, d, k, n):
  q = 1-math.exp(-mu*d)
  E_H_n_r = 0
  for r in range(k+1):
    E_H_n_r += H(n-r) * binomial(k, r) * q**r * (1-q)**(k-r)
  
  return d - 1/mu*B(k+1,0,u_l=q) + \
         1/mu*(E_H_n_r - H(n-k) )
        # 1/mu*(math.log((n-k*q)/(n-k) ) )
        # 1/mu*(math.log((n-k*q)/(n-k) ) + 1/2/(n-k*q) - 1/2/(n-k) )

def E_C_exp_k_n(mu, d, k, n, w_cancel=False):
  return E_C_exp_k_l_n(mu, d, k, k, n, w_cancel=w_cancel)

# ####################  X_i ~ D/k + Exp(mu), (l, k, n, \Delta)  ####################### #
def E_T_shiftedexp_k_l_n(D, mu, d, k, l, n):
  return D/k + E_T_exp_k_l_n(mu, d, k, l, n)

def E_C_shiftedexp_k_l_n(D, mu, d, k, l, n, w_cancel=False):
  # return n*D/k + E_C_exp_k_l_n(mu, d, k, l, n, w_cancel=w_cancel)
  if l == k:
    if not w_cancel:
      q = (d > D/k)*(1 - math.exp(-mu*(d-D/k) ) )
      # return E_C_exp_k_l_n(mu, d, k, l, n, w_cancel=False) + D + (1-q**k)*(n-k)*D/k
      return k*(1/mu + D/k)*q**k + n*(1/mu + D/k)*(1-q**k)
    else:
      if d == 0:
        return n/k*D + k/mu
      F_d = 1 - math.exp(-mu*d)
      F_d__D_over_k = (d > D/k)*(1 - math.exp(-mu*(d-D/k) ) )
      
      q = F_d__D_over_k
      a = 1 - math.exp(-mu*D/k)
      k_ = k - k*q
      E_T_ = 1/mu * a**(-k_) * B(k_+1,0,u_l=a)
      return k*(1/mu + D/k)*F_d__D_over_k**k + n*(1/mu + D/k)*(1-F_d__D_over_k**k) \
             - (n-k)*((1-F_d__D_over_k**k)/mu + E_T_*(F_d**k-F_d__D_over_k**k) )
      # - ((1-F_d__D_over_k**k)*(n-k)/mu + (n-k)*(D/2/k)*(F_d**k-F_d__D_over_k**k) )
    

# ####################  X_i ~ D/k + Exp(mu), (l=k, k, n, \Delta)  ####################### #
def E_T_shiftedexp_k_n(D, mu, d, k, n):
  return D/k + E_T_exp_k_n(mu, d, k, n)

def d_E_T_shiftedexp_k_n_dk(D, mu, d, k, n):
  q_ = 1 - math.exp(-mu*(d+D/k) )
  # return D/k**2 * (-2 + q_**k - k*(1-q_)/(n-k*q_) )
  
  # Beta_q_k_0 = mpmath.quad(lambda x: x**(k-1) * 1/(1-x), [0, q_] )
  # rhs = mu*D/k**2 * q_**k - k*Beta_q_k_0 + (mu*D/k*(1-q_) - q_)/(n-k*q_) + 1/(n-k)
  # return -2*D/k**2 + 1/mu*rhs
  
  Beta_q_ = mpmath.quad(lambda x: x**k * 1/(1-x), [0, q_] )
  rhs = (mu*D*(k+1)/k**2 * (1-q_)/q_ - math.log(q_) )*Beta_q_ + (mu*D/k*(1-q_) - q_)/(n-k*q_) + 1/(n-k)
  r = -2*D/k**2 + 1/mu*rhs
  print("k= {}, r= {}".format(k, r) )
  return r

  # q_k = 1 - math.exp(-mu*(d+D/k) )
  # q_k_1 = 1 - math.exp(-mu*(d+D/(k-1) ) )
  # # B_diff = mpmath.quad(lambda x: x**(k-1) * 1/(1-x), [0, q_k_1] ) - \
  # #         mpmath.quad(lambda x: x**k * 1/(1-x), [0, q_k] )
  # B_diff = q_k**k/k # q_k_1**k/k
  # r = 2*D*(1/k - 1/(k-1) ) + 1/mu*(B_diff + math.log((n-k*q_k)/(n-(k-1)*q_k_1) ) + 1/(n-k+1) )
  # print("k= {}, r= {}".format(k, r) )
  # return r

def E_C_shiftedexp_k_n(D, mu, d, k, n, w_cancel=False):
  return E_C_shiftedexp_k_l_n(mu, d, k, k, n, w_cancel)

# ####################  X_i ~ Pareto(loc, a), (l=k, k, n=k, \Delta)  ####################### #
def E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=True):
  if d == 0:
    # return loc*math.factorial(n)/math.factorial(n-k)*G(n-k+1-1/a)/G(n+1-1/a)
    return loc*(G(n+1)/G(n-k+1) )*(G(n-k+1-1/a)/G(n+1-1/a) )
  """
  q = (d > loc)*(1 - (loc/d)**a)
  if n == k:
    if not w_relaunch:
      return loc*G(1-1/a)*G(k+1)/G(k+1-1/a)
    else:
      if d <= loc:
        return d + E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=False)
      else:
        return d*(1-q**k) + \
               (loc/d-1)*loc*G(1-1/a)*G(k+1)/G(k+1-1/a)*I(1-q,1-1/a,k) + \
               E_T_pareto_k_n(loc, a, d, k, n, w_relaunch=False)
  else:
    if w_relaunch:
      sum_ = 0
      for r in range(k):
        sum_ += (loc*G(n-r+1)*G(n-k+1-1/a)/G(n-r+1-1/a)/G(n-k+1) - \
                 loc*G(k-r+1)*G(1-1/a)/G(k-r+1-1/a) ) * binomial(k,r) * q**r * (1-q)**(k-r)
      return sum_ + E_T_pareto_k_n(loc, a, d, k, n=k)
  """

def E_T_pareto_k_n_approx(loc, a, d, k, n, w_relaunch=True):
  q = (d > loc)*(1 - (loc/d)**a)
  if n == k:
    pass
  else:
    if w_relaunch:
      sum_ = 0
      for r in range(k):
        sum_ += (loc*G(n-k+1-1/a)/G(n-k+1)*(n-r+1)**(1/a) - \
                loc*G(1-1/a)*(k-r+1)**(1/a) ) * binomial(k,r) * q**r * (1-q)**(k-r)
      # E_R = k*q
      # sum_ = loc*G(n-k+1-1/a)/G(n-k+1)*(n-E_R+1)**(1/a) - \
      #       loc*G(1-1/a)*(k-E_R+1)**(1/a)
      return sum_ + E_T_pareto_k_n(loc, a, d, k, n=k)

def E_C_pareto_k_n(loc, a, d, k, n, w_relaunch=True, w_cancel=True):
  if w_cancel and d == 0:
    return loc*n/(a-1) * (a - (G(n)/G(n-k) )*(G(n-k+1-1/a)/G(n+1-1/a) ) )
    
    # E_C = 0
    # a_ = 1/a
    # Q = loc*(G(n+1)/G(n+1-a_) )
    # for i in range(1, k+1):
    #   E_C += loc*(G(n+1)/G(n+1-i) )*(G(n-i-1/a+1)/G(n-1/a+1) )
    # return E_C + (n-k)*Q*(G(n-k+1-a_)/G(n-k+1) )

def plot_pareto_zerodelay_red():
  loc = 3
  d = 0
  def reduction_in_E_T_for_same_E_C_w_red(red_type, loc, a, k):
    E_T_wo_red, E_T = 0, 0
    if red_type == 'coded':
      E_T_wo_red = E_T_pareto_k_n(loc, a, d, k, n=k)
      E_C_wo_red = E_C_pareto_k_n(loc, a, d, k, n=k)
      
      n_ = k+1
      while 1:
        E_C = E_C_pareto_k_n(loc, a, d, k, n_)
        # print("n_= {}, E_C_wo_red= {}, E_C= {}".format(n_, E_C_wo_red, E_C) )
        if math.isnan(E_C):
          return reduction_in_E_T_for_same_E_C_w_red__approx(red_type, loc, a, k)
        elif E_C >= E_C_wo_red:
          # print("breaking at n_= {}".format(n_) )
          break
        # print("n_= {}".format(n_) )
        n_ += 1
      E_T = E_T_pareto_k_n(loc, a, d, k, n_-1)
    elif red_type == 'reped':
      E_T_wo_red = E_T_pareto_k_c(loc, a, d, k, c=0)
      E_C_wo_red = E_C_pareto_k_c(loc, a, d, k, c=0)
      
      c_ = 1
      while 1:
        E_C = E_C_pareto_k_c(loc, a, d, k, c_)
        # print("c_= {}, E_C_wo_red= {}, E_C= {}".format(c_, E_C_wo_red, E_C) )
        if math.isnan(E_C) or math.isinf(E_C):
          return None
        elif E_C >= E_C_wo_red:
          # print("breaking at c_= {}".format(c_) )
          break
        c_ += 1
      E_T = E_T_pareto_k_c(loc, a, d, k, c_-1)
    # return E_T_wo_red - E_T
    return (E_T_wo_red - E_T)/E_T_wo_red
  
  def reduction_in_E_T_for_same_E_C_w_red__approx(red_type, loc, a, k):
    if red_type == 'coded':
      E_T_wo_red = E_T_pareto_k_n(loc, a, d, k, n=k)
      # return E_T_wo_red - loc*a
      return max(E_T_wo_red - loc*a, 0)/E_T_wo_red
    elif red_type == 'reped':
      E_T_wo_red = E_T_pareto_k_c(loc, a, d, k, c=0)
      # E_T = loc*math.factorial(k)*G(1/a)/G(k+1/a)
      c_ = max(math.floor(1/(a-1) - 1), 0)
      E_T = E_T_pareto_k_c(loc, a, d, k, c_) # exact
      return max(E_T_wo_red - E_T, 0)/E_T_wo_red
  
  # def threshold_on_tail_for_latency_reduction_at_no_cost(red_type, a):
  #   if red_type == 'reped':
  #     return 1.5
  #   elif red_type == 'coded':
  #     return (a + math.factorial(k)*G(1-1/a)/G(k+1-1/a) )**(1/k)
  
  def plot_(red_type, k):
    x_l, y_l, y_approx_l = [], [], []
    for a in numpy.arange(1.05, 3.5, 0.05):
      x_l.append(a)
      y_l.append(reduction_in_E_T_for_same_E_C_w_red(red_type, loc, a, k) )
      # y_approx_l.append(reduction_in_E_T_for_same_E_C_w_red__approx(red_type, loc, a, k) )
    # plot.axvline(x=threshold_on_tail_for_latency_reduction_at_no_cost(red_type),
    #             color=color, alpha=0.5, linestyle='--')
    legend = "Replicated" if red_type == 'reped' else "Coded"
    plot.plot(x_l, y_l, label=r'{},$k={}$'.format(legend, k), color=next(dark_color), marker=next(marker), mew=2, zorder=0, linestyle=':')
    # plot.plot(x_l, y_approx_l, label=r'{},$k={}$, approx'.format(legend, k), color=next(dark_color), marker=next(marker), mew=2, zorder=1, linestyle='')
  
  plot_('reped', k=10)
  # plot_('reped', k=20)
  # plot_('reped', k=40)
  plot_('reped', k=50)
  
  plot_('coded', k=10)
  # plot_('coded', k=20)
  # plot_('coded', k=40)
  plot_('coded', k=50)
  
  # plot_('coded', k=100)
  
  plot.legend()
  plot.title(r'$X \sim Pareto(\lambda={}, \alpha)$'.format(loc) )
  plot.xlabel(r'Tail index $\alpha$', fontsize=12)
  plot.ylabel('\n'.join(textwrap.wrap(r'Maximum percetange reduction in $E[T]$ at no added cost', 40) ),
              fontsize=12)
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  plot.savefig("plot_pareto_zerodelay_red.pdf", bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done.")
  

# ####################  Wrappers, (l, k, n, \Delta)  ####################### #
def E_T_k_l_n(task_t, D, mu, loc, a, d, k, l, n):
  if task_t == "Exp":
    if l == k: return E_T_exp_k_n(mu, d, k, n)
    else: return E_T_exp_k_l_n(mu, d, k, l, n)
  elif task_t == "SExp":
    if l == k: return E_T_shiftedexp_k_n(D, mu, d, k, n)
    else: return E_T_shiftedexp_k_l_n(D, mu, d, k, l, n)
  elif task_t == "Pareto":
    return E_T_pareto_k_n(loc, a, d, k, n)

def E_C_k_l_n(task_t, D, mu, loc, a, d, k, l, n, w_cancel, approx=False):
  if task_t == "Exp":
    if l == k: return E_C_exp_k_n(mu, d, k, l, n, w_cancel)
    else: return E_C_exp_k_l_n(mu, d, k, l, n, w_cancel)
  elif task_t == "SExp":
    return E_C_shiftedexp_k_l_n(D, mu, d, k, l, n, w_cancel)
  elif task_t == "Pareto":
    return E_C_pareto_k_n(loc, a, d, k, n, w_cancel=w_cancel)

# ####################  X_i ~ Exp(mu), (k, \Delta, c)  ####################### #
def E_T_exp_k_c(mu, d, k, c):
  if d == 0:
    return 1/mu/(c+1)*H(k)
  q = 1 - math.exp(-mu*d)
  
  E = 0
  for r in range(k+1):
    E += H(k-r)* binomial(k,r) * q**r * (1-q)**(k-r)
  return 1/mu*(H(k) - c/(c+1)*E)

def E_C_exp_k_c(mu, d, k, c, w_cancel=False):
  if w_cancel:
    return k/mu
  q = 1 - math.exp(-mu*d)
  return k/mu*(c*(1-q) + 1)

def E_T_exp_k_c_approx(mu, d, k):
  q = 1 - math.exp(-mu*d)
  # return 1/2/mu*(H(k) - math.log(1-q) )
  
  # return 1/2/mu*(2*H(k) - H(math.ceil(k-k*q) ) )
  # H_k_kq = math.log(n) + 0.5772156649 + 1/2/n
  H_k_kq = H_cont(k-k*q)
  return 1/mu*(H(k) - c/(c+1)*H_k_kq)

# ####################  X_i ~ D/k + Exp(mu), (k, \Delta, c)  ####################### #
def E_T_shiftedexp_k_c(D, mu, d, k, c):
  return D/k + E_T_exp_k_c(mu, d, k, c)

def E_C_shiftedexp_k_c(D, mu, d, k, c, w_cancel=False):
  q = (d > D/k)*(1 - math.exp(-mu*(d - D/k) ) )
  if not w_cancel:
    return (c*(1-q)+1)*(D + k/mu)
  else:
    if d == 0:
      return (c+1)*D + k/mu
    elif d <= D/k:
      c_ = c/(c+1)
      return (c+1)*(D + k*(1/mu*(1-c_*math.exp(-mu*d)) - c_*d) )
      # 2*D + k*((2-math.exp(-mu*d) )/mu - d)
    else:
      return D + k/mu*(1 + c*(1-q-math.exp(-mu*d)) )
      # D + k/mu*(2-q-math.exp(-mu*d) )

# ####################  X_i ~ Pareto(loc, a), (k, \Delta, c)  ####################### #
def E_T_pareto_k_c(loc, a, d, k, c):
  if d == 0:
    # return loc*math.factorial(k)*G(1-1/(c+1)/a)/G(k+1-1/(c+1)/a)
    return loc*G(k+1)*G(1-1/(c+1)/a)/G(k+1-1/(c+1)/a)

def E_C_pareto_k_c(loc, a, d, k, c, w_cancel=True):
  if w_cancel and d == 0:
    # log(WARNING, "loc= {}, a= {}, d= {}, k= {}, c= {}".format(loc, a, d, k, c) )
    return k*(c+1) * a*(c+1)*loc/(a*(c+1)-1)
  
  q = 0 if d <= loc else 1 - (loc/d)**a
  if not w_cancel:
    return (c*(1-q)+1)*(a*loc/(a-1) )
  else:
    def tail(x):
      if x <= loc: return 1
      else: return (loc/x)**a
    def proxy(x): return tail(x)*tail(x+d)/tail(d)
    if c == 1:
      if d <= loc:
        E_X__X_leq_d = 0
        E_Y = mpmath.quad(proxy, [0, mpmath.inf] )
        return k*(q*E_X__X_leq_d + (1-q)*(2*E_Y + d) )
      else:
        E_X__X_leq_d = loc*a/(a-1)*(1 - (loc/d)**(a-1) )/(1 - (loc/d)**a)
        E_Y = mpmath.quad(proxy, [0, mpmath.inf] )
        return k*(q*E_X__X_leq_d + (1-q)*(2*E_Y + d) )

def E_C_pareto_k_c_approx(loc, a, d, k, w_cancel=True):
  q = 0 if d <= l else 1 - (l/d)**a
  def tail(x):
    if x <= l: return 1
    else: return (l/x)**a
  def tail_delayed(x):
    return (d/(x+d) )**a
  def proxy(x): return tail(x)*tail(x+d)/tail(d)
  def proxy_u(x): return tail(x)**2
  def proxy_l(x): return tail(x+d)**2/tail(d)**2
  if c == 1:
    if d <= l:
      E_X__X_leq_d = 0
      # E_Y = (a*l + d)/2/(a-1)
      # E_Y = l-d + l**a*(l**(1-a) - (l+d)**(1-a) )/(a-1) \
      #       + l**(2*a) * (l*(l+d))**(a-1/2) / (2*a-1)
      
      # E_Y = l-d + l**a*mpmath.quad(lambda x: (x+d)**(-a), [l-d, l] ) \
      #       + l**(2*a)*mpmath.quad(lambda x: x**(-a) * (x+d)**(-a), [l, mpmath.inf] )
      # E_Y = l-d + l**a*mpmath.quad(lambda x: (x+d)**(-a), [l-d, l] ) \
      #       + (l**(2*a)*mpmath.quad(lambda x: x**(-a) * x**(-a), [l, mpmath.inf] ) \
      #         + l**(2*a)*mpmath.quad(lambda x: (x+d)**(-a) * (x+d)**(-a), [l, mpmath.inf] ) )/2
      
      # E_Y = (mpmath.quad(proxy_u, [0, mpmath.inf] ) + mpmath.quad(proxy_l, [0, mpmath.inf] ) )/2
      E_Y = mpmath.quad(tail_delayed, [0, mpmath.inf] )
      
      return k*(q*E_X__X_leq_d + (1-q)*(2*E_Y + d) )
    else:
      E_X__X_leq_d = l*a/(a-1)*(1 - (l/d)**(a-1) )/(1 - (l/d)**a)
      # E_Y = (a*l + d)/2/(a-1)
      # E_Y = mpmath.quad(proxy, [0, mpmath.inf] )
      
      # E_Y = d**a*mpmath.quad(lambda x: (x+d)**(-a), [0, l] ) \
      #       + (l*d)**a*mpmath.quad(lambda x: x**(-a) * (x+d)**(-a), [l, mpmath.inf] )
      # E_Y = d**a*mpmath.quad(lambda x: (x+d)**(-a), [0, l] ) \
      #       + ((l*d)**a*mpmath.quad(lambda x: x**(-a) * x**(-a), [l, mpmath.inf] )
      #         + (l*d)**a*mpmath.quad(lambda x: (x+d)**(-a) * (x+d)**(-a), [l, mpmath.inf] ) )/2
      
      # E_Y = (mpmath.quad(proxy_u, [0, mpmath.inf] ) + mpmath.quad(proxy_l, [0, mpmath.inf] ) )/2
      E_Y = mpmath.quad(tail, [0, mpmath.inf] )
      
      return k*(q*E_X__X_leq_d + (1-q)*(2*E_Y + d) )

# ####################  Wrappers, (k, \Delta, c)  ####################### #
def E_T_k_c(task_t, D, mu, loc, a, d, k, c):
  if task_t == "Exp":
    return E_T_exp_k_c(mu, d, k, c)
  elif task_t == "SExp":
    return E_T_shiftedexp_k_c(D, mu, d, k, c)
  elif task_t == "Pareto":
    return E_T_pareto_k_c(loc, a, d, k, c)

def E_C_k_c(task_t, D, mu, loc, a, d, k, c, w_cancel, approx=False):
  if task_t == "Exp":
    return E_C_exp_k_c(mu, d, k, c, w_cancel)
  elif task_t == "SExp":
    return E_C_shiftedexp_k_c(D, mu, d, k, c, w_cancel)
  elif task_t == "Pareto":
    if approx:
      return E_C_pareto_k_c_approx(loc, a, d, k, c, w_cancel)
    return E_C_pareto_k_c(loc, a, d, k, c, w_cancel)

if __name__ == "__main__":
  # plot_Pr_T_g_t_G_1red()
  plot_pareto_zerodelay_red()