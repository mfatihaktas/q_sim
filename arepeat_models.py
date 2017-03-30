import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, random, math, numpy, simpy, getopt, itertools, sympy

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

def approx_prob_T_k_n_geq_t(mu, d, k, n, t):
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

def approx_E_T_n_k_w_drop_given_succ(mu, gamma, n, k):
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
  E_X_k_k = sympy.mpmath.quad(lambda x: x * k*task_t_rv.cdf(x)**(k-1)*task_t_rv.pdf(x), [0, sympy.mpmath.inf] )
  return E_X_k_k - \
        sympy.mpmath.quad(lambda x: k*task_t_rv.cdf(x-d)*(1-task_t_rv.cdf(x))*task_t_rv.cdf(x)**(k-1), [d, sympy.mpmath.inf] )
  
  # def Pr_T_g_t(t):
  #   t_ = max(d, t)
  #   return (t <= d)*(task_t_rv.cdf(d)**k - task_t_rv.cdf(t)**k) \
  #           + 1 - task_t_rv.cdf(t_)**(k-1) * (k*task_t_rv.cdf(t-d)*(1 - task_t_rv.cdf(t_) ) + task_t_rv.cdf(t_) )
  # return sympy.mpmath.quad(Pr_T_g_t, [0, sympy.mpmath.inf] )

def E_T_G_1red_approx(task_t_rv, d, k):
  # mu = task_t_rv.mean()
  # return task_t_rv.cdf(d+mu) + task_t_rv.var()/2 * task_t_rv.dpdf_dx(d+mu)
  E_X_k_k = sympy.mpmath.quad(lambda x: x * k*task_t_rv.cdf(x)**(k-1)*task_t_rv.pdf(x), [0, sympy.mpmath.inf] )
  # return E_X_k_k - \
  #       sympy.mpmath.quad(lambda x: k*task_t_rv.cdf(x-d)*(1-task_t_rv.cdf(x))*task_t_rv.cdf(x)**(k-1), [d, sympy.mpmath.inf] )

def E_C_G_1red(task_t_rv, d, k, w_cancel=True):
  mu = task_t_rv.mean()
  if w_cancel:
    def Pr_T_g_t(t):
      t_ = max(d, t)
      return (t <= d)*(task_t_rv.cdf(d)**k - task_t_rv.cdf(t)**k) \
              + 1 - task_t_rv.cdf(t_)**(k-1) * (k*task_t_rv.cdf(t-d)*(1 - task_t_rv.cdf(t_) ) + task_t_rv.cdf(t_) )
    
    def E_X_k__k_minus_1(k_):
      return k_*mu - sympy.mpmath.quad(lambda x: x * k_*task_t_rv.cdf(x)**(k_-1)*task_t_rv.pdf(x), [0, sympy.mpmath.inf] )
    
    # E_X_k_k = sympy.mpmath.quad(lambda x: x * k*task_t_rv.cdf(x)**(k-1)*task_t_rv.pdf(x), [0, sympy.mpmath.inf] )
    E_T = E_T_G_1red(task_t_rv, d, k)
    # return k*mu + 2*E_T_G_1red(task_t_rv, d, k) - E_X_k_k - sympy.mpmath.quad(Pr_T_g_t, [0, d] )
    
    Pr_X_k_k__l__d_plus_X = sympy.mpmath.quad(lambda x: task_t_rv.cdf(d+x)**k*task_t_rv.pdf(x), [0, sympy.mpmath.inf] )
    Pr_X_k_k_m_1__l__d_plus_X = sympy.mpmath.quad(lambda x: (k*task_t_rv.cdf(d+x)**(k-1)*(1-task_t_rv.cdf(d+x)) + task_t_rv.cdf(d+x)**k)*task_t_rv.pdf(x), [0, sympy.mpmath.inf] )
    
    Pr_T_g_d = 1 - task_t_rv.cdf(d)**k
    sum_ = E_X_k__k_minus_1(k) + E_T \
           + sympy.mpmath.quad(Pr_T_g_t, [d, sympy.mpmath.inf] ) # - max(E_T-d-mu, 0)*(1-Pr_X_k_k__l__d_plus_X)
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
  
# ####################  X_i ~ Exp(mu), (l, k, n=k+1, \Delta)  ####################### #
def E_C_exp_k_l_n(mu, d, k, l, n, w_cancel=False):
  if w_cancel:
    return k/mu
  q = 1 - math.exp(-mu*d)
  a = I(q,k,l-k+1)
  return l/mu*a + n/mu*(1-a)

def E_T_exp_k_l_n(mu, d, k, l, n):
  q = 1 - math.exp(-mu*d)
  
  # E_T = d - sympy.mpmath.quad(lambda x: sum([binomial(l,i)*(1-math.exp(-mu*x) )**i * math.exp(-mu*x)**(l-i) for i in range(k, l+1) ] ), [0, d] )
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
  
  # E_T = d - sympy.mpmath.quad(lambda x: sum([binomial(l,i)*(1-math.exp(-mu*x) )**i * math.exp(-mu*x)**(l-i) for i in range(k, l+1) ] ), [0, d] )
  
  # E_T = d - 1/mu*sum([binomial(l,i)*B(i+1, l-i, u_l=q) for i in range(k, l+1) ] )
  # E_T = d + 1/mu*math.log(1-q)*I(q,k,l-k+1) - 1/mu/B(k,l-k+1)*sympy.mpmath.quad(lambda x: math.log(1-x)*x**(k-1) * (1-x)**(l-k), [0, q] )
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

# ####################  X_i ~ Exp(mu), (l=k, k, n=k+1, \Delta)  ####################### #
def E_C_exp_k_n(mu, d, k, n, w_cancel=False):
  if w_cancel:
    return k/mu
  return E_C_exp_k_l_n(mu, d, k, k, n)

def E_T_exp_k_n(mu, d, k, n):
  q = 1 - math.exp(-mu*d)
  
  E_H_n_r = 0
  for r in range(k+1):
    E_H_n_r += H(n-r) * binomial(k,r) * q**r * (1-q)**(k-r)
  
  return d - sympy.mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
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

# ####################  X_i ~ D/k + Exp(mu), (l=k, k, n=k+1, \Delta)  ####################### #
def E_T_shiftedexp_k_n(D, mu, d, k, n):
  return D/k + E_T_exp_k_n(mu, d, k, n)

def d_E_T_shiftedexp_k_n_dk(D, mu, d, k, n):
  q_ = 1 - math.exp(-mu*(d+D/k) )
  # return D/k**2 * (-2 + q_**k - k*(1-q_)/(n-k*q_) )
  
  # Beta_q_k_0 = sympy.mpmath.quad(lambda x: x**(k-1) * 1/(1-x), [0, q_] )
  # rhs = mu*D/k**2 * q_**k - k*Beta_q_k_0 + (mu*D/k*(1-q_) - q_)/(n-k*q_) + 1/(n-k)
  # return -2*D/k**2 + 1/mu*rhs
  
  Beta_q_ = sympy.mpmath.quad(lambda x: x**k * 1/(1-x), [0, q_] )
  rhs = (mu*D*(k+1)/k**2 * (1-q_)/q_ - math.log(q_) )*Beta_q_ + (mu*D/k*(1-q_) - q_)/(n-k*q_) + 1/(n-k)
  r = -2*D/k**2 + 1/mu*rhs
  print("k= {}, r= {}".format(k, r) )
  return r

  # q_k = 1 - math.exp(-mu*(d+D/k) )
  # q_k_1 = 1 - math.exp(-mu*(d+D/(k-1) ) )
  # # B_diff = sympy.mpmath.quad(lambda x: x**(k-1) * 1/(1-x), [0, q_k_1] ) - \
  # #         sympy.mpmath.quad(lambda x: x**k * 1/(1-x), [0, q_k] )
  # B_diff = q_k**k/k # q_k_1**k/k
  # r = 2*D*(1/k - 1/(k-1) ) + 1/mu*(B_diff + math.log((n-k*q_k)/(n-(k-1)*q_k_1) ) + 1/(n-k+1) )
  # print("k= {}, r= {}".format(k, r) )
  # return r

# ####################  X_i ~ Exp(mu), (k, \Delta)  ####################### #
def E_C_exp_k(mu, d, k, w_cancel=False):
  if w_cancel:
    return k/mu
  q = 1 - math.exp(-mu*d)
  return k/mu*(2-q)

def E_T_exp_k(mu, d, k):
  q = 1 - math.exp(-mu*d)
  
  E = 0
  for r in range(k+1):
    E += H(k-r)* binomial(k,r) * q**r * (1-q)**(k-r)
  return H(k)/mu - 1/2/mu*E

def E_T_exp_k_approx(mu, d, k):
  q = 1 - math.exp(-mu*d)
  # return 1/2/mu*(H(k) - math.log(1-q) )
  
  # return 1/2/mu*(2*H(k) - H(math.ceil(k-k*q) ) )
  n = k-k*q
  # H_k_kq = math.log(n) + 0.5772156649 + 1/2/n
  H_k_kq = sympy.mpmath.quad(lambda x: (1-x**n)/(1-x), [0, 1] )
  return 1/2/mu*(2*H(k) - H_k_kq)

# ####################  X_i ~ D/k + Exp(mu), (k, \Delta)  ####################### #
def E_T_shiftedexp_k(mu, d, k):
  return D/k + E_T_exp_k(mu, d, k)

if __name__ == "__main__":
  plot_Pr_T_g_t_G_1red()