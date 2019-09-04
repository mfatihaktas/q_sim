import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import sys, pprint, random, math, numpy, getopt, itertools, mpmath, textwrap

from rvs import *
from patch import *

def plot_compare_tails():
  exp_rv = Exp(D=1, mu=0.5)
  pareto_rv = Pareto(loc=1, a=5)
  
  x_l = []
  exp_tail_l, pareto_tail_l = [], []
  for x in numpy.logspace(0, 2, 100):
    x_l.append(x)
    exp_tail_l.append(exp_rv.tail(x) )
    pareto_tail_l.append(pareto_rv.tail(x) )
  plot.plot(x_l, exp_tail_l, label=r'$Exp$', color='red', lw=2, linestyle='-')
  plot.plot(x_l, pareto_tail_l, label=r'$Pareto$', color='green', lw=2, linestyle='-')
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  plot.legend(fontsize=14)
  
  # plot.xlim([1, 21] )
  plot.xlabel(r'Task lifetime', fontsize=14)
  plot.ylabel(r'Tail distribution', fontsize=14)
  plot.xscale('log')
  plot.yscale('log')
  fig.tight_layout()
  ax = plot.gca()
  # ax.text(2, 10**-2, r'Exp: $e^{-\mu x}$', fontsize=16, color='red')
  # ax.text(2, 10**-3, r'Pareto: $(\lambda/x)^{\alpha}$ for $x \geq \lambda$', fontsize=16, color='green')
  plot.savefig("plot_compare_tails.pdf")
  fig.clear()
  log(WARNING, "done.")

# ###########################  Tail distribution  ######################### #
def Pr_T_g_t_Exp_k_n(mu, d, k, n, t):
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
        sum_ += binom(n-r,j) * q_**j * (1-q_)**(n-r-j)
      return sum_
    sum_ = 0
    for r in range(k):
      # print("prob_X_n_r__k_r_leq_tau(r= {})= {}".format(r, prob_X_n_r__k_r_leq_tau(r) ) )
      sum_ += prob_X_n_r__k_r_leq_tau(r) * binom(k,r) * q**r * (1-q)**(k-r)
    return (1 - q**k - sum_)
  # print("lhs= {}, rhs= {}".format(lhs(), rhs() ) )
  return lhs() + rhs()

def Pr_T_g_t_Exp_k_n_approx(mu, d, k, n, t):
  q = 1 - math.exp(-mu*d)
  q_ = 1 - math.exp(-mu*t)
  tau = max(0, t - d)
  q__ = 1 - math.exp(-mu*tau)
  
  sum_ = 0
  for r in range(k+1):
    sum_ += I(q__, k-r, n-k+1) * binom(k,r) * q**r * (1-q)**(k-r)
  
  return (t < d)*(q**k - q_**k) \
         + (1-q**k) - sum_ + q**k*I(q__, 0, n-k+1)
        # + (1-q**k) - I(q__, k*(1-q), n-k+1) + q**k*I(q__, 0, n-k+1) # This does not work very well unfortunately

def Pr_T_g_t_k_n(task_t, task_dist_m, d, k, n, t):
  if task_t == "Exp":
    return Pr_T_g_t_Exp_k_n(mu, d, k, n, t)
  elif task_t == "SExp":
    return None
  elif task_t == "Pareto":
    return None

# ##################  Send n initially any k is enough, X ~ Exp(mu), each packet drops ~ Exp(gamma)  ################ #
def Pr_succ_n_k_w_drop(mu, gamma, n, k):
  p = mu/(gamma+mu)
  def prob_F_f(f):
    return binom(f+k-1, f) * p**k * (1-p)**f
  prob = 0
  for f in range(n-k+1):
    prob += prob_F_f(f)
  return prob

def E_T_n_k_w_drop_given_succ(mu, gamma, n, k):
  p = mu/(gamma+mu)
  def prob_F_f(f):
    return binom(f+k-1, f) * p**k * (1-p)**f
  
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
    return binom(f+k-1, f) * p**k * (1-p)**f
  # for f in range(n-k+1):
  #   rhs += H(n-k-f-1)*prob_F_f(f)
  
  for i in range(1, n-k):
    prob = 0
    for f in range(n-k-i):
      prob += prob_F_f(f)
    rhs += 1/i * prob
    # rhs += 1/i * sum([binom(n,j) * p**j * (1-p)**(n-j) for j in range(k+i+1, n+1) ] )
  
  for i in range(1, n-k):
    rhs += 1/i * sum([binom(n,j) * p**j * (1-p)**(n-j) for j in range(k+i+1, n+1) ] )
  
  # for j in range(k+2, n+1):
  #   rhs += H(j-k-1) * binom(n,j) * p**j * (1-p)**(n-j)
  
  E_T -= 1/mu*rhs
  
  return E_T/Pr_succ_n_k_w_drop(mu, gamma, n, k)

# ##################  X ~ G, (l=k, k, n=k+1, \Delta)  ################ #
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
  
# ********************************  CODING  ******************************** #
# ####################  X ~ Exp(mu), (l, k, n, \Delta)  ####################### #
def E_C_exp_k_l_n(mu, d, k, l, n, w_cancel=False):
  if w_cancel:
    return k/mu
  q = 1 - math.exp(-mu*d)
  a = I(q,k,l-k+1)
  return l/mu*a + n/mu*(1-a)

def E_T_exp_k_l_n(mu, d, k, l, n):
  q = 1 - math.exp(-mu*d)
  
  # E_T = d - mpmath.quad(lambda x: sum([binom(l,i)*(1-math.exp(-mu*x) )**i * math.exp(-mu*x)**(l-i) for i in range(k, l+1) ] ), [0, d] )
  # for r in range(k):
  #   E_T += (H(n-r) - H(n-k) ) * binom(l,r) * q**r * (1-q)**(l-r)
  # return E_T
  E_T = H(l)-H(l-k) + I(1-q, l-k+1, k)*(H(l-k)-H(n-k) )
  for r in range(k):
    E_T += (H(n-r) - H(l-r) ) * binom(l,r) * q**r * (1-q)**(l-r)
  return E_T

def E_T_exp_k_l_n_approx(mu, d, k, l, n):
  q = 1 - math.exp(-mu*d)
  
  # E_T = 0
  # for r in range(k+1):
  #   E_T += (H(n) - H(n-r) ) * binom(k,r) * q**r * (1-q)**(k-r)
  
  # E_T = d - mpmath.quad(lambda x: sum([binom(l,i)*(1-math.exp(-mu*x) )**i * math.exp(-mu*x)**(l-i) for i in range(k, l+1) ] ), [0, d] )
  
  # E_T = d - 1/mu*sum([binom(l,i)*B(i+1, l-i, u_l=q) for i in range(k, l+1) ] )
  # E_T = d + 1/mu*math.log(1-q)*I(q,k,l-k+1) - 1/mu/B(k,l-k+1)*mpmath.quad(lambda x: math.log(1-x)*x**(k-1) * (1-x)**(l-k), [0, q] )
  # E_T = d + 1/mu*math.log(1-q)*I(q,k,l-k+1) - 1/mu/B(k,l-k+1)*sum([-1/i * B(k+i,l-k+1,u_l=q) for i in range(1, 100) ] )
  # E_T = d + 1/mu*math.log(1-q)*I(q,k,l-k+1) - 1/mu/B(k,l-k+1)*sum([-1/i * B(k+i,l-k+1,u_l=q) for i in range(1, 2) ] )
  
  # for r in range(k):
  #   E_T += (H(n-r) - H(n-k) ) * binom(l,r) * q**r * (1-q)**(l-r)
  # sum_ = 0
  # for r in range(k):
  #   sum_ += H(n-r) * binom(l,r) * q**r * (1-q)**(l-r)
  # for r in range(k):
  #   sum_ += H(n-r) * binom(k,r) * q**r * (1-q)**(k-r) * (1-q)**(l-k)
  # for r in range(k):
  #   sum_ -= H(n-k) * binom(l,r) * q**r * (1-q)**(l-r)
  # return E_T + sum_
  
  # for r in range(k):
  #   E_T += (H(n-r) - H(n-k) ) * binom(k,r) * q**r * (1-q)**(k-r)
  
  # return E_T + (H(n) - H(n-k) )*I(1-q,l-k,k+1) + \
  #         -(l*q*1/n*I(1-q,l-k,k) )
  # #       -(l*q*(1/n+1/2/n**2)*I(1-q,l-k,k) + binom(l,2)*(q/n)**2 * I(1-q, l-k, k-1) )
  # #       # -(l*q*1/n*I(1-q,l-k,k) )
  
  # return E_T + math.log((n-l*q)/(n-k) ) # works well for large l, k
  
  # Did not turn out to be very good
  # for r in range(k):
  #   sum_ += math.log((n-r)/(l-r) ) * binom(l,r) * q**r * (1-q)**(l-r)
  # return math.log(l/(l-k) ) + sum_ + I(1-q, l-k+1, k)*math.log((l-k)/(n-k) )
  
  # return (math.ceil(l*q) == l)*(H(l) - H(l-k) ) + (math.ceil(l*q) != l)*(d + 1/mu*math.log(math.exp(-mu*d)*l/(n-k) + (n-l)/(n-k) ) )
  # return q**k * (H(l) - H(l-k) ) + (1-q**k)*(d + 1/mu*math.log(math.exp(-mu*d)*l/(n-k) + (n-l)/(n-k) ) )
  # return E_T + I(1-q, l-k+1, k)*(1/mu*math.log(math.exp(-mu*d)*l/(n-k) + (n-l)/(n-k) ) )
  
  E_T = H(l)-H(l-k) + I(1-q, l-k+1, k)*(H(l-k)-H(n-k) )
  # for r in range(k):
  #   E_T += math.log((n-r)/(n-k) ) * binom(l,r) * q**r * (1-q)**(l-r)
  r_ = k*q # math.ceil(k*q)
  # E_T += H(n-r_) - H(l-r_)
  E_T += math.log(n-r_) - math.log(l-r_)
  return E_T

# ####################  X ~ Exp(mu), (l=k, k, n, \Delta)  ####################### #
def E_T_exp_k_n(mu, d, k, n):
  if d == 0:
    return 1/mu*(H(n) - H(n-k) )
  q = 1 - math.exp(-mu*d)
  
  E_H_n_r = 0
  for r in range(k+1):
    E_H_n_r += H(n-r) * binom(k,r) * q**r * (1-q)**(k-r)
  
  return d - mpmath.quad(lambda x: (1-math.exp(-mu*x) )**k, [0, d] ) + \
         1/mu*(E_H_n_r - H(n-k) )
        # 1/mu*sum_

def E_T_exp_k_n_approx(mu, d, k, n):
  q = 1-math.exp(-mu*d)
  E_H_n_r = 0
  for r in range(k+1):
    E_H_n_r += H(n-r) * binom(k, r) * q**r * (1-q)**(k-r)
  
  return d - 1/mu*B(k+1,0,u_l=q) + \
         1/mu*(E_H_n_r - H(n-k) )
        # 1/mu*(math.log((n-k*q)/(n-k) ) )
        # 1/mu*(math.log((n-k*q)/(n-k) ) + 1/2/(n-k*q) - 1/2/(n-k) )

def E_C_exp_k_n(mu, d, k, n, w_cancel=False):
  return E_C_exp_k_l_n(mu, d, k, k, n, w_cancel=w_cancel)

# ####################  X ~ D/k + Exp(mu), (l, k, n, \Delta)  ####################### #
def E_T_shiftedexp_k_l_n(D, mu, d, k, l, n):
  return D/k + E_T_exp_k_l_n(mu, d, k, l, n)

def E_C_shiftedexp_k_l_n(D, mu, d, k, l, n, w_cancel=False):
  s = D/k
  # return n*s + E_C_exp_k_l_n(mu, d, k, l, n, w_cancel=w_cancel)
  if l == k:
    if not w_cancel:
      q = (d > s)*(1 - math.exp(-mu*(d-s) ) )
      # return E_C_exp_k_l_n(mu, d, k, l, n, w_cancel=False) + D + (1-q**k)*(n-k)*s
      return k*(1/mu + s)*q**k + n*(1/mu + s)*(1-q**k)
    else:
      if d == 0:
        return n/k*D + k/mu
      elif d < s:
        q = 1 - math.exp(-mu*d)
        return n*s - q**k*(d - H(k)/mu) + k/mu
      else:
        F_d = 1 - math.exp(-mu*d)
        F_d__D_over_k = (d > s)*(1 - math.exp(-mu*(d-s) ) )
        
        q = F_d__D_over_k
        a = 1 - math.exp(-mu*s)
        k_ = k - k*q
        # print("B(k_+1, 0, u_l=a)= {}".format(B(k_+1, 0, u_l=a) ) )
        E_T_ = 1/mu * a**(-k_) * B(k_+1, 0, u_l=a)
        # print("E_T_= {}".format(E_T_) )
        return k*(1/mu + s)*F_d__D_over_k**k + n*(1/mu + s)*(1-F_d__D_over_k**k) \
               - (n-k)*((1-F_d__D_over_k**k)/mu + E_T_*(F_d**k-F_d__D_over_k**k) )
        # - ((1-F_d__D_over_k**k)*(n-k)/mu + (n-k)*(D/2/k)*(F_d**k-F_d__D_over_k**k) )

# ####################  X ~ D/k + Exp(mu), (l=k, k, n, \Delta)  ####################### #
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
  return E_C_shiftedexp_k_l_n(D, mu, d, k, k, n, w_cancel)

# ********************************  REPLICATION  ******************************** #
# ####################  X ~ Exp(mu), (k, \Delta, c)  ####################### #
def E_T_exp_k_c(mu, d, k, c):
  if d == 0:
    return 1/mu/(c+1)*H(k)
  q = 1 - math.exp(-mu*d)
  
  E = 0
  for r in range(k+1):
    E += H(k-r)* binom(k,r) * q**r * (1-q)**(k-r)
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

# ####################  X ~ D/k + Exp(mu), (k, \Delta, c)  ####################### #
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

# ####################  X ~ Pareto(loc, a), (k, \Delta, c)  ####################### #
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

def E_C_pareto_k_c_approx(l, a, d, k, c, w_cancel=True):
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

# ##############################  E[T^2], E[C^2]  ############################## #
# ### X ~ SExp ### #
def E_exp_X_i_j(mu, n, i, j):
  if i > j:
    _j = j
    j = i
    i = _j
  return (H_2(n) - H_2(n-i) + (H(n) - H(n-i))*(H(n) - H(n-j)) )/mu**2

def E_T_2_sexp_k_c(D, mu, k, c):
  return (D/k + H(k)/(c+1)/mu)**2 + H_2(k)/(c+1)**2/mu**2

def E_C_2_sexp_k_c(D, mu, k, c):
  # return ((c+1)*D + k/mu)**2 + (k/mu)**2
  mu_ = (c+1)*mu
  E_C_2 = D**2 + 2*D*k/mu_
  for i in range(1, k+1):
    for j in range(1, k+1):
      E_C_2 += E_exp_X_i_j(mu_, k, i, j)
  
  return (c+1)**2 * E_C_2

def E_T_2_sexp_k_n(D, mu, k, n):
  return (H_2(n) - H_2(n-k))/mu**2 + (D/k + (H(n) - H(n-k))/mu)**2

def E_C_2_sexp_k_n(D, mu, k, n):
  E_C_2 = (n*D/k)**2 + 2*n*D/mu + (n-k)**2*E_exp_X_i_j(mu, n, k, k)
  for i in range(1, k+1):
    E_C_2 += 2*(n-k)*E_exp_X_i_j(mu, n, i, k)
  for i in range(1, k+1):
    for j in range(1, k+1):
      E_C_2 += E_exp_X_i_j(mu, n, i, j)
  
  return E_C_2

# ### X ~ Pareto ### #
def E_pareto_X_i_j(loc, a, n, i, j):
  if i > j:
    _j = j
    j = i
    i = _j
  if a <= max(2/(n-i+1), 1/(n-j+1) ):
    return 0 # None
  return loc**2*G(n+1)/G(n+1-2/a) * G(n-i+1-2/a)/G(n-i+1-1/a) * G(n-j+1-1/a)/G(n-j+1)

def E_T_2_pareto_k_c(loc, a, k, c):
  a_ = (c+1)*a
  if a_ > 1:
    return E_pareto_X_i_j(loc, a_, k, k, k)
  else:
    return None

def E_C_2_pareto_k_c(loc, a, k, c):
  a_ = (c+1)*a
  # if a_ > 2:
  #   return (k*(c+1))**2 * loc**2*a_/(a_-2)
  # else:
  #   None
  E_C_2 = 0
  for i in range(1, k+1):
    for j in range(1, k+1):
      E_C_2 += E_pareto_X_i_j(loc, a_, k, i, j)
  
  return (c+1)**2 * E_C_2

def E_T_2_pareto_k_n(loc, a, k, n):
  return E_pareto_X_i_j(loc, a, n, k, k)

def E_C_2_pareto_k_n(loc, a, k, n):
  E_C_2 = (n-k)**2*E_pareto_X_i_j(loc, a, n, k, k)
  for i in range(1, k+1):
    E_C_2 += 2*(n-k)*E_pareto_X_i_j(loc, a, n, i, k)
  for i in range(1, k+1):
    for j in range(1, k+1):
      E_C_2 += E_pareto_X_i_j(loc, a, n, i, j)
  
  return E_C_2

# ### Wrappers ### #
def E_T_2_k_c(task_t, task_dist_m, k, c, load_m=None):
  if task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    return E_T_2_sexp_k_c(D, mu, k, c)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    return E_T_2_pareto_k_c(loc, a, k, c)

def E_T_2_k_n(task_t, task_dist_m, k, n, load_m=None):
  if task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    return E_T_2_sexp_k_n(D, mu, k, n)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    return E_T_2_pareto_k_n(loc, a, k, n)

def E_C_2_k_c(task_t, task_dist_m, k, c, load_m=None):
  if task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    return E_C_2_sexp_k_c(D, mu, k, c)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    return E_C_2_pareto_k_c(loc, a, k, c)

def E_C_2_k_n(task_t, task_dist_m, k, n, load_m=None):
  if task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    return E_C_2_sexp_k_n(D, mu, k, n)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    return E_C_2_pareto_k_n(loc, a, k, n)

# ***************************  X ~ Pareto, (k, n/c, \Delta) with Relaunch  *************************** #
def Pr_T_g_t_pareto_k_wrelaunch(loc, a, d, k, t):
  q_t = 1 - (loc/t)**a if t > loc else 0
  if d == 0:
    return 1 - q_t**k
  
  q_d = 1 - (loc/d)**a if d > loc else 0
  q1 = 1 - (d/t)**a if t > d else 0
  q2 = 1 - (loc/(t-d))**a if t-d > loc else 0
  # p_T_geq_t = 1 - q_t**k + (t > d)*(q_d + q1*(1-q_d))**k - (t-d > loc)*(q_d + q2*(1-q_d))**k - q_d**k*((t > d) - (t-d > loc))
  sum_ = 0
  for r in range(k):
  # for r in range(k+1):
    sum_ += (q1**(k-r) - q2**(k-r)) * binom(k, r) * q_d**r * (1-q_d)**(k-r)
  p_T_geq_t = 1 - q_t**k + sum_
  
  return p_T_geq_t

def Pr_T_g_t_pareto_k_wrelaunch_approx(loc, a, d, k, t):
  q_t = 1 - (loc/t)**a if t > loc else 0
  if d == 0:
    return 1 - q_t**k
  
  q_d = 1 - (loc/d)**a if d > loc else 0
  q1 = 1 - (d/t)**a if t > d else 0
  q2 = 1 - (loc/(t-d))**a if t-d > loc else 0
  # p_T_geq_t = 1 - q_t**k + (t > d)*(q_d + q1*(1-q_d))**k - (t-d > loc)*(q_d + q2*(1-q_d))**k # - q_d**k*((t > d) - (t-d > loc))
  p_T_geq_t = 1 - q_t**k + (q_d + q1*(1-q_d))**k - (q_d + q2*(1-q_d))**k
  # sum_ = 0
  # for r in range(k+1):
  #   sum_ += (q1**(k-r) - q2**(k-r)) * binom(k, r) * q_d**r * (1-q_d)**(k-r)
  # p_T_geq_t = 1 - q_t**k + sum_
  return p_T_geq_t

def Delta_for_min_p_T_g_t_pareto_k_wrelaunch(loc, a, k):
  t = 20*loc
  d_min = None
  p_T_g_t_min = float('Inf')
  for d in numpy.linspace(0, 2*t, 1000):
    p_T_g_t = Pr_T_g_t_pareto_k_wrelaunch(loc, a, d, k, t)
    if p_T_g_t < p_T_g_t_min:
      p_T_g_t_min = p_T_g_t
      d_min = d
  return d_min

def E_X_n_k_pareto(loc, a, n, k):
  if k == 0:
    return 0
  elif n == k and n > 170:
    return loc*(k+1)**(1/a) * G(1-1/a)
  elif n > 170:
    return loc*((n+1)/(n-k+1))**(1/a)
  return loc*G(n+1)/G(n-k+1)*G(n-k+1-1/a)/G(n+1-1/a)

def E_T_pareto_k_n(loc, a, d, k, n):
  if d == 0 or n == k:
    return E_X_n_k_pareto(loc, a, n, k)
  else:
    log(ERROR, "Cannot be expressed analytically!")
    return None

def Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, k):
  return loc*math.sqrt(G(k+1)*G(1-1/a)/G(k+1-1/a) )

def E_T_pareto_k_n_wrelaunch(loc, a, d, k, n):
  if d == 0:
    return E_X_n_k_pareto(loc, a, n, k)
  q = (d > loc)*(1 - (loc/d)**a)
  
  if d <= loc:
    return d + E_X_n_k_pareto(loc, a, n, k)
  elif n == k:
    def g(k, a):
      if k > 170:
        return loc*(k+1)**(1/a) * G(1-1/a)
      return loc*G(1-1/a)*G(k+1)/G(k+1-1/a)
    return d*(1-q**k) + g(k, a)*(1 + (loc/d-1)*I(1-q,1-1/a,k) )
  elif n > k:
    if d <= loc:
      return d + E_T_pareto_k_n(loc, a, 0, k, n)
    else:
      sum_ = 0
      for r in range(k):
        sum_ += (E_X_n_k_pareto(loc, a, n-r, k-r) - E_X_n_k_pareto(loc, a, k-r, k-r) ) * binom(k,r) * q**r * (1-q)**(k-r)
      return sum_ + E_T_pareto_k_n_wrelaunch(loc, a, d, k, n=k)
      # return d*(1-q**k) + loc*(B(n-k*q+1, -1/a)/B(n-k+1, -1/a) + k*B(k, 1-1/a, u_l=q) - q**k)

def E_T_pareto_k_n_wrelaunch_approx(loc, a, d, k, n):
  q = (d > loc)*(1 - (loc/d)**a)
  
  if n == k:
    return E_T_pareto_k_n_wrelaunch(loc, a, d, k, n)
  elif n > k:
    if d <= loc:
      return d + E_T_pareto_k_n(loc, a, 0, k, n)
    else:
      return d*(1-q**k) + loc*(B(n-k*q+1, -1/a)/B(n-k+1, -1/a) + k*B(k, 1-1/a, u_l=q) - q**k)

def E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel=True):
  if d == 0 and w_cancel:
    if n > 170:
      return loc/(a-1) * (a*n - (n-k)*((n+1)/(n-k+1))**(1/a) )
    return loc*n/(a-1) * (a - G(n)/G(n-k)*G(n-k+1-1/a)/G(n+1-1/a) )
  
  q = (d > loc)*(1 - (loc/d)**a)
  if w_cancel:
    if d <= loc:
      return k*d + E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=True)
    else:
      # EX_given_X_leq_d = a/(a-1)/q * (loc - d*(1-q) )
      # E = 0
      # for r in range(k+1):
      #   E += G(n-r+1)/G(n-r+1-1/a) * binom(k,r) * q**r * (1-q)**(k-r)
      
      # return k*q*EX_given_X_leq_d + k*(1-q)*d + a/(a-1)*loc*(n - k*q) \
      #       - loc/(a-1) * G(n-k+1-1/a)/G(n-k)*E \
      #       - q**k * loc*(n-k)
      
      return a/(a-1)*(k*(1-q)*(loc-d) + n*loc) + k*(1-q)*d - loc*(n-k)*q**k \
            - loc/(a-1) * (n-k)*B(n-k*q+1, -1/a)/B(n-k+1, -1/a) # G(n-k+1-1/a)/G(n-k)*E
  else:
    if d <= loc:
      return k*d + n*loc/(1-1/a)
    else:
      return a/(a-1)*(k*loc*(1-q+q**k) + n*loc*(1-q**k) ) - k*d*(1-q)/(a-1)

def E_C_pareto_k_n_wrelaunch_approx(loc, a, d, k, n, w_cancel=True):
  q = (d > loc)*(1 - (loc/d)**a)
  if w_cancel:
    if d <= loc:
      return k*d + E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=True)
    else:
      return a/(a-1)*(k*(1-q)*(loc-d) + n*loc) + k*(1-q)*d - loc*(n-k)*q**k \
              - loc/(a-1) * (n-k)*B(n-k*q+1, -1/a)/B(n-k+1, -1/a)

def E_T_pareto_k_c_wrelaunch(loc, a, d, k, c):
  if c == 0:
    return E_T_pareto_k_n_wrelaunch(loc, a, d, k, n=k)
  
  ## Following function defs are required to get the expected result
  # def I(u_l, m, n):
  #   return scipy.special.betainc(m, n, u_l)
  # def B(m, n, u_l=1):
  #   if u_l == 1:
  #     return scipy.special.beta(m, n)
  #   else:
  #     return I(u_l, m, n)*B(m, n)
  # def G(z):
  #   return scipy.special.gamma(z)
  
  a_ = (c+1)*a # 1/(c+1)/a
  if d <= loc:
    return d + loc*G(k+1)*G(1-1/a_)/G(k+1-1/a_)
  
  q = (d > loc)*(1 - (loc/d)**a)
  if d == 0:
    return d + loc*G(k+1)*G(1-a_)/G(k+1-a_)
  else:
    f = lambda a: loc*G(1-1/a)/G(-1/a)*B(k-k*q+1, -1/a)
    # print("E_T_pareto_k_c_wrelaunch(loc, a, d, k, 0)= {}".format(E_T_pareto_k_c_wrelaunch(loc, a, d, k, 0) ) )
    return f(a_) - f(a) \
           + E_T_pareto_k_c_wrelaunch(loc, a, d, k, 0)

def E_C_pareto_k_c_wrelaunch(loc, a, d, k, c, w_cancel=True):
  q = (d > loc)*(1 - (loc/d)**a) if d else 0
  if not w_cancel:
    if d <= loc:
      return k*d + k*loc*(c+1)*a/(a-1)
    else:
      return k*(loc - d*(1-q) )*a/(a-1) + k*d*(1-q) \
             + k*loc*(c+1)*(1-q)*a/(a-1)
  else: # w_cancel
    if d <= loc:
      return k*d + k*loc*(c+1)*(c+1)*a/((c+1)*a - 1)
    else:
      return k*(loc - d*(1-q) )*a/(a-1) + k*d*(1-q) \
             + k*loc*(c+1)*(1-q)*(c+1)*a/((c+1)*a-1)

# ****************  Launch all at the beginning, relaunch incomplete at \Delta  ****************** #
def Delta_for_min_E_T_pareto_k_n0_wrelaunch(loc, a, k, n):
  # return loc*math.sqrt(G(k+1)*G(1-1/a)/G(k+1-1/a) )
  return math.sqrt(loc*E_T_pareto_k_n(loc, a, 0, k, n) )

def E_T_pareto_k_cd0_wrelaunch(loc, a, d, k, c):
  q = (d > loc)*(1 - (loc/d)**((c+1)*a) )
  if d <= loc:
    return d + E_T_pareto_k_c(loc, a, 0, k, c)
  else:
    return d*(1 - q**k) + \
      E_T_pareto_k_c(loc, a, 0, k, c)*(1 + (loc/d - 1)*I(1-q,1-1/a,k) )

def E_T_pareto_k_nd0_wrelaunch(loc, a, d, k, n):
  if d == 0:
    return E_X_n_k_pareto(loc, a, n, k)
  q = (d > loc)*(1 - (loc/d)**a)
  
  if d <= loc:
    return d + E_X_n_k_pareto(loc, a, n, k)
  elif n == k:
    def g(k, a):
      if k > 170:
        return loc*(k+1)**(1/a) * G(1-1/a)
      return loc*G(1-1/a)*G(k+1)/G(k+1-1/a)
    
    return d*(1-q**k) + g(k, a)*((loc/d-1)*I(1-q,1-1/a,k) + 1)
  elif n > k:
    if d <= loc:
      return d + E_T_pareto_k_n(loc, a, 0, k, n)
    else:
      # s = 0
      # for r in range(k):
      #   s += (d + E_X_n_k_pareto(loc, a, n-r, k-r) - E_X_n_k_pareto(d, a, n-r, k-r) ) * binom(n,r) * q**r * (1-q)**(n-r)
      # return s + E_X_n_k_pareto(loc, a, n, k)
      return d*I(1-q, n-k+1, k) + ((loc/d - 1)*I(1-q, n-k+1-1/a, k) + 1)*E_T_pareto_k_n(loc, a, 0, k, n)

def E_C_pareto_k_nd0_wrelaunch(loc, a, d, k, n, w_cancel=True):
  if d == 0 and w_cancel:
    if n > 170:
      return loc/(a-1) * (a*n - (n-k)*((n+1)/(n-k+1))**(1/a) )
    return loc*n/(a-1) * (a - G(n)/G(n-k)*G(n-k+1-1/a)/G(n+1-1/a) )
  
  q = (d > loc)*(1 - (loc/d)**a)
  if w_cancel:
    if d <= loc:
      return n*d + E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=True)
    else:
      def p_R_r(r): return binom(n, r) * q**r * (1-q)**(n-r)
      def inner_sum(r):
        s = 0
        for i in range(1, k-r+1):
          s += E_X_n_k_pareto(loc, a, n-r, i) - E_X_n_k_pareto(d, a, n-r, i)
        return s
      
      S = 0
      for r in range(k):
        S += ((n-r)*d + inner_sum(r) + (n-k)*(E_X_n_k_pareto(loc, a, n-r, k-r) - E_X_n_k_pareto(d, a, n-r, k-r) ) ) * p_R_r(r)
      return S + E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n, w_cancel=True)

# **********************  Retain k task at t=\Delta  ********************* #
def E_T_pareto_k_n_retainl_atd(loc, a, k, n, d):
  q = 1 - (loc/d)**a if d > loc else 0
  
  if d <= loc:
    # return E_X_n_k_pareto(loc, a, k, k)
    return None
  else:
    E_T = E_X_n_k_pareto(loc, a, n, k)
    
    sum_ = 0
    # for r in range(k):
    #   sum_ += G(r+1-1/a)/G(r+1) * binom(k, r) * q**r * (1-q)**(k-r)
    # E_T += d*G(k+1)/G(k+1-1/a) * sum_
    for r in range(k):
      l = min(n-r, k)
      # l = n - r # min(n-r, k)
      sum_ += E_X_n_k_pareto(d, a, l, k-r) * binom(k, r) * q**r * (1-q)**(k-r)
    E_T += sum_
    
    sum_ = 0
    # for r in range(k):
    #   sum_ += G(n-r+1)/G(n-r+1-1/a) * binom(k, r) * q**r * (1-q)**(k-r)
    # E_T -= d*G(n-k+1-1/a)/G(n-k+1) * sum_
    for r in range(k):
      sum_ += E_X_n_k_pareto(d, a, n-r, k-r) * binom(k, r) * q**r * (1-q)**(k-r)
    E_T -= sum_
    
    return E_T

def approx_E_T_pareto_k_n_retainl_atd(loc, a, k, n, d):
  q = 1 - (loc/d)**a if d > loc else 0
  return E_X_n_k_pareto(loc, a, n, k) + \
         d*B(k*q+1-1/a, 1/a)/B(k+1-1/a, 1/a) - d*B(n-k*q+1, -1/a)/B(n-k+1, -1/a)

def E_C_pareto_k_n_retainl_atd(loc, a, k, n, d):
  if d <= loc:
    # return (n-k)*d + k*E_X_n_k_pareto(loc, a, k, k)
    return None
  else:
    q = 1 - (loc/d)**a if d > loc else 0
    # E_X__X_leq_d = loc*a/(a-1)*(1 - (loc/d)**(a-1) )/(1 - (loc/d)**a)
    # X = Pareto(loc, a)
    # E_X__X_leq_d = mpmath.quad(lambda x: x*X.pdf(x), [0, d] )/X.cdf(d)
    
    E_C = 0
    for r in range(k):
      l = min(n-r, k)
      # l = n - r
      # l = k - r
      
      sum_ = d*(n-r-l)
      for i in range(1, k-r+1):
        sum_ += E_X_n_k_pareto(d, a, l, i) - E_X_n_k_pareto(d, a, n-r, i)
      sum_ += (l-k+r)*E_X_n_k_pareto(d, a, l, k-r) - (n-k)*E_X_n_k_pareto(d, a, n-r, k-r)
      
      E_C += sum_ * binom(k, r) * q**r * (1-q)**(k-r)
    
    return E_C + E_C_pareto_k_n_wrelaunch(loc, a, 0, k, n)

# ### X ~ TPareto ### #
def E_TPareto_X_n_i(l, u, a, n, i):
  # return l* G(n+1)/G(i)/G(n-i+1) * mpmath.quad(lambda x: (1 - (1-u**(-a))*x)**(-1/a) * x**(i-1)*(1-x)**(n-i), [0, 1] )
  return l* G(n+1)/G(i)/G(n-i+1) * \
    scipy.integrate.quad(lambda x: (1 - (1-u**(-a))*x)**(-1/a) * x**(i-1)*(1-x)**(n-i), 0, 1)[0]

def E_T_k_n_TPareto(l, u, a, k, n):
  return E_TPareto_X_n_i(l, u, a, n, k)

def E_C_k_n_TPareto(l, u, a, k, n):
  E_C = 0
  for i in range(k+1):
    E_C += E_TPareto_X_n_i(l, u, a, n, i)
  E_C += (n-k)*E_TPareto_X_n_i(l, u, a, n, k)
  return E_C

# ### When redundancy changes the tail ### #
def a_wred(ro_0, a_0, ro):
  # K = ro/(1 - ro)*(1-ro_0)/ro_0
  # return K*a_0/((K-1)*a_0 + 1)
  # return a_0*ro_0/ro
  # return a_0*(1-ro)/(1-ro_0)
  
  K = ro_0 + math.sqrt(a_0 - 1)
  return a_0 - (K*ro-ro_0)**2

def a_wred_(a_0, r):
  # K = (1 + math.sqrt(a_0 - 1))/10
  # return a_0 - (K*r-1)**2
  
  # a_last = 1 # 0.2 # 1
  # r_mid, r_max = 3, 10
  # if r <= r_mid:
  #   return a_0
  # m = (a_0 - a_last)/(r_mid - r_max)
  # n = a_0 - m*r_mid
  # return m*r + n
  
  # r_min, r_mid, r_max = 1, 3, 10
  # a_min, a_mid, a_max = 1, 0.8*a_0, a_0
  r_1, r_2, r_3 = 1, 5, 10 # 1, 5, 10
  a_1, a_2, a_3 = a_0, 0.9*a_0, 1
  return a_1*(r - r_2)/(r_1 - r_2)*(r - r_3)/(r_1 - r_3) + \
         a_2*(r - r_1)/(r_2 - r_1)*(r - r_3)/(r_2 - r_3) + \
         a_3*(r - r_1)/(r_3 - r_1)*(r - r_2)/(r_3 - r_2)

def plot_a_wred():
  ro_0 = 0.1
  a_0 = 10
  
  ro_l, a_l = [], []
  for ro in numpy.linspace(ro_0, 0.95, 25):
    ro_l.append(ro)
    a_l.append(a_wred(ro_0, a_0, ro) )
  plot.plot(ro_l, a_l, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.xlabel(r'$\rho$', fontsize=13)
  plot.ylabel(r'$\alpha$', fontsize=13)
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_a_wred.png")
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_Tnp1_over_Tn():
  ro_0, a_0 = 0.3, 1.2
  def Tnp1_over_Tn(k, n):
    a_n = a_wred(ro_0, a_0, ro_0*n/k)
    a_np1 = a_wred(ro_0, a_0, ro_0*(n+1)/k)
    R = (n+1)/(n-k+1) * G(n+2-1/a_np1-k)/G(n+2-1/a_np1) * G(n+1-1/a_n)/G(n+1-1/a_n-k)
    # return math.log(R) if R > 0 else None
    return R
  
  def approx_Tnp1_over_Tn(k, n):
    a_n = a_wred(ro_0, a_0, ro_0*n/k)
    a_np1 = a_wred(ro_0, a_0, ro_0*(n+1)/k)
    # return (n+1)/(n-k+1) * ((n+1-1/a_n)/(n+2-1/a_np1))**k
    
    # b_1 = n-k+1-1/a_n
    # b_2 = n-k+2-1/a_np1
    # R = (n+1)/(n-k+1) * k**(b_1 - b_2) * G(b_2)/G(b_1)
    
    R = ((n+2)/(n-k+2))**(1/a_np1) * ((n+1)/(n-k+1))**(-1/a_n)
    # return math.log(R) if R > 0 else None
    return R
  
  def plot_lb_a_n_over_a_np1(k):
    def lb_a_n_over_a_np1(k, n):
      return math.log(1 + k/(n-k+1))/math.log(1 + k/(n-k+2))
    
    n_l, lb_a_n_over_a_np1_l = [], []
    for n in numpy.arange(k, 5*k, 1):
      n_l.append(n)
      lb_a_n_over_a_np1_l.append(lb_a_n_over_a_np1(k, n) )
    plot.plot(n_l, lb_a_n_over_a_np1_l, label=r'$k= {}$'.format(k), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    
  # plot_lb_a_n_over_a_np1(k=15)
  # plot_lb_a_n_over_a_np1(k=100)
  # plot.legend()
  # plot.xlabel(r'$n$', fontsize=13)
  # plot.ylabel(r'$\alpha_n/\alpha_{n+1} \geq$', fontsize=13)
  # fig = plot.gcf()
  # fig.tight_layout()
  # plot.savefig("plot_lb_a_n_over_a_np1.png")
  # plot.gcf().clear()
  
  k = 100
  n_l, r_l, r_approx_l = [], [], []
  for n in numpy.arange(k, 2*k, 1):
    n_l.append(n)
    r_l.append(Tnp1_over_Tn(k, n) )
    r_approx_l.append(approx_Tnp1_over_Tn(k, n) )
  plot.plot(n_l, r_l, label='Exact', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(n_l, r_approx_l, label='Approx', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.legend()
  plot.xlabel(r'$n$', fontsize=13)
  plot.ylabel(r'$E[T_{n+1}]/E[T_n]$', fontsize=13)
  plot.title(r'$k= {}$'.format(k) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_Tnp1_over_Tn.png")
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_Trj_over_Tri():
  ro_0, a_0 = 0.3, 1.2
  def Trj_over_Tri(ri, rj):
    task_t = 'Pareto'
    ni = math.floor(k*ri)
    task_dist_m = {'loc': 1, 'a': a_wred(ro_0, a_0, ro_0*ri) }
    Tri = E_T_k_l_n(task_t, task_dist_m, 0, k, k, ni)
    
    nj = math.floor(k*rj)
    task_dist_m = {'loc': 1, 'a': a_wred(ro_0, a_0, ro_0*rj) }
    Trj = E_T_k_l_n(task_t, task_dist_m, 0, k, k, nj)
    return Trj/Tri
  
  def suffcond_ai_over_aj_for_ETgain(k, ri, rj):
    ni = math.floor(k*ri)
    nj = math.floor(k*rj)
    if ni == k or nj == k: return None
    return math.log(ni/(ni-k+1) )/math.log((nj+1)/(nj-k) )
  # def approx_suffcond_ai_over_aj_for_ETgain(k, ri, rj):
  #   return math.log(k*ri/(k*(ri-1)+1) )/math.log((k*rj+1)/(k*(rj-1) ) )
  
  def suffcond_ai_over_aj_for_ETpain(k, ri, rj):
    ni = math.floor(k*ri)
    nj = math.floor(k*rj)
    if ni == k or nj == k: return None
    return math.log((ni+1)/(ni-k) )/math.log(nj/(nj-k+1) )
  
  def approxcond_an_over_anp1_for_ETpain(k, n):
    return math.log(1 + k/(n-k+1) )/math.log(1 + k/(n-k+2) )
  
  # def looser_suffcond_ai_over_aj_for_ETpain(k, ri, rj):
  #   return (k*ri+1)*rj/(ri-1)/(k-1)
  
  def plot_cond_an_over_anp1(k):
    n_l, suff_for_ETgain_l, suff_for_ETpain_l = [], [], []
    approx_for_ETpain_l = []
    for n in numpy.arange(k, 4*k, 1):
      n_l.append(n)
      ri, rj = n/k, (n+1)/k
      suff_for_ETgain_l.append(suffcond_ai_over_aj_for_ETgain(k, ri, rj) )
      suff_for_ETpain_l.append(suffcond_ai_over_aj_for_ETpain(k, ri, rj) )
      approx_for_ETpain_l.append(approxcond_an_over_anp1_for_ETpain(k, n) )
    plot.plot(n_l, suff_for_ETgain_l, label=r'$k= {}$, Gain if {}'.format(k, r'$\alpha_n/\alpha_{n+1} \leq$'), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(n_l, suff_for_ETpain_l, label=r'$k= {}$, Pain if {}'.format(k, r'$\alpha_n/\alpha_{n+1} \geq$'), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(n_l, approx_for_ETpain_l, label=r'$k= {}$, Pain if {}'.format(k, r'$\alpha_n/\alpha_{n+1} \gtrapprox$'), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot_cond_an_over_anp1(k=15)
  plot_cond_an_over_anp1(k=15)
  plot.legend()
  # plot.xscale('log')
  # plot.yscale('log')
  plot.xlabel(r'$n$', fontsize=13)
  plot.ylabel(r'', fontsize=13)
  fig = plot.gcf()
  fig.tight_layout()
  plot.savefig("plot_an_over_anp1.png")
  plot.gcf().clear()
  
  # n_l, r_l, r_approx_l = [], [], []
  # for n in numpy.arange(k, 2*k, 1):
  #   n_l.append(n)
  #   r_l.append(Tnp1_over_Tn(n) )
  #   r_approx_l.append(approx_Tnp1_over_Tn(n) )
  # plot.plot(n_l, r_l, label='Exact', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot.plot(n_l, r_approx_l, label='Approx', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot.legend()
  # plot.xlabel(r'$n$', fontsize=13)
  # plot.ylabel(r'$E[T_{n+1}]/E[T_n]$', fontsize=13)
  # plot.title(r'$k= {}$'.format(k) )
  # fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  # fig.tight_layout()
  # plot.savefig("plot_Tnp1_over_Tn.png")
  # plot.gcf().clear()
  log(WARNING, "done.")

def n_max_before_ETpain(ro_0, a_0, k):
  def approxcond_an_over_anp1_for_ETpain(k, n):
    return math.log(1 + k/(n-k+1) )/math.log(1 + k/(n-k+2) )
  n = k
  a_n = a_wred(ro_0, a_0, ro=1)
  while True:
    # if ro >= 0.99: return None
    a_np1 = a_wred(ro_0, a_0, ro=ro_0*(n+1)/k)
    if a_n/a_np1 >= approxcond_an_over_anp1_for_ETpain(k, n):
      break
    a_n = a_np1
    n += 1
  return n

# ******************************  Wrappers  ****************************** #
def E_T_k_l_n(task_t, task_dist_m, d, k, l, n, load_m=None):
  if task_t == "Exp":
    mu = task_dist_m["mu"]
    if l == k: return E_T_exp_k_n(mu, d, k, n)
    else: return E_T_exp_k_l_n(mu, d, k, l, n)
  elif task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    if l == k: return E_T_shiftedexp_k_n(D, mu, d, k, n)
    else: return E_T_shiftedexp_k_l_n(D, mu, d, k, l, n)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    if load_m is not None:
      # ro_0, a_0 = load_m['ro_0'], load_m['a_0']
      # ro = ro_0*n/k
      # if ro >= 0.99: return None
      # a = a_wred(ro_0, a_0, ro)
      r = n/k
      a = a_wred_(load_m['a_0'], r)
      log(WARNING, "r= {}, a= {}".format(r, a) )
    return E_T_pareto_k_n(loc, a, d, k, n)
  elif task_t == "TPareto":
    l, u, a = task_dist_m["l"], task_dist_m["u"], task_dist_m["a"]
    return E_T_k_n_TPareto(l, u, a, k, n)

def E_C_k_l_n(task_t, task_dist_m, d, k, l, n, w_cancel, load_m=None):
  if task_t == "Exp":
    mu = task_dist_m["mu"]
    if l == k: return E_C_exp_k_n(mu, d, k, l, n, w_cancel)
    else: return E_C_exp_k_l_n(mu, d, k, l, n, w_cancel)
  elif task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    return E_C_shiftedexp_k_l_n(D, mu, d, k, l, n, w_cancel)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    if load_m is not None:
      # ro_0, a_0 = load_m['ro_0'], load_m['a_0']
      # ro = ro_0*n/k
      # if ro >= 0.99: return None
      # a = a_wred(ro_0, a_0, ro)
      r = n/k
      a = a_wred_(load_m['a_0'], r)
      log(WARNING, "r= {}, a= {}".format(r, a) )
    return E_C_pareto_k_n_wrelaunch(loc, a, d, k, n, w_cancel=w_cancel)
  elif task_t == "TPareto":
    l, u, a = task_dist_m["l"], task_dist_m["u"], task_dist_m["a"]
    return E_C_k_n_TPareto(l, u, a, k, n)

def E_T_k_c(task_t, task_dist_m, d, k, c, load_m=None):
  if task_t == "Exp":
    mu = task_dist_m["mu"]
    return E_T_exp_k_c(mu, d, k, c)
  elif task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    return E_T_shiftedexp_k_c(D, mu, d, k, c)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    if load_m is not None:
      # ro_0, a_0 = load_m['ro_0'], load_m['a_0']
      # ro = ro_0*(c+1)
      # if ro >= 0.99: return None
      # a = a_wred(ro_0, a_0, ro)
      r = c+1
      a = a_wred_(load_m['a_0'], r)
      log(WARNING, "r= {}, a= {}".format(r, a) )
    return E_T_pareto_k_c(loc, a, d, k, c)

def E_C_k_c(task_t, task_dist_m, d, k, c, w_cancel, load_m=None, approx=False):
  if task_t == "Exp":
    mu = task_dist_m["mu"]
    return E_C_exp_k_c(mu, d, k, c, w_cancel)
  elif task_t == "SExp":
    D, mu = task_dist_m["D"], task_dist_m["mu"]
    return E_C_shiftedexp_k_c(D, mu, d, k, c, w_cancel)
  elif task_t == "Pareto":
    loc, a = task_dist_m["loc"], task_dist_m["a"]
    if load_m is not None:
      # ro_0, a_0 = load_m['ro_0'], load_m['a_0']
      # ro = ro_0*(c+1)
      # if ro >= 0.99: return None
      # a = a_wred(ro_0, a_0, ro)
      r = c+1
      a = a_wred_(load_m['a_0'], r)
      log(WARNING, "r= {}, a= {}".format(r, a) )
    if approx:
      return E_C_pareto_k_c_approx(loc, a, d, k, c, w_cancel)
    return E_C_pareto_k_c(loc, a, d, k, c, w_cancel)

def plot_deneme():
  # # To check if sum of ratios of gamma formula works for non-integer sum indices (Works)
  # def compare(k, a):
  #   actual = 0
  #   for i in range(2, k+1):
  #     actual += G(i-2/a)/G(i-1/a)
  #   approx = G(k+1-2/a)/(1-1/a)/G(k-1/a) - G(2-2/a)/(1-1/a)/G(1-1/a)
  #   print("actual= {}, approx= {}".format(actual, approx) )
  
  # compare(k=10, a=2)
  # compare(k=10, a=3)
  # compare(k=10, a=4)
  
  k = 100
  loc = 3
  a_l, d_l = [], []
  for a in numpy.linspace(1, 40, 100):
    a_l.append(a)
    d_l.append(Delta_for_min_E_T_pareto_k_wrelaunch(loc, a, k) )
  plot.plot(a_l, d_l, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.xlabel(r'$\alpha$', fontsize=12)
  plot.ylabel(r'$\Delta^*$', fontsize=12)
  plot.savefig("plot_deneme.png")
  log(WARNING, "done.")

def deneme():
  p = 0.99
  def prob_fast(k, n):
    return sum([binom(n, i) * p**i * (1-p)**(n-i) for i in range(k, n+1) ] )
  k = 100
  for n in range(k, 2*k):
    pf = prob_fast(k, n)
    print("n= {}, pf= {}".format(n, pf) )

if __name__ == "__main__":
  # plot_Pr_T_g_t_G_1red()
  # plot_pareto_zerodelay_red()
  # plot_deneme()
  
  # plot_compare_tails()
  deneme()
  
  # plot_a_wred()
  # plot_Tnp1_over_Tn()
  # plot_Trj_over_Tri()
