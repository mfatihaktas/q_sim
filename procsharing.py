import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

import sys, pprint, math, numpy, scipy, simpy, sympy, getopt, itertools

from sim import *
from arepeat_models import *
from arepeat_sim import *

class Proc(object):
  def __init__(self, _id, size, remaining):
    self._id = _id
    self.size = size
    self.remaining = remaining
    
    self.entrance_time = None
  
  def __repr__(self):
    return "Proc[id= {}, size= {}, remaining= {}]".format(self._id, self.size, self.remaining)
  
  def deep_copy(self):
    p = Proc(self._id, self.size, self.remaining)
    p.entrance_time = self.entrance_time
    return p

# *******************************  PacketGenerator  ****************************** #
class PG(object): # Packet Generator
  def __init__(self, env, ar, psize_dist):
    self.env = env
    self.ar = ar
    self.psize_dist = psize_dist
    
    self.n_sent = 0
    self.out = None
  
  def init(self):
    self.env.process(self.run() )
  
  def run(self):
    while 1:
      yield self.env.timeout(random.expovariate(self.ar) )
      self.n_sent += 1
      
      s = self.psize_dist.gen_sample()
      self.out.put(Proc(self.n_sent, s, s) )

class PSQ(object): # Process Sharing Queue
  def __init__(self, env, h):
    self.env = env
    self.h = h
    
    self.p_l = []
    self.sinterrupt_flag = False
    self.sinterrupt = None
    self.got_busy = None
    
    self.restime_l = []
    self.slowdown_l = []
    
    self.store = simpy.Store(env)
    self.action = env.process(self.serv_run() )
    self.action = env.process(self.put_run() )
  
  def __repr__(self):
    return "PSQ[h= {}]".format(self.h)
  
  def serv_run(self):
    while True:
      p_l = self.p_l[:self.h]
      if len(p_l) == 0:
        sim_log(DEBUG, self.env, self, "idle; waiting for arrival", None)
        self.got_busy = self.env.event()
        yield (self.got_busy)
        sim_log(DEBUG, self.env, self, "got busy!", None)
        continue
      r_l = [p.remaining for p in p_l]
      t = min(r_l)
      i_min = r_l.index(t)
      serv_size = len(p_l)
      sim_log(DEBUG, self.env, self, "back to serv; t= {}, i_min= {}, serv_size= {}".format(t, i_min, serv_size), None)
      start_t = self.env.now
      
      self.sinterrupt = self.env.event()
      yield (self.sinterrupt | self.env.timeout(t) )
      serv_t = (self.env.now - start_t)/serv_size
      for i in range(serv_size):
        try:
          self.p_l[i].remaining -= serv_t
        except IndexError:
          break
      
      if self.sinterrupt_flag:
        sim_log(DEBUG, self.env, self, "serv interrupted", None)
        self.sinterrupt_flag = False
        self.sinterrupt = None
      else:
        p = self.p_l.pop(i_min)
        sim_log(DEBUG, self.env, self, "serv done", p)
        
        lifetime = self.env.now - p.entrance_time
        self.restime_l.append(lifetime)
        self.slowdown_l.append(lifetime/p.size)
  
  def put_run(self):
    while True:
      p = (yield self.store.get() )
      _l = len(self.p_l)
      self.p_l.append(p)
      if _l == 0:
        self.got_busy.succeed()
      elif _l < self.h:
        self.sinterrupt_flag = True
        self.sinterrupt.succeed()
  
  def put(self, p, preempt=False):
    p.entrance_time = self.env.now
    sim_log(DEBUG, self.env, self, "recved", p)
    
    return self.store.put(p.deep_copy() )

class QMonitor(object):
  def __init__(self, env, q, poll_interval):
    self.q = q
    self.env = env
    self.poll_interval = poll_interval
    
    # self.pollt_l = []
    # self.qlength_l = []
    self.qbusy_l = []
    self.action = env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_interval)
      
      self.qbusy_l.append(len(q.p_l) != 0)
      # self.pollt_l.append(self.env.now)
      # self.qlength_l.append(self.q.length() )

# *********************************  Sim  *********************************** #
def sim_psq(num_f_run, ar, h, psize_dist):
  sum_ = 0
  for f in range(num_f_run):
    log(WARNING, "ar= {}".format(ar) )
    
    env = simpy.Environment()
    pg = PG(env, ar, psize_dist)
    psq = PSQ(env, h)
    pg.out = psq
    pg.init()
    env.run(until=50000)
    
    l = psq.slowdown_l
    if len(l): sum_ += float(sum(l) )/len(l)
  E_sl = sum_/num_f_run
  print(">> E_sl= {}".format(E_sl) )
  if E_sl > 50: return None
  return E_sl

def plot_psq():
  # D, mu = 1, 0.5
  # psize_dist = Exp(mu, D=D)
  # proc_in_latex = "Exp(D={}, \mu={})".format(D, mu)
  
  loc, a = 1, 4
  psize_dist = Pareto(loc, a)
  proc_in_latex = r'Pareto(\lambda={}, \alpha={})'.format(loc, a)
  
  num_f_run = 1
  def plot_(h):
    ar_l, E_sl_l, E_T_l = [], [], []
    ar_ub = 1/psize_dist.mean() # /h
    log(WARNING, "h= {}, ar_ub= {}, psize_dist= {}".format(h, ar_ub, psize_dist) )
    
    # for ar in numpy.linspace(0.05, ar_ub, 5):
    for ar in numpy.linspace(0.05, 1.5*ar_ub, 10):
      ar_l.append(ar)
    
      E_sl = sim_psq(num_f_run, ar, h, psize_dist)
      E_sl_l.append(E_sl)
    print("E_sl_l= \n{}".format(pprint.pformat(E_sl_l) ) )
    plot.plot(ar_l, E_sl_l, label="h= {}".format(h), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot_(h=1)
  plot_(h=2)
  # plot_(h=3)
  plot_(h=4)
  plot_(h=1000)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'$\lambda$', fontsize=12)
  plot.ylabel(r'$E[Slowdown]$', fontsize=12)
  plot.title(r'$P \sim {}$'.format(proc_in_latex) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_psq.png")
  log(WARNING, "done; psize_dist= {}".format(psize_dist) )

# ****************************************  Fitting the tail  ************************************ #
def fit_pareto(sample_l):
  n = len(sample_l)
  
  fit_upper_tail = False # True
  if not fit_upper_tail:
    l = sample_l[-1]
    D = 0
    for s in sample_l:
      D += math.log(s) - math.log(l)
    a = (n-1)/D
  elif fit_upper_tail:
    l = sample_l[-1]
    i = int(math.sqrt(n) ) # int(n*0.3)
    sample_l = sample_l[:i]
    l_ = sample_l[-1]
    D = 0
    for s in sample_l:
      D += math.log(s) - math.log(l_)
    a = i/D
  log(WARNING, "done; l= {}, a= {}".format(l, a) )
  return l, a

def fit_tpareto(sample_l):
  # sample_l is ordered in descending order
  n = len(sample_l)
  log(WARNING, "n= {}".format(n) )
  fit_upper_tail = False # True
  def solve_a(eq):
    a = 0.01
    _a = None
    while True:
      if eq(a) > 0:
        _a = a
        a += 0.01
      else:
        return _a if _a is not None else 0.01
  
  u = sample_l[0]
  if not fit_upper_tail:
    l = sample_l[-1]
    r = l/u
    # Did not work somehow
    # a = sympy.Symbol('a')
    # a = sympy.solve(n/a + n*r**a*math.log(r)/(1-r**a) - sum([math.log(x/l) for x in sample_l] ) )
    a = solve_a(lambda a: n/a + n*r**a*math.log(r)/(1-r**a) - sum([math.log(x/l) for x in sample_l] ) )
  else:
    i = int(math.sqrt(n) ) # int(n*0.3)
    X_ip1 = sample_l[i+1]
    r = X_ip1/u
    a = solve_a(lambda a: i/a + i*r**a*math.log(r)/(1-r**a) - sum([math.log(x) - math.log(X_ip1) for x in sample_l[:i+1] ] ) )
    l = i**(1/a) * X_ip1*(n - (n-i)*(X_ip1/u)**a)**(-1/a)
  log(WARNING, "done; l= {}, u= {}, a= {}".format(l, u, a) )
  return l, u, a

def plot_psq_tail():
  # D, mu = 1, 0.5
  # psize_dist = Exp(mu, D=D)
  # proc_in_latex = "Exp(D={}, \mu={})".format(D, mu)
  
  # l, a = 1, 2 # 4
  # psize_dist = Pareto(l, a)
  # proc_in_latex = r'Pareto(\lambda={}, \alpha={})'.format(l, a)
  
  l, u, a = 1, 10**10, 1.1
  psize_dist = TPareto(l, u, a)
  proc_in_latex = r'TPareto(l={}, u={}, \alpha={})'.format(l, u, a)
  
  # l, u, p = 1, 10, 0.2
  # psize_dist = Bern(l, u, p)
  # proc_in_latex = r'Bern(l={}, u={}, p={})'.format(l, u, p)
  
  # l, u = 1, 100
  # psize_dist = DUniform(l, u)
  # proc_in_latex = r'DUniform(l={}, u={})'.format(l, u)
  log(WARNING, "psize_dist= {}".format(psize_dist) )
  
  ar_ub = 1/psize_dist.mean()
  
  ar_l, ro_l = []
  tpar_l_l, tpar_u_l, tpar_a_l = [], [], []
  par_l_l, par_a_l = [], []
  def plot_(num_frun, ar, h):
    log(WARNING, "ar= {}, h= {}".format(ar, h) )
    
    ro_sum = 0
    tpar_l_sum, tpar_u_sum, tpar_a_sum = 0, 0, 0
    par_l_sum, par_a_sum = 0, 0
    for f in range(num_frun):
      env = simpy.Environment()
      pg = PG(env, ar, psize_dist)
      psq = PSQ(env, h)
      pg.out = psq
      pg.init()
      qm = QMonitor(env, psq, poll_interval=0.1)
      env.run(until=50000*20)
      
      sl_l = numpy.sort(psq.slowdown_l)
      E_Sl = float(sum(sl_l) )/len(sl_l)
      print(">>> E_Sl= {}".format(E_Sl) )
      if E_Sl > 1000*10:
        return None
      x_l = sl_l[::-1]
      # y_sim_l = numpy.arange(sl_l.size)/sl_l.size
      # plot.plot(x_l, y_sim_l, label="ar= {}".format(ar), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
      i_ = None
      for i in range(len(x_l)-1, 0, -1):
        if x_l[i] > 1.01: i_ = i; break
      x_l = x_l[:i_]
      # y_l = numpy.arange(x_l.size)/x_l.size
      # plot.plot(x_l, y_l, label="ar= {}".format(ar), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
      
      ro_sum += sum(qm.busy_l)/len(qm.busy_l)
      l, u, a = fit_tpareto(x_l)
      tpar_l_sum += l
      tpar_u_sum += u
      tpar_a_sum += a
      l, a = fit_pareto(x_l)
      par_l_sum += l
      par_u_sum += u
    ar_l.append(ar)
    ro_l.append(ro_sum/num_frun)
    
    l, u, a = tpar_l_sum/num_frun, tpar_u_sum/num_frun, tpar_a_sum/num_frun
    print("l= {}, u= {}, a= {}".format(l, u, a) )
    tpar_l_l.append(l)
    tpar_u_l.append(u)
    tpar_a_l.append(a)
    # rv = TPareto(l, u, a)
    # y_l = []
    # for x in x_l: y_l.append(rv.tail(x) )
    # plot.plot(x_l, y_l, label=r'$TPareto(l= %.2f, u= %.2f, \alpha= %.2f), \lambda= %.2f$' % (l, u, a, ar), color=next(dark_color), linestyle='-')
    
    l, a = tpar_l_sum/num_frun, tpar_a_sum/num_frun
    par_l_l.append(l)
    par_a_l.append(a)
    # rv = Pareto(l, a)
    # y_l = []
    # for x in x_l: y_l.append(rv.tail(x) )
    # plot.plot(x_l, y_l, label=r'$Pareto(l= %.2f, \alpha= %.2f), \lambda= %.2f$' % (l, a, ar), color=next(dark_color), linestyle='-')
    
    # plot.legend()
    # plot.xscale('log')
    # plot.yscale('log')
    # plot.xlabel(r'Slowdown', fontsize=13)
    # plot.ylabel(r'Tail distribution', fontsize=13)
    # plot.title(r'$P \sim {}$, $h= {}$, $\lambda= {}$'.format(proc_in_latex, h, ar) )
    # plot.savefig("plot_psq_h_{}_ar_{}.png".format(h, ar) )
    # plot.gcf().clear()
    return 0
  num_frun = 3
  h = 8
  for ar in numpy.arange(0.01, 5*ar_ub, 0.025):
  # for ar in numpy.linspace(0.05, ar_ub, 10):
    if plot_(num_frun, ar, h) is None:
      break
  
  plot.title(r'$P \sim {}, h= {}$'.format(proc_in_latex, h) )
  fig, axes = plot.subplots(2, 1, sharex=True)
  axes[0].plot(ro_l, tpar_a_l, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  axes[0].set_ylabel(r'$\alpha$', fontsize=13)
  axes[1].set_yscale('log')
  axes[1].plot(ro_l, tpar_u_l, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  axes[1].set_ylabel(r'$u$', fontsize=13)
  axes[1].set_xlabel(r'$\rho', fontsize=13)
  plot.savefig("plot_fitted_tpar_h_{}.png".format(h) )
  plot.gcf().clear()
  
  plot.title(r'$P \sim {}, h= {}$'.format(proc_in_latex, h) )
  plot.plot(ro_l, par_a_l, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.xlabel(r'$\rho$', fontsize=13)
  plot.ylabel(r'$\alpha$', fontsize=13)
  plot.savefig("plot_fitted_par_h_{}.png".format(h) )
  plot.gcf().clear()
  log(WARNING, "done.")

def MG1_T():
  l, u, a = 1, 100, 1.5
  B = TPareto(l, u, a)
  print("B.mean= {}".format(B.mean() ) )
  
  def Lt_B(s):
    # Pareto(l, a)
    # return a*(l*s)**a * math.exp(l*s) * G(-a)*scipy.special.gammaincc(-a, l*s)
    
    K = a*l**a/(1 - (l/u)**a)
    # return K*(l**(-a)*scipy.special.expn(a+1, s*l) - u**(-a)*scipy.special.expn(a+1, s*u) )
    # return 1/(s+1)
    
    # def complex_quad(func, a, b, **kwargs):
    #   func_real = lambda x: scipy.real(func(x))
    #   real_integral = scipy.integrate.quad(func_real, a, b, **kwargs)
    #   imag_integral = scipy.integrate.quad(lambda x: scipy.imag(func(x)), a, b, **kwargs)
    #   # return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])
    #   return real_integral[0] + 1j*imag_integral[0]
    # return K*complex_quad(lambda x: x**(-a-1)*scipy.exp(-s*x), l, u)
    
    return K*mpmath.quad(lambda x: x**(-a-1)*math.exp(-s*x), [l, u] )
    # print("s= {}".format(s) )
    # return K*scipy.integrate.quad(lambda x: x**(-a-1)*scipy.exp(-s*x), l, u)
  
  def Lt_T(ar, s):
    ro = ar/B.mean()
    if (ar*Lt_B(s) + s - ar) == 0:
      return None
    return (1-ro)*Lt_B(s)*s/(ar*Lt_B(s) + s - ar)
  
  def chernoff_bound(ar, s, t):
    return Lt_T(ar, -s)*math.exp(-s*t)
  # s_min = scipy.optimize.minimize_scalar(lambda s: chernoff_bound(ar, s, t=100), bounds=(0, 0.6), method='bounded')
  # print("s_min=\n{}".format(s_min) )
  # s_min = float(s_min['x'] )
  # print("s_min= {}".format(s_min) )
  
  def s_for_minchernoff(ar, t):
    _s = 0.001
    _c = chernoff_bound(ar, _s, t)
    while True:
      s = _s + 0.001
      c = chernoff_bound(ar, s, t)
      if c > _c:
        return _s
      _s, _c = s, c
  
  ar = 0.2
  t_l, chernoff_l = [], []
  # for t in numpy.linspace(l, 10000*l, 100):
  for t in numpy.logspace(1, 4, 100):
    # Does not work, don't try again!
    # f_ = mpmath.invertlaplace(Lt_B, t, method='talbot')
    
    s_min = s_for_minchernoff(ar, t)
    print("s_min= {}".format(s_min) )
    cb = chernoff_bound(ar, s_min, t)
    if cb < 0: break
    t_l.append(t)
    print("cb= {}".format(cb) )
    chernoff_l.append(cb)
  plot.plot(t_l, chernoff_l, label=r'Chernoff bound, ar= {}'.format(ar), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.legend()
  plot.xlabel(r'$t$', fontsize=12)
  plot.xscale('log')
  plot.yscale('log')
  
  # s_l, Lt_B_l, Lt_T_l, cb_l = [], [], [], []
  # for s in numpy.linspace(0, 1, 20):
  #   s_l.append(s)
  #   # Lt_B_l.append(Lt_B(-s) )
  #   Lt_T_l.append(Lt_T(ar, -s) )
  #   # cb_l.append(chernoff_bound(ar, s, t=100) )
  # # plot.plot(s_l, Lt_B_l, label=r'$B(s)$', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot.plot(s_l, Lt_T_l, label=r'$T(s)$, $\lambda$= {}'.format(ar), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # # plot.plot(s_l, cb_l, label=r'$Chernoff$, $\lambda$= {}'.format(ar), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot.legend()
  # plot.xlabel(r'$s$', fontsize=12)
  
  plot.title(r'$B \sim {}$'.format(B) )
  plot.savefig("MG1_T.png")
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_EC_vs_ET_wsim():
  l, u, a = 1, 10**10, 1.1 # 1, 100, 1.5
  psize_dist = TPareto(l, u, a)
  proc_in_latex = r'TPareto(l={}, u={}, \alpha={})'.format(l, u, a)
  log(WARNING, "psize_dist= {}".format(psize_dist) )
  
  ar_ub = 1/psize_dist.mean()
  def sim(ar, h):
    env = simpy.Environment()
    pg = PG(env, ar, psize_dist)
    psq = PSQ(env, h)
    pg.out = psq
    pg.init()
    env.run(until=50000*20)
    
    E_Sl = float(sum(psq.slowdown_l) )/len(psq.slowdown_l)
    print(">>> E_Sl= {}".format(E_Sl) )
    if E_Sl > 1000*10: return None # 200 1000*1
    return psq.restime_l # psq.slowdown_l
  
  def plot_EC_vs_ET(num_frun, h, ar, k):
    x_sim_l, y_sim_l = [], []
    x_tpar_l, y_tpar_l = [], []
    x_par_l, y_par_l = [], []
    n_k = 0
    done = False
    while not done:
      n = k + n_k
      ar = n/k * ar
      # ar = (1 + n_k/k/2) * ar
      E_T_sim, E_C_sim, E_T_tpar, E_C_tpar, E_T_par, E_C_par = 0, 0, 0, 0, 0, 0
      for f in range(num_frun):
        s_l = sim(ar, h)
        if s_l is None: done = True; break
        
        taskt_rv = SimRV(s_l)
        stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(taskt_rv, 0, k, k, n, num_run=10000*10)
        E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
        E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
        print("Sim: E_T= {}, E_C= {}".format(E_T, E_C) )
        E_T_sim += E_T
        E_C_sim += E_C
        
        s_l = numpy.sort(s_l)[::-1]
        l, u, a = fit_tpareto(s_l)
        task_t = "TPareto"
        task_dist_m = {'l': l, 'u': u, 'a': a}
        E_T = E_T_k_l_n(task_t, task_dist_m, 0, k, k, n)
        E_C = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n, w_cancel=True)
        print("Fitted TPareto: E_T= {}, E_C= {}".format(E_T, E_C) )
        E_T_tpar += E_T
        E_C_tpar += E_C
        
        l, a = fit_pareto(s_l)
        task_t = "Pareto"
        task_dist_m = {'loc': l, 'a': a}
        E_T = E_T_k_l_n(task_t, task_dist_m, 0, k, k, n)
        E_C = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n, w_cancel=True)
        print("Fitted Pareto: E_T= {}, E_C= {}".format(E_T, E_C) )
        E_T_par += E_T
        E_C_par += E_C
      if not done:
        x_sim_l.append(E_T_sim/num_frun)
        y_sim_l.append(E_C_sim/num_frun)
        x_tpar_l.append(E_T_tpar/num_frun)
        y_tpar_l.append(E_C_tpar/num_frun)
        x_par_l.append(E_T_par/num_frun)
        y_par_l.append(E_C_par/num_frun)
      n_k += 1
    plot.plot(x_sim_l[0], y_sim_l[0], label=r'No redundancy', zorder=2, marker='x', color='blue', mew=3, ms=9)
    plot.plot(x_sim_l, y_sim_l, label=r'Simulation', zorder=1, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    # plot.legend()
    # plot.xscale('log')
    # plot.yscale('log')
    # plot.xlabel(r'$E[T]$', fontsize=13)
    # plot.ylabel(r'$E[C]$', fontsize=13)
    # plot.title(r'$T \sim {}$, $k= {}$'.format(proc_in_latex, k) )
    # plot.savefig("plot_EC_vs_ET_sim.png" )
    # plot.gcf().clear()
    
    plot.plot(x_tpar_l[0], y_tpar_l[0], zorder=2, marker='x', color='blue', mew=3, ms=9)
    plot.plot(x_tpar_l, y_tpar_l, zorder=0, label=r'Using fitted Truncated-Pareto', color=next(dark_color), linestyle='-')
    plot.plot(x_par_l[0], y_par_l[0], zorder=2, marker='x', color='blue', mew=3, ms=9)
    plot.plot(x_par_l, y_par_l, zorder=0, label=r'Using fitted Pareto', color=next(dark_color), linestyle='-.')
    # plot.legend()
    # plot.xscale('log')
    # plot.yscale('log')
    # plot.xlabel(r'$E[T]$', fontsize=13)
    # plot.ylabel(r'$E[C]$', fontsize=13)
    # plot.title(r'$T \sim {}$, $k= {}$'.format(proc_in_latex, k) )
    # plot.savefig("plot_EC_vs_ET_model.png" )
    # plot.gcf().clear()
  
  num_frun = 3 # 5
  h = 8
  k = 100
  ar = 0.01 # 0.05 # 0.1*ar_ub
  plot_EC_vs_ET(num_frun, h, ar, k)
  
  plot.legend()
  plot.xscale('log')
  plot.yscale('log')
  plot.xlabel(r'$E[T]$', fontsize=13)
  plot.ylabel(r'$E[C]$', fontsize=13)
  # plot.title(r'$T \sim {}$, $k= {}$'.format(proc_in_latex, k) )
  plot.title(r'$k= {}$'.format(k) )
  plot.savefig("plot_EC_vs_ET.png")
  plot.gcf().clear()
  log(WARNING, "done; k= {}".format(k) )

if __name__ == "__main__":
  # plot_psq()
  plot_psq_tail()
  # MG1_T()
  
  # plot_EC_vs_ET_wsim()
