import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

from sim import *
import sys, pprint, math, numpy, simpy, sympy, getopt, itertools

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
    
    self.qt_l = []
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
        # self.qt_l.append(lifetime)
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

def fit_tpareto(sample_l):
  # sample_l = numpy.sort(sample_l)[::-1]
  n = len(sample_l)
  print("n= {}".format(n) )
  l = sample_l[-1]
  u = sample_l[0]
  
  # a = sympy.Symbol('a')
  # a_ = sympy.solve(n/a + n*r**a*math.log(r)/(1-r**a) - sum([math.log(x/l) for x in sample_l] ) )
  def eq(a):
    # Fitting whole tail
    r = l/u
    return n/a + n*r**a*math.log(r)/(1-r**a) - sum([math.log(x/l) for x in sample_l] )
    
    # # Fitting upper tail
    # i = int(n*0.8)
    # X_ip1 = sample_l[i+1]
    # def l(a):
    #   return i**(1/a) * X_ip1*(n - (n-i)*(X_ip1/u)**a)**(-1/a)
    # r = X_ip1/u
    # return i/a + i*r**a*math.log(r)/(1-r**a) - sum([math.log(x/X_ip1) for x in sample_l[:i+1] ] )
  
  a_ = None
  a_l, eq_l = [], []
  for a in numpy.linspace(0.05, 100, 1000):
    # a_l.append(a)
    # eq_l.append(eq(a) )
    if eq(a) <= 0:
      a_ = a
      break
  # plot.plot(a_l, eq_l, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot.xlabel(r'$a$', fontsize=12)
  # plot.ylabel(r'$eq$', fontsize=12)
  # plot.savefig("fit_tpareto.png")
  log(WARNING, "done.")
  
  return l, u, a_

def plot_psq_tail():
  # D, mu = 1, 0.5
  # psize_dist = Exp(mu, D=D)
  # proc_in_latex = "Exp(D={}, \mu={})".format(D, mu)
  
  # l, a = 1, 1.5 # 4
  # psize_dist = Pareto(l, a)
  # proc_in_latex = r'Pareto(\lambda={}, \alpha={})'.format(l, a)
  
  l, u, a = 1, 100, 1.5
  psize_dist = TPareto(l, u, a)
  proc_in_latex = r'TPareto(l={}, u={}, \alpha={})'.format(l, u, a)
  
  # l, u = 1, 100
  # psize_dist = DUniform(l, u)
  # proc_in_latex = r'DUniform(l={}, u={})'.format(l, u)
  log(WARNING, "psize_dist= {}".format(psize_dist) )
  
  ar_ub = 1/psize_dist.mean()
  
  ar_l, fitted_a_l = [], []
  def plot_(ar, h):
    log(WARNING, "ar= {}, h= {}".format(ar, h) )
    
    env = simpy.Environment()
    pg = PG(env, ar, psize_dist)
    psq = PSQ(env, h)
    pg.out = psq
    pg.init()
    env.run(until=50000*10)
    
    sl_l = numpy.sort(psq.slowdown_l)
    x_l = sl_l[::-1]
    # y_sim_l = numpy.arange(sl_l.size)/sl_l.size
    # plot.plot(x_l, y_sim_l, label="ar= {}".format(ar), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    
    i_ = None
    for i in range(len(x_l)-1, 0, -1):
      if x_l[i] > 1.01:
        i_ = i
        break
    x_l = x_l[:i_]
    y_sim_l = numpy.arange(x_l.size)/x_l.size
    plot.plot(x_l, y_sim_l, label="ar= {}".format(ar), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    
    l, u, a = fit_tpareto(x_l)
    print("l= {}, u= {}, a= {}".format(l, u, a) )
    ar_l.append(ar)
    fitted_a_l.append(a)
    # rv = TPareto(l, u, a)
    # y_l = []
    # for x in x_l:
    #   y_l.append(rv.tail(x) )
    # plot.plot(x_l, y_l, label="fitted, ar= {}".format(ar), color=next(dark_color), linestyle='-')
    
    plot.legend()
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel(r'Slowdown', fontsize=13)
    plot.ylabel(r'Tail distribution', fontsize=13)
    plot.title(r'$P \sim {}$, $\lambda$= {}'.format(proc_in_latex, ar) )
    plot.savefig("plot_psq_h_{}_ar_{}.png".format(h, ar) )
    plot.gcf().clear()
  
  h = 4
  for ar in numpy.linspace(0.05, 1.5*ar_ub, 5):
    plot_(ar, h)
  
  plot.plot(ar_l, fitted_a_l, label="h= {}".format(h), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.xlabel(r'$\lambda$', fontsize=13)
  plot.ylabel(r'$\alpha$', fontsize=13)
  plot.title(r'$P \sim {}$'.format(proc_in_latex) )
  plot.savefig("plot_fitted_tail_h_{}.png".format(h) )
  plot.gcf().clear()
  
  log(WARNING, "done; psize_dist= {}".format(psize_dist) )

if __name__ == "__main__":
  # plot_psq()
  plot_psq_tail()
