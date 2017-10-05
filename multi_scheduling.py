import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

import sys, pprint, math, numpy, simpy, getopt, itertools

from sim import *
from arepeat_models import *

class Job(object):
  def __init__(self, _id, k):
    self._id = _id
    self.k = k
    self.prev_hop_id = None
  
  def __repr__(self):
    return "Job[id= {}, k= {}]".format(self._id, self.k)
  
  def deep_copy(self):
    j = Job(self._id, self.k)
    j.prev_hop_id = self.prev_hop_id
    return j

class JG(object): # Job Generator
  def __init__(self, env, ar, k_dist):
    self.env = env
    self.ar = ar
    self.k_dist = k_dist
    
    self.nsent = 0
    self.out = None
    
    self.action = None
  
  def init(self):
    self.action = self.env.process(self.run() )
  
  def run(self):
    while 1:
      yield self.env.timeout(random.expovariate(self.ar) )
      
      self.nsent += 1
      k = self.k_dist.gen_sample()
      self.out.put(Job(self.nsent, k) )

class JQ(object):
  def __init__(self, env, in_qid_l):
    self.env = env
    self.in_qid_l = in_qid_l
    
    self.jid__t_l_map = {}
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
  
  def __repr__(self):
    return "JQ[in_qid_l= {}]".format(self.in_qid_l)
  
  def run(self):
    while True:
      t = (yield self.store.get() )
      if t._id not in self.jid__t_l_map:
        self.jid__t_l_map[t._id] = []
      self.jid__t_l_map[t._id].append(t.deep_copy() )
      # check for completion
      # print("Checking completion; jid__t_l_map=\n {}".format(self.jid__t_l_map) )
      t_l = self.jid__t_l_map[t._id]
      if len(t_l) > t.k:
        log(ERROR, "len(t_l)= {} > k= {}".format(len(t_l), t.k) )
      elif len(t_l) < t.k:
        continue
      else:
        self.jid__t_l_map.pop(t._id, None)
        self.out_c.put_c({'jid': t._id, 'm': 'done', 'departed_from': [t.prev_hop_id for t in t_l] } )
  
  def put(self, t):
    sim_log(DEBUG, self.env, self, "recved", t)
    return self.store.put(t)

class TServer(object): # Task
  def __init__(self, _id, env, serv, serv_dist_m, out, out_c):
    self._id = _id
    self.env = env
    self.out = out
    self.out_c = out_c
    
    if serv == "SExp":
      self.serv_time = Exp(serv_dist_m['mu'], D=serv_dist_m['D'] )
    elif serv == "Pareto":
      self.serv_time = Pareto(serv_dist_m['loc'], serv_dist_m['a'] )
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    
    self.t_inserv = None
    self.cancel = None
    self.cancel_flag = False
  
  def __repr__(self):
    return "TServer[id = {}]".format(self._id)
  
  def idle(self):
    return (self.t_inserv is None)
  
  def run(self):
    while True:
      self.t_inserv = (yield self.store.get() )
      
      self.cancel = self.env.event()
      yield (self.cancel | self.env.timeout(self.serv_time.gen_sample() ) )
      if self.cancel_flag:
        sim_log(DEBUG, self.env, self, "cancelled", self.t_inserv)
        self.cancel_flag = False
      else:
        self.t_inserv.prev_hop_id = self._id
        self.out.put(self.t_inserv)
      self.out_c.put_c({'qid': self._id, 'm': 'idle'} )
      self.t_inserv = None
  
  def put(self, t):
    sim_log(DEBUG, self.env, self, "recved", t)
    return self.store.put(t)
  
  def put_c(self, m):
    if m['m'] == 'done':
      if self.t_inserv.jid != m['jid']:
        log(ERROR, "t_inserv.jid= {} != m['jid']= {}".format(self.t_inserv.jid, m['jid']) )
    self.cancel_flag = True
    self.cancel.succeed()

class MultiQ(object):
  def __init__(self, env, N, sching_m, serv, serv_dist_m):
    self.env = env
    self.N = N
    self.sching_m = sching_m
    
    self.sching_t = self.sching_m["t"]
    
    self.jq = JQ(env, range(self.N) )
    self.jq.out_c = self
    self.s_l = [TServer(i, env, serv, serv_dist_m, out=self.jq, out_c=self) for i in range(self.N) ]
    
    self.j_l = []
    self.nidleq = self.N
    self.waiting_for_nidleq = None
    self.jid_qid_l = {}
    self.action = env.process(self.run() )
    
    self.jid_enttime = {}
    self.qt_l = []
    self.nj_departed = 0
  
  def __repr__(self):
    return "MultiQ[N = {}]".format(self.N)
  
  def length(self):
    return len(self.j_l)
  
  def idle_qid_l(self):
    l = []
    for i, q in enumerate(self.s_l):
      if q.idle():
        l.append(i)
    return l
  
  def run(self):
    while True:
      if len(self.j_l) == 0:
        self.got_busy = self.env.event()
        yield (self.got_busy)
        self.got_busy = None
        sim_log(DEBUG, self.env, self, "got busy!", None)
        continue
      
      j = self.j_l.pop(0)
      if self.sching_t == "red":
        nidleq_needed = j.k + self.sching_m["n-k"]
        if self.nidleq < nidleq_needed:
          self.waiting_for_nidleq = nidleq_needed
          sim_log(DEBUG, self.env, self, "waiting for nidleq= {} for".format(nidleq_needed), j)
          self.got_nidleq = self.env.event()
          yield (self.got_nidleq)
          self.waiting_for_nidleq = None
          sim_log(DEBUG, self.env, self, "got nidleq= {} for".format(nidleq_needed), j)
        toi_l = self.idle_qid_l()[:nidleq_needed]
      elif self.sching_t == "red-to-idle":
        pass
      self.jid_qid_l[j._id] = toi_l
      self.nidleq -= nidleq_needed
      for i in toi_l:
        self.s_l[i].put(j.deep_copy() )
  
  def put(self, j):
    sim_log(DEBUG, self.env, self, "recved", j)
    self.jid_enttime[j._id] = self.env.now
    
    _l = len(self.j_l)
    self.j_l.append(j)
    if _l == 0 and self.got_busy is not None:
      self.got_busy.succeed()
  
  def put_c(self, m):
    sim_log(DEBUG, self.env, self, "recved", m)
    if 'qid' in m and m['m'] == 'idle':
      self.nidleq += 1
      if self.waiting_for_nidleq is not None and self.waiting_for_nidleq == self.nidleq:
        self.got_nidleq.succeed()
    elif 'jid' in m and m['m'] == 'done':
      jid = m['jid']
      self.qt_l.append(self.env.now - self.jid_enttime[jid] )
      self.jid_enttime.pop(jid, None)
      self.nj_departed += 1
      
      departed_from_qid_l = m['departed_from']
      for i in self.jid_qid_l[jid]:
        if i not in departed_from_qid_l:
          self.s_l[i].put_c({'jid': jid, 'm': 'cancel'} )
      self.jid_qid_l.pop(jid)

def W_M_G_c(c, ar, E_S, E_S_2):
  # from https://en.wikipedia.org/wiki/M/G/k_queue
  def W_M_M_c(c, ar, mu):
    # from https://en.wikipedia.org/wiki/M/M/c_queue
    ro = ar/mu/c
    denom = 1 + (1-ro)*G(c+1)/(c*ro)**c * sum([(c*ro)**k/G(k+1) for k in range(c) ] )
    C = 1/denom
    return C/(c*mu - ar)
  coeff_var = math.sqrt(E_S_2 - E_S**2)/E_S
  return (coeff_var**2 + 1)/2 * W_M_M_c(c, ar, 1/E_S)

def E_T_multiq(ar, N, k, n, serv, serv_dist_m):
  E_S = E_T_k_l_n(serv, serv_dist_m, 0, k, k, n)
  E_S_2 = E_T_2_k(serv, serv_dist_m, k, n=n)
  E_X = E_T_k_l_n(serv, serv_dist_m, 0, k, k, N)
  E_X_2 = E_T_2_k(serv, serv_dist_m, k, n=N)
  
  num_group = int(N/n)
  # ar = ar/num_group
  # return E_S + ar*E_S_2/2/(1 - ar*E_S)
  
  # return E_S + W_M_G_c(num_group, ar, E_S, E_S_2)
  # return E_S + W_M_G_c(num_group, ar, E_X, E_X_2)
  
  # E_S = n/N * E_C_k_l_n(serv, serv_dist_m, 0, k, k, n, w_cancel=True)
  # E_S_2 = (n/N)**2 * E_C_2_k(serv, serv_dist_m, k, n=n)
  # return E_S + ar*E_S_2/2/(1 - ar*E_S)
  
  return E_S + ar*E_X_2/2/(1 - ar*E_X)

# *********************************  Sim  *********************************** #
def sim_multiq(num_f_run, ar, N, sching_m, k_dist, serv, serv_dist_m):
  E_T_sum = 0
  for f in range(num_f_run):
    log(WARNING, "ar= {}, N= {}, sching_m= {}".format(ar, N, sching_m) )
    
    env = simpy.Environment()
    jg = JG(env, ar, k_dist)
    mq = MultiQ(env, N, sching_m, serv, serv_dist_m)
    jg.out = mq
    jg.init()
    qm = QMonitor(env, mq, poll_interval= 0.1)
    env.run(until=50000)
    
    l = mq.qt_l
    if len(l): E_T_sum += float(sum(l) )/len(l)
    print("jg.nsent= {}, mq.nj_departed= {}".format(jg.nsent, mq.nj_departed) )
    
    # plot.plot(qm.pollt_l, qm.qlength_l, label='ar= {}'.format(ar), color=next(dark_color), marker=next(marker), linestyle='--', zorder=0, mew=1)
    # plot.legend()
    # plot.xlabel(r'time', fontsize=12)
    # plot.ylabel(r'Q length', fontsize=12)
    # plot.title('sching_m= {}'.format(sching_m) )
    # plot.savefig("plot_sim_multiq_n_k_{}.png".format(sching_m['n-k'] ), bbox_inches='tight')
    # plot.gcf().clear()
  E_T = E_T_sum/num_f_run
  print(">> E_T_sim= {}".format(E_T) )
  return E_T

def plot_multiq():
  N = 100
  # k_dist = DUniform(4, 4)
  k_dist = DUniform(2, 20)
  serv = "Pareto" # "SExp" # "Pareto"
  D, mu = 0, 1
  loc, a = 3, 1.5
  if serv == "SExp":
    serv_dist_m = {'D': D, 'mu': mu}
    serv_in_latex = "Exp(\mu={})".format(mu)
    serv_time = Exp(mu)
  elif serv == "Pareto":
    serv_dist_m = {'loc': loc, 'a': a}
    serv_in_latex = r'Pareto(\lambda={}, \alpha={})'.format(loc, a)
    serv_time = Pareto(loc, a)
  ar_ub = 0.7 * N/serv_time.mean()/k_dist.mean()
  log(WARNING, "N= {}, ar_ub= {}, k_dist= {}, serv_dist_m= {}".format(N, ar_ub, k_dist, serv_dist_m) )
  
  x_l, y_l, y_sim_l = [], [], []
  
  sim = False
  if serv == "SExp":
    if N == 10:
      pass
  else:
    sim = True
  
  num_f_run = 1
  def plot_E_T_vs_ar(n_k):
    for ar in numpy.linspace(0.05, ar_ub, 6):
    # for ar in numpy.linspace(0.1, 0.1, 1):
      x_l.append(ar)
      if sim:
        sching_m = {"t": "red", "n-k": n_k}
        E_T = sim_multiq(num_f_run, ar, N, sching_m, k_dist, serv, serv_dist_m)
        y_sim_l.append(E_T)
      
      k = int(k_dist.mean() )
      E_T = E_T_multiq(ar, N, k, k + n_k, serv, serv_dist_m)
      print("E_T= {}".format(E_T) )
      y_l.append(E_T)
      
    print("y_sim_l= \n{}".format(pprint.pformat(y_sim_l) ) )
    plot.plot(x_l, y_sim_l, label=r'sim, $n-k$= {}'.format(n_k), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(x_l, y_l, label=r'$n-k$= {}'.format(n_k), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$\lambda$', fontsize=12)
  
  def plot_E_T_vs_n_k(ar):
    for n_k in range(N-k_dist.u_l):
      sching_m = {"t": "red", "n-k": n_k}
      E_T = sim_multiq(num_f_run, ar, N, sching_m, k_dist, serv, serv_dist_m)
      if E_T > 100:
        break
      x_l.append(n_k)
      y_l.append(E_T)
    plot.plot(x_l, y_l, label=r'Always red', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$n-k$', fontsize=12)
  
  # plot_E_T_vs_ar(n_k=1)
  plot_E_T_vs_n_k(ar=ar_ub*0.75)
  
  plot.ylabel(r'$E[T]$', fontsize=12)
  plot.legend(prop={'size':11} )
  plot.title(r'$S \sim {}$, $N={}$'.format(serv_in_latex, N) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_multiq_N_{}.png".format(N) )
  log(WARNING, "done; N= {}".format(N) )

if __name__ == "__main__":
  plot_multiq()
  