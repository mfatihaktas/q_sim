import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

import sys, pprint, math, numpy, simpy, getopt, itertools, operator

from rvs import *
from arepeat_models import *
from arepeat_sim import *

# ************************  Multiple Qs for Jobs with multiple Tasks  **************************** #
class Job(object):
  def __init__(self, _id, k, tsize, n=0):
    self._id = _id
    self.k = k
    self.tsize = tsize
    self.n = n
    self.prev_hop_id = None
  
  def __repr__(self):
    return "Job[id= {}, k= {}, n= {}]".format(self._id, self.k, self.n)
  
  def deep_copy(self):
    j = Job(self._id, self.k, self.tsize, self.n)
    j.prev_hop_id = self.prev_hop_id
    return j

class JG(object): # Job Generator
  def __init__(self, env, ar, k_dist, tsize_dist):
    self.env = env
    self.ar = ar
    self.k_dist = k_dist
    self.tsize_dist = tsize_dist
    
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
      tsize = self.tsize_dist.gen_sample()
      self.out.put(Job(self.nsent, k, tsize) )

class Task(object):
  def __init__(self, jid, k, size, remaining):
    self.jid = jid
    self.k = k
    self.size = size
    self.remaining = remaining
    
    self.prev_hop_id = None
    self.ent_time = None
  
  def __repr__(self):
    return "Task[jid= {}, k= {}, size= {}, remaining= {}]".format(self.jid, self.k, self.size, self.remaining)
  
  def deep_copy(self):
    t = Task(self.jid, self.k, self.size, self.remaining)
    t.prev_hop_id = self.prev_hop_id
    t.ent_time = self.ent_time
    return t

class PSQ(object): # Process Sharing Queue
  def __init__(self, _id, env, h, out):
    self._id = _id
    self.env = env
    self.h = h
    self.out = out
    
    self.t_l = []
    self.tinserv_l = []
    self.got_busy = None
    self.sinterrupt = None
    self.add_to_serv = False
    self.cancel = False
    self.cancel_jid = None
    
    self.store = simpy.Store(env)
    self.action = env.process(self.serv_run() )
    self.action = env.process(self.put_run() )
    
    self.lt_l = []
    self.sl_l = []
  
  def __repr__(self):
    return "PSQ[id= {}]".format(self._id)
  
  def length(self):
    return len(self.t_l)
  
  def serv_run(self):
    while True:
      self.tinserv_l = self.t_l[:self.h]
      if len(self.tinserv_l) == 0:
        # sim_log(DEBUG, self.env, self, "idle; waiting for arrival", None)
        self.got_busy = self.env.event()
        yield (self.got_busy)
        # sim_log(DEBUG, self.env, self, "got busy!", None)
        continue
      # TODO: This seems wrong
      # t_justmovedHoL = self.tinserv_l[-1]
      # self.out.put_c({'m': 'HoL', 'jid': t_justmovedHoL.jid, 'k': t_justmovedHoL.k, 'qid': self._id} )
      
      serv_size = len(self.tinserv_l)
      r_l = [self.tinserv_l[i].remaining for i in range(serv_size) ]
      time = min(r_l)
      i_min = r_l.index(time)
      sim_log(DEBUG, self.env, self, "back to serv; time= {}, serv_size= {}".format(time, serv_size), None)
      start_t = self.env.now
      
      self.sinterrupt = self.env.event()
      yield (self.sinterrupt | self.env.timeout(time) )
      serv_t = (self.env.now - start_t)/serv_size
      for i in range(serv_size):
        try:
          self.t_l[i].remaining -= serv_t
        except IndexError:
          break
      
      if self.add_to_serv:
        sim_log(DEBUG, self.env, self, "new task added to serv", None)
        self.sinterrupt = None
        self.add_to_serv = False
      elif self.cancel:
        # sim_log(DEBUG, self.env, self, "will cancel task in serv", t)
        for t in self.t_l:
          if t.jid == self.cancel_jid:
            sim_log(DEBUG, self.env, self, "cancelled task in serv", t)
            self.t_l.remove(t)
        self.sinterrupt = None
        self.cancel = False
      else:
        t = self.t_l.pop(i_min)
        sim_log(DEBUG, self.env, self, "serv done", t)
        
        lt = self.env.now - t.ent_time
        self.lt_l.append(lt)
        self.sl_l.append(lt/t.size)
        
        t.prev_hop_id = self._id
        self.out.put(t)
  
  def put_run(self):
    while True:
      t = (yield self.store.get() )
      _l = len(self.t_l)
      self.t_l.append(t)
      if _l == 0:
        self.got_busy.succeed()
      elif _l < self.h:
        self.add_to_serv = True
        self.sinterrupt.succeed()
  
  def put(self, t):
    sim_log(DEBUG, self.env, self, "recved", t)
    t.ent_time = self.env.now
    return self.store.put(t) # .deep_copy()
  
  def put_c(self, m):
    sim_log(DEBUG, self.env, self, "recved; tinserv_l= {}".format(self.tinserv_l), m)
    
    if m['m'] == 'cancel':
      jid = m['jid']
      if jid in [t.jid for t in self.tinserv_l]:
        self.cancel = True
        self.cancel_jid = jid
        self.sinterrupt.succeed()
      else:
        for t in self.t_l:
          if t.jid == jid:
            self.t_l.remove(t)

class FCFS(object):
  def __init__(self, _id, env, slowdown_dist, out):
    self._id = _id
    self.env = env
    self.slowdown_dist = slowdown_dist
    self.out = out
    
    self.t_l = []
    self.t_inserv = None
    self.got_busy = None
    self.cancel_flag = False
    self.cancel = None
    
    self.sl_l = []
    
    self.action = env.process(self.serv_run() )
  
  def __repr__(self):
    return "FCFS[_id= {}]".format(self._id)
  
  def length(self):
    return len(self.t_l) + (self.t_inserv is not None)
  
  def serv_run(self):
    while True:
      if len(self.t_l) == 0:
        self.got_busy = self.env.event()
        yield (self.got_busy)
        self.got_busy = None
        # sim_log(DEBUG, self.env, self, "got busy!", None)
      self.t_inserv = self.t_l.pop(0)
      
      self.cancel = self.env.event()
      clk_start_time = self.env.now
      st = self.t_inserv.size * self.slowdown_dist.gen_sample()
      sim_log(DEBUG, self.env, self, "starting {}s-clock on ".format(st), self.t_inserv)
      yield (self.cancel | self.env.timeout(st) )
      
      if self.cancel_flag:
        sim_log(DEBUG, self.env, self, "cancelled clock on ", self.t_inserv)
        self.cancel_flag = False
      else:
        sim_log(DEBUG, self.env, self, "serv done in {}s on ".format(self.env.now-clk_start_time), self.t_inserv)
        self.sl_l.append((self.env.now - self.t_inserv.ent_time)/self.t_inserv.size)
      
        self.t_inserv.prev_hop_id = self._id
        self.out.put(self.t_inserv)
      self.t_inserv = None
  
  def put(self, t):
    sim_log(DEBUG, self.env, self, "recved", t)
    _l = len(self.t_l)
    t.ent_time = self.env.now
    self.t_l.append(t) # .deep_copy()
    if self.got_busy is not None and _l == 0:
      self.got_busy.succeed()
  
  def put_c(self, m):
    sim_log(DEBUG, self.env, self, "recved", m)
    
    if m['m'] == 'cancel':
      jid = m['jid']
      for t in self.t_l:
        if t.jid == jid:
          self.t_l.remove(t)
      if jid == self.t_inserv.jid:
        self.cancel_flag = True
        self.cancel.succeed()

class JQ(object):
  def __init__(self, env, in_qid_l):
    self.env = env
    self.in_qid_l = in_qid_l
    
    self.jid__t_l_map = {}
    self.deped_jid_l = []
    
    self.jid_HoLqid_l_map = {}
    self.movedHoL_jid_l = []
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    
    self.store_c = simpy.Store(env)
    self.action = env.process(self.run_c() )
  
  def __repr__(self):
    return "JQ[in_qid_l= {}]".format(self.in_qid_l)
  
  def run(self):
    while True:
      t = (yield self.store.get() )
      if t.jid in self.deped_jid_l: # Redundant tasks of a job may be received
        continue
      
      if t.jid not in self.jid__t_l_map:
        self.jid__t_l_map[t.jid] = []
      self.jid__t_l_map[t.jid].append(t.deep_copy() )
      
      t_l = self.jid__t_l_map[t.jid]
      if len(t_l) > t.k:
        log(ERROR, "len(t_l)= {} > k= {}".format(len(t_l), t.k) )
      elif len(t_l) < t.k:
        continue
      else:
        self.jid__t_l_map.pop(t.jid, None)
        self.deped_jid_l.append(t.jid)
        self.out_c.put_c({'jid': t.jid, 'm': 'jdone', 'deped_from': [t.prev_hop_id for t in t_l] } )
  
  def put(self, t):
    sim_log(DEBUG, self.env, self, "recved", t)
    return self.store.put(t)
  
  def run_c(self):
    while True:
      m = (yield self.store_c.get() )
      if m['m'] == 'HoL':
        jid, k, qid = m['jid'], m['k'], m['qid']
        if m['jid'] in self.movedHoL_jid_l: # Redundant tasks may move HoL simultaneously
          continue
        
        if jid not in self.jid_HoLqid_l_map:
          self.jid_HoLqid_l_map[jid] = []
        self.jid_HoLqid_l_map[jid].append(qid)
        
        HoLqid_l = self.jid_HoLqid_l_map[jid]
        if len(HoLqid_l) > k:
          log(ERROR, "len(HoLqid_l)= {} > k= {}".format(len(HoLqid_l), k) )
        elif len(HoLqid_l) < k:
          continue
        else:
          self.movedHoL_jid_l.append(jid)
          HoLqid_l = self.jid_HoLqid_l_map[jid]
          self.out_c.put_c({'m': 'jHoL', 'jid': jid, 'at': HoLqid_l} )
          self.jid_HoLqid_l_map.pop(jid, None)
  
  def put_c(self, m):
    sim_log(DEBUG, self.env, self, "recved", m)
    return self.store_c.put(m)

class MultiQ(object):
  def __init__(self, env, N, sching_m):
    self.env = env
    self.N = N
    self.sching_m = sching_m
    
    self.jq = JQ(env, range(self.N) )
    self.jq.out_c = self
    self.q_l = [PSQ(i, env, h=8, out=self.jq) for i in range(self.N) ]
    # slowdown_dist = DUniform(1, 1) # Dolly()
    # self.q_l = [FCFS(i, env, slowdown_dist, out=self.jq) for i in range(self.N) ]
    
    self.jid_info_m = {}
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    
    self.jsl_l = []
  
  def __repr__(self):
    return "MultiQ[N= {}]".format(self.N)
  
  def tlt_l(self):
    l = []
    for q in self.q_l: l.extend(q.lt_l)
    return l
  
  def get_sorted_qids(self):
    qid_length_m = {q._id: q.length() for q in self.q_l}
    # print("qid_length_m= {}".format(qid_length_m) )
    qid_length_l = sorted(qid_length_m.items(), key=operator.itemgetter(1) )
    # print("qid_length_l= {}".format(qid_length_l) )
    return [qid_length[0] for qid_length in qid_length_l]
  
  def run(self):
    while True:
      j = (yield self.store.get() )
      # toi_l = random.sample(range(self.N), j.n)
      toi_l = self.get_sorted_qids()[:j.n]
      
      for i in toi_l:
        self.q_l[i].put(Task(j._id, j.k, j.tsize, j.tsize) )
      self.jid_info_m[j._id] = {'ent_time': self.env.now, 'tsize': j.tsize, 'qid_l': toi_l}
  
  def put(self, j):
    sim_log(DEBUG, self.env, self, "recved", j)
    if self.sching_m['t'] == 'coded':
      # j.n = j.k + self.sching_m['n-k']
      j.n = math.floor(j.k*self.sching_m['r'] )
    return self.store.put(j.deep_copy() )
  
  def put_c(self, m):
    sim_log(DEBUG, self.env, self, "recved", m)
    if m['m'] == 'jHoL':
      jid = m['jid']
      jinfo = self.jid_info_m[jid]
      
      for i in jinfo['qid_l']:
        if i not in m['at']:
          self.q_l[i].put_c({'m': 'cancel', 'jid': jid} )
    elif m['m'] == 'jdone':
      jid = m['jid']
      jinfo = self.jid_info_m[jid]
      self.jsl_l.append((self.env.now - jinfo['ent_time'])/jinfo['tsize'] )
      
      for i in jinfo['qid_l']:
        if i not in m['deped_from']:
          self.q_l[i].put_c({'m': 'cancel', 'jid': jid} )
      
      self.jid_info_m.pop(jid, None)

# class MultiQ_RedToIdle(MultiQ):
#   def __init__(self, env, N, sching_m):
#     super().__init__(env, N, sching_m)
    
#     self.qid_redjid_l = {i: [] for i in range(N) }
  
#   def __repr__(self):
#     return "MultiQ_RedToIdle[N = {}]".format(self.N)
  
#   def get_sorted_qid_length_l(self):
#     qid_length_m = {q._id: q.length() for q in self.q_l}
#     qid_length_l = sorted(qid_length_m.items(), key=operator.itemgetter(1) )
#     return qid_length_l
  
#   def run(self):
#     while True:
#       j = (yield self.store.get() )
      
#       qid_length_l = self.get_sorted_qid_length_l()
#       toqid_l = []
#       for i, qid_length in enumerate(qid_length_l):
#         if i < j.k:
#           toqid_l.append(qid_length[0] )
#         elif i < j.n:
#           if qid_length[1] == 0:
#             qid = qid_length[0]
#             toqid_l.append(qid)
#             self.qid_redjid_l[qid].append(j._id)
      
#       for i, qid in enumerate(toqid_l):
#         if i < j.k:
#           for jid in self.qid_redjid_l[qid]:
#             self.q_l[qid].put_c({'m': 'cancel', 'jid': jid} )
#             self.qid_redjid_l[qid].remove(jid)
#         self.q_l[qid].put(Task(j._id, j.k, j.tsize, j.tsize) )
#       self.jid_info_m[j._id] = {'ent_time': self.env.now, 'tsize': j.tsize, 'qid_l': toqid_l}

# *********************************************  Sim  ******************************************** #
def sim_multiq(num_srun, ar, N, sching_m, k_dist, tsize_dist):
  k = 100
  n = math.floor(k*sching_m['r'] )
  E_T_sim, E_C_sim, E_T_tpar, E_C_tpar, E_T_par, E_C_par = 0, 0, 0, 0, 0, 0
  E_Sl_sum = 0
  for s in range(num_srun):
    log(WARNING, "s= {}, ar= {}, N= {}, sching_m= {}".format(s, ar, N, sching_m) )
    
    env = simpy.Environment()
    jg = JG(env, ar, k_dist, tsize_dist)
    if sching_m['t'] == 'red-to-idle':
      mq = MultiQ_RedToIdle(env, N, sching_m)
    else:
      mq = MultiQ(env, N, sching_m)
    jg.out = mq
    jg.init()
    env.run(until=50000) # *20
    
    E_jsl = sum(mq.jsl_l)/len(mq.jsl_l)
    print(">>> E_jsl= {}, jg.nj_sent= {}, mq.nj_departed= {}".format(E_jsl, jg.nsent, len(mq.jq.deped_jid_l) ) )
    if E_jsl > 8*1000: return None
    E_Sl_sum += E_jsl
    
    # Plotting the tail of task lifetimes
    s_l = numpy.sort(mq.tlt_l() )[::-1]
    i_ = None
    for i in range(len(s_l)-1, 0, -1):
      if s_l[i] > 1.01:
        i_ = i
        break
    s_l = s_l[:i_]
    y_sim_l = numpy.arange(s_l.size)/s_l.size
    # plot.plot(s_l, y_sim_l, label=r'n-k= {}'.format(sching_m['n-k'] ), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    
    taskt_rv = SimRV(s_l)
    stat_id__trial_sampleavg_l_m = sim_arepeat_k_l_n(taskt_rv, 0, k, k, n, num_run=10000*10)
    E_T = sum(stat_id__trial_sampleavg_l_m['T'] )/len(stat_id__trial_sampleavg_l_m['T'] )
    E_C = sum(stat_id__trial_sampleavg_l_m['C_wc'] )/len(stat_id__trial_sampleavg_l_m['C_wc'] )
    print("Sim: E_T= {}, E_C= {}".format(E_T, E_C) )
    E_T_sim += E_T
    E_C_sim += E_C
    
    l, u, a = fit_tpareto(s_l)
    # rv = TPareto(l, u, a)
    # y_l = []
    # for x in s_l:
    #   y_l.append(rv.tail(x) )
    # plot.plot(s_l, y_l, label="TPareto, n-k= {}, a= {}".format(sching_m['n-k'], a), color=next(dark_color), linestyle='-')
    task_t = "TPareto"
    task_dist_m = {'l': l, 'u': u, 'a': a}
    E_T = E_T_k_l_n(task_t, task_dist_m, 0, k, k, n)
    E_C = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n, w_cancel=True)
    print("Fitted TPareto: E_T= {}, E_C= {}".format(E_T, E_C) )
    E_T_tpar += E_T
    E_C_tpar += E_C
    
    l, a = fit_pareto(s_l)
    # rv = Pareto(l, a)
    # y_l = []
    # for x in s_l:
    #   y_l.append(rv.tail(x) )
    # plot.plot(s_l, y_l, label="Pareto, n-k= {}, a= {}".format(sching_m['n-k'], a), color=next(dark_color), linestyle='-')
    task_t = "Pareto"
    task_dist_m = {'loc': l, 'a': a}
    E_T = E_T_k_l_n(task_t, task_dist_m, 0, k, k, n)
    E_C = E_C_k_l_n(task_t, task_dist_m, 0, k, k, n, w_cancel=True)
    print("Fitted Pareto: E_T= {}, E_C= {}".format(E_T, E_C) )
    E_T_par += E_T
    E_C_par += E_C
    
    # plot.legend()
    # plot.xscale('log')
    # plot.yscale('log')
    # plot.xlabel(r'Lifetime', fontsize=13)
    # plot.ylabel(r'Tail distribution', fontsize=13)
    # plot.title(r'$T \sim {}$, $\lambda= {}$'.format(tsize_dist, ar) )
    # plot.savefig("plot_psq_n_k_{}.png".format(sching_m['n-k'] ) )
    # plot.gcf().clear()
  E_Sl = E_Sl_sum/num_srun
  E_T_sim, E_C_sim = E_T_sim/num_srun, E_C_sim/num_srun
  E_T_tpar, E_C_tpar = E_T_tpar/num_srun, E_C_tpar/num_srun
  E_T_par, E_C_par = E_T_par/num_srun, E_C_par/num_srun
  print(">> E_Sl= {}".format(E_Sl) )
  return {'E_Sl': E_Sl,
          'E_T_sim': E_T_sim, 'E_C_sim': E_C_sim,
          'E_T_tpar': E_T_tpar, 'E_C_tpar': E_C_tpar,
          'E_T_par': E_T_par, 'E_C_par': E_C_par}

def plot_multiq():
  N = 5 # 10
  k_dist = DUniform(1, 3) # DUniform(1, 4)
  l, u, a = 1, 10**10, 1.1
  tsize_dist = TPareto(l, u, a)
  tsize_in_latex = r'TPareto(l= {}, u= {}, \alpha= {})'.format(l, u, a)
  # ar_ub = 0.7 * N/tsize_dist.mean()/k_dist.mean()
  log(WARNING, "N= {}, k_dist= {}, tsize_dist= {}".format(N, k_dist, tsize_dist) )
  
  num_srun = 1
  def plot_ESl_vs_ar(sching_t, n_k):
    x_l, y_l = [], []
    sching_m = {'t': sching_t, 'n-k': n_k}
    ar = 0.05
    while True:
      sim_m = sim_multiq(num_srun, ar, N, sching_m, k_dist, tsize_dist)
      if sim_m is None: break
      E_Sl = sim_m['E_Sl']
      x_l.append(ar)
      y_l.append(E_Sl)
      ar += 0.1
    
    print("y_l=\n{}".format(pprint.pformat(y_l) ) )
    plot.plot(x_l, y_l, label=r'$n-k= {}$'.format(n_k), marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$\lambda$', fontsize=12)
    plot.ylabel(r'$E[Slowdown]$', fontsize=12)
  
  def plot_ESl_vs_n_k(sching_t, ar):
    x_l, y_l = [], []
    for n_k in range(N - k_dist.u_l):
      sching_m = {'t': sching_t, 'n-k': n_k}
      sim_m = sim_multiq(num_srun, ar, N, sching_m, k_dist, tsize_dist)
      if sim_m is None: break
      E_Sl = sim_m['E_Sl']
      x_l.append(n_k)
      y_l.append(E_Sl)
    plot.plot(x_l, y_l, label=sching_t, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.xlabel(r'$n-k$', fontsize=12)
    plot.ylabel(r'$E[Slowdown]$', fontsize=12)
  def plot_EC_vs_ET(ar):
    r_l = []
    x_sim_l, y_sim_l, x_tpar_l, y_tpar_l, x_par_l, y_par_l = [], [], [], [], [], []
    for r in numpy.arange(1, 10, 0.5):
      sching_m = {'t': 'coded', 'r': r}
      sim_m = sim_multiq(num_srun, ar, N, sching_m, k_dist, tsize_dist)
      if sim_m is None: break
      r_l.append(r)
      x_sim_l.append(sim_m['E_T_sim'] )
      y_sim_l.append(sim_m['E_C_sim'] )
      x_tpar_l.append(sim_m['E_T_tpar'] )
      y_tpar_l.append(sim_m['E_C_tpar'] )
      x_par_l.append(sim_m['E_T_par'] )
      y_par_l.append(sim_m['E_C_par'] )
    plot.plot(x_sim_l[0], y_sim_l[0], label=r'No redundancy', zorder=2, marker='x', color='blue', mew=3, ms=9)
    plot.plot(x_sim_l, y_sim_l, label=r'Simulation', zorder=1, marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(x_tpar_l, y_tpar_l, zorder=0, label=r'Using fitted Truncated-Pareto', color=next(dark_color), linestyle='-.', lw=2)
    plot.plot(x_par_l, y_par_l, zorder=0, label=r'Using fitted Pareto', color=next(dark_color), linestyle='-', lw=2)
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel(r'$E[T]$', fontsize=13)
    plot.ylabel(r'$E[C]$', fontsize=13)
  # plot_ESl_vs_ar('coded', n_k=0)
  # plot_ESl_vs_ar('coded', n_k=2)
  
  # plot_ESl_vs_ar('red-to-idle', n_k=2)
  # plot_ESl_vs_ar('coded', n_k=2)
  
  # plot_ESl_vs_n_k('red-to-idle', ar=0.2)
  # plot_ESl_vs_n_k('coded', ar=0.1)
  
  plot_EC_vs_ET(ar=0.4)
  # 
  plot.legend(prop={'size':11} )
  plot.title(r'$N= {}$, $k \sim {}$, $T \sim {}$'.format(N, k_dist, tsize_in_latex) )
  fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_multiq_N_{}.png".format(N) )
  fig.clear()
  log(WARNING, "done; N= {}".format(N) )

if __name__ == "__main__":
  plot_multiq()
  