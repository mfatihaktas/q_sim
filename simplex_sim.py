import simpy, random, copy, pprint
# from simpy.core import BoundClass
# from simpy.resources import base

from sim import *
from deprecated import *
from mds_sim import MDSQ
from patch import *

class MT_PG(PG):
  def __init__(self, env, _id, ar, sym_l, hot_sym=None, hot_ar=None, flow_id=0, fixed=False):
    super().__init__(env, _id, ar, flow_id)
    self.sym_l = sym_l
    self.hot_sym = hot_sym
    self.hot_ar = hot_ar
    self.fixed = fixed
    
    self.sym__n_sent = {}
  
  def init(self):
    self.env.process(self.run() )
    if self.hot_sym is not None:
      self.env.process(self.run_pop() )
  
  def send(self, p):
    self.out.put(p)
    
    self.n_sent += 1
    if p.sym not in self.sym__n_sent:
      self.sym__n_sent[p.sym] = 0
    self.sym__n_sent[p.sym] += 1
  
  def run(self):
    sym_l_ = list(self.sym_l)
    if self.hot_sym is not None:
      sym_l_.remove(self.hot_sym)
    
    while 1:
      yield self.env.timeout(random.expovariate(self.ar) )
      if not self.fixed:
        s = sym_l_[random.randint(0, len(sym_l_)-1) ]
      else:
        s = sym_l_[0]
      
      p = Packet(time=self.env.now, size=1, _id=self.n_sent, sym=s, flow_id=self.flow_id)
      self.send(p)
  
  def run_pop(self):
    while 1:
      yield self.env.timeout(random.expovariate(self.hot_ar) )
      
      p = Packet(time=self.env.now, size=1, _id=self.n_sent, sym=self.hot_sym, flow_id=self.flow_id)
      self.send(p)

class MT_AV_JQ(object):
  def __init__(self, _id, env, input_qid_l, sym__rgroup_l_m):
    self._id = _id
    self.env = env
    self.input_qid_l = input_qid_l
    self.sym__rgroup_l_m = sym__rgroup_l_m
    
    self.jid__p_l_m = {}
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.out = None
    self.out_c = None
    self.out_m = None
    self.action = env.process(self.run() )
    self.action = env.process(self.run_c() )
  
  def __repr__(self):
    return "MT_AV_JQ[_id= {}, input_qid_l= {}]".format(self._id, self.input_qid_l)
  
  def check_for_job_completion(self, p):
    p_l = self.jid__p_l_m[p.job_id]
    recved_from_qid_l = [p.prev_hop_id for p in p_l]
    rgroup_l = self.sym__rgroup_l_m[p.sym]
    success = False
    for rg in rgroup_l:
      success = True
      for qid in rg:
        if qid not in recved_from_qid_l:
          success = False
      if success:
        break
    if success:
      if len(rgroup_l) > 1:
        cp = CPacket(_id=p.job_id, sym=p.sym, prev_hop_id=self._id, departed_qid_l=recved_from_qid_l)
        # sim_log(WARNING, self.env, self, "p= {}, sending completion for".format(p), cp)
        self.out_c.put_c(cp)
      p.winner_id = p.prev_hop_id
      p.prev_hop_id = self._id
      self.out.put(p)
      if self.out_m is not None:
        self.out_m.put_m(MPacket(_id=p.job_id, event_str=MPACKET_JOB_DEPARTED) )
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      if p.job_id not in self.jid__p_l_m:
        self.jid__p_l_m[p.job_id] = []
      self.jid__p_l_m[p.job_id].append(p)
      self.check_for_job_completion(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.ref_time = self.env.now
    return self.store.put(p.deep_copy() )
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      self.jid__p_l_m.pop(p.job_id, None)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp.deep_copy() )

# **************************************  Availability Q  **************************************** #
"""
  Implements (n, k, r, t)-LRC Q
  k objects are encoded into n codewords.
  The (r, t)-availability ensures that each systematic node can be regenerated using one of the t
  disjoint repair groups of other servers, each of size at most r (typically << k)
  [When do the Availability Codes Make the Stored Data More Available?]
"""
class AVQ(object): # Availability
  def __init__(self, _id, env, t, r, k, servdist_m, sching, w_sys=True, out=None):
    self._id = _id
    self.env = env
    self.t = t
    self.r = r
    self.k = k
    self.sching = sching
    self.w_sys = w_sys
    self.out = out
    
    self.jsink = JSink(_id, env)
    self.jsink.out = out
    self.id_q_map = {}
    
    num_q = int(1 + t*r) if w_sys else int(t*r)
    qid_l = [i for i in range(num_q) ]
    self.qid_l = []
    if w_sys:
      for g in range(1, t + 1 + 1):
        if g == 1:
          self.qid_l.append(qid_l[0] )
        else:
          li = 1+(g-2)*r
          ri = li + r
          self.qid_l.append("mds{}".format(qid_l[li:ri] ) )
    else:
      for g in range(1, t + 1):
        li = (g-1)*r
        ri = li + r
        self.qid_l.append("mds{}".format(qid_l[li:ri] ) )
    self.join_q = JQ(_id, env, k=1, input_qid_l=self.qid_l)
    self.join_q.out = self.jsink
    self.join_q.out_c = self
    self.join_q.out_m = None # can be set by the caller if desired
    if w_sys:
      for g in range(1, t + 1 + 1):
        q = None
        if g == 1:
          q = FCFS(self.qid_l[g-1], env, servdist_m)
          log(DEBUG, "g= {}, q= {}".format(g, q) )
          q.out = self.join_q
        else:
          li = 1+(g-2)*r
          ri = li + r
          q = MDSQ("mds{}".format(list(qid_l[li:ri] ) ), env, k, qid_l[li:ri], servdist_m, out=self.join_q)
          log(DEBUG, "g= {}, q= {}".format(g, q) )
        self.id_q_map[g] = q
    else:
      for g in range(1, t + 1):
        li = (g-1)*r
        ri = li + r
        q = MDSQ("mds{}".format(list(qid_l[li:ri] ) ), env, k, qid_l[li:ri], servdist_m, out=self.join_q)
        log(DEBUG, "g= {}, q= {}, qmu_l[li:ri]= {}".format(g, q, qmu_l[li:ri] ) )
        self.id_q_map[g] = q
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.action = env.process(self.run() ) # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
    self.servtype__num_m = (t+1)*[0]
  
  def __repr__(self):
    return "AVQ[k= {}, r= {}, t= {}]".format(self.k, self.r, self.t)
    # return "AVQ[k= {}, r= {}, t= {}, qid_l= [{}] ]".format(self.k, self.r, ",".join(self.qid_l) )
  
  def state(self, job_id_to_exclude=[]):
    state = []
    for g, q in self.id_q_map.items():
      state += q.state(job_id_to_exclude)
    
    return state
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      if self.sching == "rep-to-all":
        for i, q in self.id_q_map.items():
          q.put(p.deep_copy() )
      else: # select-one at random
        id_l = [i for i,q in self.id_q_map.items() ]
        i = id_l[random.randint(0, len(id_l)-1) ]
        self.id_q_map[i].put(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.job_id_counter += 1
    p.job_id = self.job_id_counter
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      #
      next_job_id = cp._id + 1
      type_ = 0
      for g, q in self.id_q_map.items():
        if g == 1:
          if not q._in(next_job_id): break
        else:
          if q.n_servers_in(next_job_id) < self.r: type_ += 1
      self.servtype__num_m[type_] += 1
      #
      for g, q in self.id_q_map.items():
        if q._id not in cp.departed_qid_l:
          q.put_c(cp.deep_copy() )
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

class AVQMonitor(object):
  def __init__(self, env, q, poll_rate):
    self.q = q
    self.env = env
    self.poll_rate = poll_rate
    
    # self.t_l = [] # Time steps that the numbers polled from the q
    # self.n_l = [] # Num of jobs in the q
    self.polled_state__counter_map = {}
    self.state__num_found_by_job_departed_map = {}
    self.start_setup__num_found_by_job_departed_map = {}
    
    self.store = simpy.Store(env)
    env.process(self.run() )
    env.process(self.run_m() )
  
  def run(self):
    while True:
      yield self.env.timeout(1/self.poll_rate)
      # self.t_l.append(self.env.now)
      q_state = self.q.state()
      # print("q_state= {}".format(pprint.pformat(q_state) ) )
      sub_qs_state = [max(q_state) for i in range(len(q_state) ) ]
      num_jobs = max(sub_qs_state)
      for i, s in enumerate(sub_qs_state):
        sub_qs_state[i] -= q_state[i]
      
      # num_job = max(state)
      # self.n_l.append(num_job)
      # rel_state = list_to_str(sub_qs_state)
      rel_state = "{},({},{})".format(num_jobs, sub_qs_state[1], sub_qs_state[2] )
      if rel_state not in self.polled_state__counter_map:
        self.polled_state__counter_map[rel_state] = 0
      self.polled_state__counter_map[rel_state] += 1
  
  def run_m(self):
    while True:
      p = (yield self.store.get() )
      if p.event_str == MPACKET_JOB_DEPARTED:
        # state = list_to_str(self.q.join_q.state() )
        state = list_to_str(self.q.state([p._id] ) )
        if state not in self.state__num_found_by_job_departed_map:
          self.state__num_found_by_job_departed_map[state] = 0
        self.state__num_found_by_job_departed_map[state] += 1
        
        state = self.q.state([p._id] )
        start_setup = "partial"
        if len(set(state) ) <= 1:
          start_setup = "complete"
        if start_setup not in self.start_setup__num_found_by_job_departed_map:
          self.start_setup__num_found_by_job_departed_map[start_setup] = 0
        self.start_setup__num_found_by_job_departed_map[start_setup] += 1
      else:
        log(ERROR, "unexpected p.event_str= {}".format(p.event_str) )
        return 1
  
  def put_m(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    return self.store.put(p)
  
# *******************************  Mixed-Traffic Availability Q  ******************************* #
class MT_AVQ(object):
  def __init__(self, _id, env, t, sym__rgroup_l_m, servdist_m, sching='rep-to-all', out=None):
    self._id = _id
    self.env = env
    self.sym__rgroup_l_m = sym__rgroup_l_m
    self.sching = sching
    self.out = out
    
    self.num_q = int(1 + t*2)
    self.qid_l = [i for i in range(self.num_q) ]
    
    self.jsink = JSink(_id, env)
    self.jsink.out = out
    self.join_q = MT_AV_JQ(_id, env, self.qid_l, sym__rgroup_l_m)
    self.join_q.out = self.jsink
    self.join_q.out_c = self
    self.join_q.out_m = None # can be set by the caller if desired
    
    self.id_q_map = {}
    for i in self.qid_l:
      q = FCFS(i, env, servdist_m)
      log(DEBUG, "q= {}".format(q) )
      q.out = self.join_q
      self.id_q_map[i] = q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    env.process(self.run() )
    env.process(self.run_c() )
    
    self.job_id_counter = 0
  
  def __repr__(self):
    return "MT_AVQ[qid_l= {}]".format(self.qid_l)
  
  def state(self, job_id_to_exclude=[]):
    state = []
    for i, q in self.id_q_map.items():
      state += q.state(job_id_to_exclude)
    return state
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      for i, q in self.id_q_map.items():
        q.put(p.deep_copy() )
      
      if self.sching == 'rep-to-all':
        for i, q in self.id_q_map.items():
          q.put(p.deep_copy() )
      elif self.sching == 'select-one':
        id_l = [i for i,q in self.id_q_map.items() ]
        i = id_l[random.randint(0, len(id_l)-1) ]
        self.id_q_map[i].put(p)
      
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.job_id_counter += 1
    p.job_id = self.job_id_counter
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      for i, q in self.id_q_map.items():
        if q._id not in cp.departed_qid_l:
          q.put_c(cp.deep_copy() )
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

class AVQMonitor(object):
  def __init__(self, env, q, poll_rate):
    self.q = q
    self.env = env
    self.poll_rate = poll_rate
    
    self.t_l = [] # Time steps that the numbers polled from the q
    # self.n_l = [] # Num of jobs in the q
    self.polled_state__counter_map = {}
    self.state__num_found_by_job_departed_map = {}
    self.start_setup__num_found_by_job_departed_map = {}
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    self.action = env.process(self.run_m() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_rate() )
      # self.t_l.append(self.env.now)
      q_state = self.q.state()
      # print("q_state= {}".format(pprint.pformat(q_state) ) )
      sub_qs_state = [max(q_state) for i in range(len(q_state) ) ]
      num_jobs = max(sub_qs_state)
      for i, s in enumerate(sub_qs_state):
        sub_qs_state[i] -= q_state[i]
      
      # num_job = max(state)
      # self.n_l.append(num_job)
      # rel_state = list_to_str(sub_qs_state)
      rel_state = "{},({},{})".format(num_jobs, sub_qs_state[1], sub_qs_state[2] )
      if rel_state not in self.polled_state__counter_map:
        self.polled_state__counter_map[rel_state] = 0
      self.polled_state__counter_map[rel_state] += 1
  
  def run_m(self):
    while True:
      p = (yield self.store.get() )
      if p.event_str == MPACKET_JOB_DEPARTED:
        # state = list_to_str(self.q.join_q.state() )
        state = list_to_str(self.q.state([p._id] ) )
        if state not in self.state__num_found_by_job_departed_map:
          self.state__num_found_by_job_departed_map[state] = 0
        self.state__num_found_by_job_departed_map[state] += 1
        
        state = self.q.state([p._id] )
        start_setup = "partial"
        if len(set(state) ) <= 1:
          start_setup = "complete"
        if start_setup not in self.start_setup__num_found_by_job_departed_map:
          self.start_setup__num_found_by_job_departed_map[start_setup] = 0
        self.start_setup__num_found_by_job_departed_map[start_setup] += 1
      else:
        log(ERROR, "unexpected p.event_str= {}".format(p.event_str) )
        return 1
  
  def put_m(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    return self.store.put(p)
  