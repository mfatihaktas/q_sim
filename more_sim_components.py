import simpy, random, copy, pprint
from simpy.core import BoundClass
from simpy.resources import base
from heapq import heappush, heappop
from patch import *
from deprecated import *

class Packet(object):
  def __init__(self, time, size, _id, sym=None, flow_id=0):
    self.time = time
    self.size = size
    self._id = _id
    self.sym = sym
    self.flow_id = flow_id
    self.ref_time = 0 # for casual use
    # for FJ and MDS Q implementation
    self.prev_hop_id = None
    self.entrance_time = None
    self.job_id = None
    self.winner_id = None
  
  def deep_copy(self):
    p = Packet(time=self.time, size=self.size, _id=self._id, sym=self.sym, flow_id=self.flow_id)
    p.ref_time = self.ref_time
    p.prev_hop_id = self.prev_hop_id
    p.entrance_time = self.entrance_time
    p.job_id = self.job_id
    p.winner_id = self.winner_id
    return p
  
  def __repr__(self):
    return "Packet[_id: {}, prev_hop_id: {}, job_id: {}]".\
      format(self._id, self.prev_hop_id, self.job_id)

class CPacket(object): # Control
  def __init__(self, _id, prev_hop_id=None, departed_qid_l=[]):
    self._id = _id
    self.prev_hop_id = prev_hop_id
    self.departed_qid_l = departed_qid_l
  
  def deep_copy(self):
   return CPacket(self._id, self.prev_hop_id)
  
  def __repr__(self):
    return "CPacket[_id= {}, prev_hop_id= {}, departed_qid_l= {}]".format(self._id, self.prev_hop_id, self.departed_qid_l)

MPACKET_JOB_DEPARTED = "job_departed"
class MPacket(object): # Monitor
  def __init__(self, _id, event_str):
    self._id = _id
    self.event_str = event_str
  
  def deep_copy(self):
   return MPacket(self._id, self.event_str)
  
  def __repr__(self):
    return "MPacket[_id= {}, event_str= {}]".format(self._id, self.event_str)

class MT_PacketGenerator(object):
  def __init__(self, env, _id, adist, sdist, sym_l=None, initial_delay=0, finish=float("inf"), flow_id=0):
    self._id = _id
    self.env = env
    self.adist = adist
    self.sdist = sdist
    self.sym_l = sym_l
    self.initial_delay = initial_delay
    self.finish = finish
    self.n_sent = 0
    self.flow_id = flow_id
    
    self.sym__n_sent = {}
    self.out = None
    self.action = env.process(self.run())  # starts the run() method as a SimPy process

  def run(self):
    yield self.env.timeout(self.initial_delay)
    while self.env.now < self.finish:
      # wait for next transmission
      yield self.env.timeout(self.adist() )
      self.n_sent += 1
      if self.sym_l is not None:
        sym = self.sym_l[random.randint(0,len(self.sym_l)-1) ]
        if sym not in self.sym__n_sent:
          self.sym__n_sent[sym] = 0
        self.sym__n_sent[sym] = self.sym__n_sent[sym] + 1
        p = Packet(time=self.env.now, size=self.sdist(), _id=self.n_sent, sym=sym, flow_id=self.flow_id)
      else:
        p = Packet(time=self.env.now, size=self.sdist(), _id=self.n_sent, flow_id=self.flow_id)
      self.out.put(p)

class Q(object):
  def __init__(self, _id, env, num_s):
    self._id = _id
    self.env = env
    self.num_s = num_s

class S1_Q(Q): # Memoryless service, 1 server
  def __init__(self, _id, env, serv_rate, serv_dist=None, rate=None, qlimit_n=None, qlimit_B=None, debug=False):
    super().__init__(_id, env, 1)
    # self._id = _id
    # self.env = env
    self.serv_rate = serv_rate
    self.serv_dist = lambda: random.expovariate(serv_rate)
    self.rate = rate
    self.qlimit_n = qlimit_n
    self.qlimit_B = qlimit_B
    self.debug = debug
    
    self.p_l = []
    self.p_in_serv = None
    self.cancel = None
    self.n_recved = 0
    self.n_dropped = 0
    self.size_n = 0  # Current size of the queue in n
    self.size_B = 0  # Current size of the queue in bytes
    self.busy = 0  # Used to track if a packet is currently being sent
    self.wt_l = []
    self.qt_l = []
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.syncer = simpy.Store(env) # simpy.Resource(env, capacity=1)
    self.out = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    self.action = env.process(self.run_helper() )
  
  def __repr__(self):
    # return "S1_Q[_id= {}, rate= {}, qlimit_n= {}, qlimit_B= {}]".format(_id, self.rate, self.qlimit_n, self.qlimit_B)
    return "S1_Q[_id= {}, mu= {}]".format(self._id, self.serv_rate)
  
  def length(self):
    return len(self.p_l) + self.busy
  
  def state(self, job_id_to_exclude=[]):
    if not job_id_to_exclude:
      return [self.length() ]
    
    p_l = []
    for p in self.p_l:
      if p.job_id not in job_id_to_exclude:
        p_l.append(p)
    state = len(p_l)
    if self.p_in_serv != None and (self.p_in_serv.job_id not in job_id_to_exclude):
      state += 1
    return [state]
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      self.p_l.append(p)
      self.syncer.put(1)
  
  def run_helper(self): # To implement del from self.store
    while True:
      (yield self.syncer.get() )
      if len(self.p_l) == 0:
        continue # log(ERROR, "self.p_l is empty!") # May happen because of task cancellations
      self.p_in_serv = self.p_l.pop(0)
      self.wt_l.append(self.env.now - self.p_in_serv.ref_time)
      self.size_n -= 1
      self.size_B -= self.p_in_serv.size
      self.busy = 1
      if self.cancel is None:
        self.cancel = self.env.event()
      exp_clock_start_time = None
      if self.serv_dist is None:
        if self.rate is None:
          log(ERROR, "self.serv_dist is None but self.rate is None too!")
          return 1
        log(WARNING, "self.serv_dist is None!")
        sim_log(WARNING, self.env, self, "starting D-clock! on ", self.p_in_serv)
        yield (self.cancel | self.env.timeout(self.p_in_serv.size/self.rate) ) # service
      else:
        exp_clock_start_time = self.env.now
        sim_log(DEBUG, self.env, self, "starting Exp-clock on ", self.p_in_serv)
        yield (self.cancel | self.env.timeout(self.serv_dist() ) ) # service
      if self.cancel is None: # task got cancelled
        sim_log(DEBUG, self.env, self, "cancelling", self.p_in_serv)
        sim_log(DEBUG, self.env, self, "cancelled Exp-clock on ", self.p_in_serv)
      else:
        sim_log(DEBUG, self.env, self, "done with Exp-clock in {}s on ".format(self.env.now-exp_clock_start_time), self.p_in_serv)
        self.qt_l.append(self.env.now - self.p_in_serv.ref_time)
        if self.out is not None:
          sim_log(DEBUG, self.env, self, "finished serv, forwarding", self.p_in_serv)
          self.p_in_serv.prev_hop_id = self._id
          self.out.put(self.p_in_serv)
        else:
          sim_log(DEBUG, self.env, self, "finished serv", self.p_in_serv)
          
      self.busy = 0
  
  def put(self, p):
    self.n_recved += 1
    p.ref_time = self.env.now
    sim_log(DEBUG, self.env, self, "recved", p)
    t_size_n = self.size_n + 1
    t_size_B = self.size_B + p.size
    if (self.qlimit_n is not None and t_size_n > self.qlimit_n) or \
       (self.qlimit_B is not None and t_size_B > self.qlimit_B):
      sim_log(DEBUG, self.env, self, "dropping", p)
      self.n_dropped += 1
      return
    else:
      self.size_n = t_size_n
      self.size_B = t_size_B
      return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      if self.p_in_serv.job_id == cp._id:
        if self.cancel is None:
          log(ERROR, "self.cancel is None!")
          return 1
        self.cancel.succeed()
        self.cancel = None
      
      for p in self.p_l:
        if p.job_id == cp._id:
          self.p_l.remove(p)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)
    
class QMonitor(object):
  def __init__(self, env, q, dist):
    self.q = q
    self.env = env
    self.dist = dist
    
    self.t_l = [] # time steps that the numbers polled from the q
    self.n_l = [] # num of packets in the q
    self.action = env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.dist() )
      
      self.t_l.append(self.env.now)
      self.n_l.append(len(self.q.store.items) + self.q.busy)

# *******************************************  MDSQ  ********************************************* #
class JSink(object): # Join
  def __init__(self, _id, env):
    self._id = _id
    self.env = env
    
    self.st_l = [] # total time in system
    self.qid__num_win_map = {}
    
    self.store = simpy.Store(env)
    self.out = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
  
  def __repr__(self):
    return "JSink[_id= {}]".format(self._id)
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      self.st_l.append(self.env.now - p.entrance_time)
      if p.winner_id not in self.qid__num_win_map:
        self.qid__num_win_map[p.winner_id] = 0
      self.qid__num_win_map[p.winner_id] += 1
      
      if self.out is not None:
        self.out.put(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    return self.store.put(p)

class JQ(object):
  def __init__(self, _id, env, input_qid_l, sym__rgroup_l_map):
    self._id = _id
    self.env = env
    self.input_qid_l = input_qid_l
    self.sym__rgroup_l_map = sym__rgroup_l_map
    
    # self.input_id__pq_map = {i: [] for i in input_qid_l}
    self.job_id__p_l_map = {}
    # self.job_id__departed_qid_l_map = {}
    self.qt_l = []
    self.length = 0 # maximum of the lengths of all pq's
    # self.state__num_found_map = {}
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.out = None
    self.out_c = None
    self.out_m = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
  
  def __repr__(self):
    return "JQ[_id= {}, input_qid_l= {}]".format(self._id, self.input_qid_l)
  
  def state(self):
    # return [len(pq) for i, pq in self.input_id__pq_map.items() ]
    return None
  
  def check_for_job_completion(self, p):
    p_l = self.job_id__p_l_map[p.job_id]
    recved_from_qid_l = [p.prev_hop_id for p in p_l]
    rgroup_l = self.sym__rgroup_l_map[p.sym]
    success = False
    for rgroup in rgroup_l:
      success = True
      for qid in rgroup:
        if qid not in recved_from_qid_l:
          success = False
      if success:
        break
    if success:
      self.out_c.put_c(CPacket(_id=p.job_id, prev_hop_id=self._id, departed_qid_l=[p.prev_hop_id for p in p_l] ) )
      
      p.winner_id = p.prev_hop_id
      p.prev_hop_id = self._id
      self.out.put(p)
      if self.out_m is not None:
        self.out_m.put_m(MPacket(_id=p.job_id, event_str=MPACKET_JOB_DEPARTED) )
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      # self.input_id__pq_map[p.prev_hop_id].append(p)
      if p.job_id not in self.job_id__p_l_map:
        self.job_id__p_l_map[p.job_id] = []
        # self.job_id__departed_qid_l_map[p.job_id] = []
      self.job_id__p_l_map[p.job_id].append(p)
      # self.job_id__departed_qid_l_map[p.job_id].append([p.prev_hop_id] )
      self.check_for_job_completion(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.ref_time = self.env.now
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      for j, pq in self.input_id__pq_map.items():
        for p in pq:
          if p.job_id == cp._id:
            pq.remove(p)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

# **************************************  Availability Q  **************************************** #
class CodedStorageQ(object):
  def __init__(self, _id, env, qid_l, qmu_l, sym__rgroup_l_map, w_sys=False, out=None):
    self._id = _id
    self.env = env
    # self.out = out
    
    self.num_q = len(qid_l)
    
    self.join_sink = JSink(_id, env)
    self.join_sink.out = out
    self.qid_q_map = {}
    
    self.qid_l = qid_l
    self.join_q = JQ(_id=_id, env=env, input_qid_l=self.qid_l, sym__rgroup_l_map=sym__rgroup_l_map)
    self.join_q.out = self.join_sink
    self.join_q.out_c = self
    self.join_q.out_m = None # can be set by the caller if desired
    for i, qid in enumerate(qid_l):
      q = S1_Q(_id=qid, env=env, serv_rate=qmu_l[i] )
      log(DEBUG, "i= {}, q= {}".format(i, q) )
      q.out = self.join_q
      self.qid_q_map[i] = q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.out = None
    self.action = env.process(self.run() )
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
  
  def __repr__(self):
    return "CSQ[qid_l= {}]".format(self.qid_l)
  
  def state(self, job_id_to_exclude=[]):
    state = []
    for i, q in self.qid_q_map.items():
      state += q.state(job_id_to_exclude)
    return state
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      for qid, q in self.qid_q_map.items():
        q.put(p.deep_copy() )
      
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.job_id_counter += 1
    p.job_id = self.job_id_counter
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      for g, q in self.qid_q_map.items():
        if g != -1 and q._id not in cp.departed_qid_l:
          q.put_c(cp.deep_copy() )
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

class AVQMonitor(object):
  def __init__(self, env, aq, poll_dist):
    self.aq = aq
    self.env = env
    self.poll_dist = poll_dist
    
    self.t_l = [] # Time steps that the numbers polled from the aq
    # self.n_l = [] # Num of jobs in the aq
    self.polled_state__counter_map = {}
    self.state__num_found_by_job_departed_map = {}
    self.start_setup__num_found_by_job_departed_map = {}
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    self.action = env.process(self.run_m() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_dist() )
      # self.t_l.append(self.env.now)
      aq_state = self.aq.state()
      # print("aq_state= {}".format(pprint.pformat(aq_state) ) )
      sub_qs_state = [max(aq_state) for i in range(len(aq_state) ) ]
      num_jobs = max(sub_qs_state)
      for i, s in enumerate(sub_qs_state):
        sub_qs_state[i] -= aq_state[i]
      
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
        # state = list_to_str(self.aq.join_q.state() )
        state = list_to_str(self.aq.state([p._id] ) )
        if state not in self.state__num_found_by_job_departed_map:
          self.state__num_found_by_job_departed_map[state] = 0
        self.state__num_found_by_job_departed_map[state] += 1
        
        state = self.aq.state([p._id] )
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
  