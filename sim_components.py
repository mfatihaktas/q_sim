import simpy, random, copy, pprint

from patch import *

# *******************************  Packet  ****************************** #
SYS_TRAFF_SYM = 'r'
MDS_TRAFF_SYM = 'm'
SIMPLEX_TRAFF_SYM = 's'

class Packet(object):
  def __init__(self, time, _id, size=1, sym=None, flow_id=0):
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
    p = Packet(time=self.time, _id=self._id, size=self.size, sym=self.sym, flow_id=self.flow_id)
    p.ref_time = self.ref_time
    p.prev_hop_id = self.prev_hop_id
    p.entrance_time = self.entrance_time
    p.job_id = self.job_id
    p.winner_id = self.winner_id
    return p
  
  def __repr__(self):
    return "Packet[flow_id: {}, _id: {}]".format(self.flow_id, self._id)

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

# *******************************  PacketGenerator  ****************************** #
class PG(object): # Packet Generator
  def __init__(self, env, _id, inter_arr_dist, flow_id=0):
    self._id = _id
    self.env = env
    self.inter_arr_dist = inter_arr_dist
    self.flow_id = flow_id
    
    self.n_sent = 0
    self.out = None
    self.action = None
  
  def init(self):
    self.action = self.env.process(self.run() )  # starts the run() method as a SimPy process
  
  def run(self):
    while 1:
      # wait for next transmission
      yield self.env.timeout(self.inter_arr_dist() )
      self.n_sent += 1
      p = Packet(time=self.env.now, _id=self.n_sent, size=1, flow_id=self.flow_id)
      self.out.put(p)

# *************************************  JSink, JQ  ************************************* #
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

#  Any k out of n signals job completion
class JQ(object): # JoinQ
  def __init__(self, _id, env, k, input_qid_l):
    self._id = _id
    self.env = env
    self.k = k
    self.input_qid_l = input_qid_l
    
    self.job_id__p_l_map = {}
    self.qt_l = []
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.out = None
    self.out_c = None
    self.out_m = None
    self.action = env.process(self.run() )
    self.action = env.process(self.run_c() )
  
  def __repr__(self):
    return "JQ[_id= {}, k= {}, input_qid_l= {}]".format(self._id, self.k, self.input_qid_l)
  
  def state(self):
    return None
  
  def check_for_job_completion(self, p):
    p_l = self.job_id__p_l_map[p.job_id]
    if len(p_l) > self.k:
      log(ERROR, "len(p_l)= {} > k= {}".format(len(p_l), self.k) )
      return 1
    elif len(p_l) < self.k:
      return 0
    self.job_id__p_l_map.pop(p.job_id, None)
    self.qt_l.append(self.env.now - p.ref_time)
    recved_from_qid_l = [p.prev_hop_id for p in p_l]
    
    self.out_c.put_c(CPacket(_id=p.job_id, prev_hop_id=self._id, departed_qid_l=recved_from_qid_l) )
    p.winner_id = p.prev_hop_id
    p.prev_hop_id = self._id
    self.out.put(p)
    if self.out_m is not None:
      self.out_m.put_m(MPacket(_id=p.job_id, event_str=MPACKET_JOB_DEPARTED) )
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      if p.job_id not in self.job_id__p_l_map:
        self.job_id__p_l_map[p.job_id] = []
      self.job_id__p_l_map[p.job_id].append(p.deep_copy() )
      self.check_for_job_completion(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.ref_time = self.env.now
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      self.job_id__p_l_map.pop(cp._id, None)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

# *****************************************  Q  ************************************* #
class Q(object):
  def __init__(self, _id, env, num_s):
    self._id = _id
    self.env = env
    self.num_s = num_s

class FCFS(Q): # First Come First Serve
  def __init__(self, _id, env, serv_rate, serv_dist=None, rate=None, q_limit=None, debug=False):
    super().__init__(_id, env, 1)
    
    self.serv_rate = serv_rate
    self.serv_dist = lambda: random.expovariate(serv_rate)
    self.rate = rate
    self.q_limit = q_limit
    self.debug = debug
    
    self.p_l = []
    self.p_in_serv = None
    self.cancel, self.preempt = None, None
    self.n_recved = 0
    self.n_dropped = 0
    self.size_n = 0  # Current size of the queue in n
    self.busy = False  # Used to track if a packet is currently being sent
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
    return "FCFS[_id= {}, mu= {}]".format(self._id, self.serv_rate)
  
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
  
  def _in(self, job_id):
    if job_id == self.p_in_serv.job_id:
      return True
    for p in self.p_l:
      if job_id == p.job_id:
        return True
    return False
  
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
      self.busy = True
      if self.cancel is None:
        self.cancel = self.env.event()
      if self.preempt is None:
        self.preempt = self.env.event()
      
      exp_clock_start_time = None
      if self.serv_dist is None:
        if self.rate is None:
          log(ERROR, "self.serv_dist is None but self.rate is None too!")
          return 1
        log(WARNING, "self.serv_dist is None!")
        sim_log(WARNING, self.env, self, "starting D-clock! on ", self.p_in_serv)
        yield (self.cancel | self.preempt | self.env.timeout(self.p_in_serv.size/self.rate) ) # service
      else:
        exp_clock_start_time = self.env.now
        sim_log(DEBUG, self.env, self, "starting Exp-clock on ", self.p_in_serv)
        yield (self.cancel | self.preempt | self.env.timeout(self.serv_dist() ) ) # service
      if self.cancel is None: # task got cancelled
        sim_log(DEBUG, self.env, self, "cancelling", self.p_in_serv)
        sim_log(DEBUG, self.env, self, "cancelled Exp-clock on ", self.p_in_serv)
      elif self.preempt is None: # task got preempted
        self.p_l.insert(1, self.p_in_serv)
        sim_log(DEBUG, self.env, self, "preempted ", self.p_in_serv)
      else:
        sim_log(DEBUG, self.env, self, "done with Exp-clock in {}s on ".format(self.env.now-exp_clock_start_time), self.p_in_serv)
        self.qt_l.append(self.env.now - self.p_in_serv.ref_time)
        if self.out is not None and self.p_in_serv.sym != SYS_TRAFF_SYM:
          sim_log(DEBUG, self.env, self, "finished serv, forwarding", self.p_in_serv)
          self.p_in_serv.prev_hop_id = self._id
          self.out.put(self.p_in_serv)
        else:
          sim_log(DEBUG, self.env, self, "finished serv", self.p_in_serv)
      self.busy = False
  
  def put(self, p, preempt=False):
    self.n_recved += 1
    p.ref_time = self.env.now
    sim_log(DEBUG, self.env, self, "recved", p)
    t_size_n = self.size_n + 1
    if (self.q_limit is not None and t_size_n > self.q_limit):
      sim_log(DEBUG, self.env, self, "dropping", p)
      self.n_dropped += 1
      return
    else:
      self.size_n = t_size_n
      if preempt and self.busy:
        self.preempt.succeed()
        self.preempt = None
        self.p_l.insert(0, p)
        self.syncer.put(1)
        return
      return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      if self.p_in_serv is not None and self.p_in_serv.job_id == cp._id:
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
