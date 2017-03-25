import simpy, random, copy, pprint

from patch import *

# *******************************  Packet  ****************************** #
SYS_TRAFF_SYM = 'r'
MDS_TRAFF_SYM = 'm'
SIMPLEX_TRAFF_SYM = 's'

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
    return "Packet[_id: {}, prev_hop_id: {}, job_id: {}, sym: {}]".\
      format(self._id, self.prev_hop_id, self.job_id, self.sym)

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
class PacketGenerator(object):
  def __init__(self, env, _id, adist, sdist = None, initial_delay=0, finish=float("inf"), flow_id=0, sym=None):
    self._id = _id
    self.env = env
    self.adist = adist
    self.sdist = sdist
    self.initial_delay = initial_delay
    self.finish = finish
    self.flow_id = flow_id
    self.sym = sym
    
    self.n_sent = 0
    self.out = None
    self.action = None
  
  def init(self):
    self.action = self.env.process(self.run() )  # starts the run() method as a SimPy process
  
  def run(self):
    yield self.env.timeout(self.initial_delay)
    while self.env.now < self.finish:
      # wait for next transmission
      yield self.env.timeout(self.adist() )
      self.n_sent += 1
      p = Packet(time=self.env.now, size=1, _id=self.n_sent, flow_id=self.flow_id, sym=self.sym)
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

class JQ(object): # JoinQ for MDS; completion of any k tasks out of n means job completion
  def __init__(self, _id, env, k, input_qid_l):
    self._id = _id
    self.env = env
    self.k = k
    self.input_qid_l = input_qid_l
    
    self.fj = False
    if self.k == len(input_qid_l):
      self.fj = True
    # self.input_id__pq_map = {i: [] for i in input_qid_l}
    
    self.job_id__p_l_map = {}
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
    return "JQ[_id= {}, k= {}, input_qid_l= {}]".format(self._id, self.k, self.input_qid_l)
  
  def state(self):
    # return [len(pq) for i, pq in self.input_id__pq_map.items() ]
    return None
  
  def check_for_job_completion(self, p):
    # now = self.env.now
    # len_l = [len(pq) for i, pq in self.input_id__pq_map.items() ]
    # num_non_zero = len([l for l in len_l if l > 0] )
    # if num_non_zero > self.k:
    #   log(ERROR, "num_non_zero= {} > k= {}".format(num_non_zero, self.k) )
    #   return 1
    # elif num_non_zero < self.k:
    #   return 0
    
    # ref_p = None
    # departed_q_id_l = []
    # for j, pq in self.input_id__pq_map.items():
    #   if len(pq) == 0:
    #     continue
    #   p = pq.pop(0)
    #   departed_q_id_l.append(p.prev_hop_id)
    #   if ref_p is None:
    #     ref_p = p
    #   else:
    #     if (p.prev_hop_id == ref_p.prev_hop_id) or \
    #       (p.entrance_time != ref_p.entrance_time) or \
    #       (p.job_id != ref_p.job_id):
    #       log(ERROR, "JQ= {}\n\tsupposed to be tasks of the same job;\n\tp= {}\n\tref_p= {}".format(self, p, ref_p) )
    #       return 1
    #   self.qt_l.append(now - p.ref_time)
    # if not self.fj:
    #   self.out_c.put_c(CPacket(_id=ref_p.job_id, prev_hop_id=self._id, departed_q_id_l=departed_q_id_l) )
    
    # ref_p.winner_id = ref_p.prev_hop_id
    # ref_p.prev_hop_id = self._id
    # self.out.put(ref_p)
    # self.length = max([len(pq) for i, pq in self.input_id__pq_map.items() ] )
    # if self.out_m is not None:
    #   self.out_m.put_m(MPacket(_id=ref_p.job_id, event_str=MPACKET_JOB_DEPARTED) )
    
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
      # if p.prev_hop_id not in self.input_qid_l:
      #   log(ERROR, "packet can NOT continue {}; packet= {}".format(self, p) )
      #   return 1
      
      # state = list_to_str(state() )
      # if state not in self.state__num_found_map:
      #   self.state__num_found_map[state] = 0
      # self.state__num_found_map[state] += 1
      
      # self.input_id__pq_map[p.prev_hop_id].append(p)
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
      # for j, pq in self.input_id__pq_map.items():
      #   for p in pq:
      #     if p.job_id == cp._id:
      #       pq.remove(p)
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
    self.cancel, self.preempt = None, None
    self.n_recved = 0
    self.n_dropped = 0
    self.size_n = 0  # Current size of the queue in n
    self.size_B = 0  # Current size of the queue in bytes
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
    t_size_B = self.size_B + p.size
    if (self.qlimit_n is not None and t_size_n > self.qlimit_n) or \
       (self.qlimit_B is not None and t_size_B > self.qlimit_B):
      sim_log(DEBUG, self.env, self, "dropping", p)
      self.n_dropped += 1
      return
    else:
      self.size_n = t_size_n
      self.size_B = t_size_B
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
