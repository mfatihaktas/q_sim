import simpy, random, copy, pprint
from simpy.core import BoundClass
from simpy.resources import base
from heapq import heappush, heappop
from patch import *
from deprecated import *

class Packet(object):
  def __init__(self, time, size, _id, src="a", dst="z", flow_id=0):
    self.time = time
    self.size = size
    self._id = _id
    self.src = src
    self.dst = dst
    self.flow_id = flow_id
    self.ref_time = 0 # for casual use
    # for FJ and MDS Q implementation
    self.prev_hop_id = None
    self.entrance_time = None
    self.job_id = None
  
  def deep_copy(self):
    p = Packet(time=self.time, size=self.size, _id=self._id, src=self.src, dst=self.dst, flow_id=self.flow_id)
    p.ref_time = self.ref_time
    p.prev_hop_id = self.prev_hop_id
    p.entrance_time = self.entrance_time
    p.job_id = self.job_id
    return p
  
  def __repr__(self):
    return "Packet[_id: {}, prev_hop_id: {}, job_id: {}]".\
      format(self._id, self.prev_hop_id, self.job_id)

class CPacket(object): # Control
  def __init__(self, _id, prev_hop_id=None):
    self._id = _id
    self.prev_hop_id = prev_hop_id
  
  def deep_copy(self):
   return CPacket(self._id, self.prev_hop_id)
  
  def __repr__(self):
    return "CPacket[_id= {}, prev_hop_id= {}]".format(self._id, self.prev_hop_id)

MPACKET_JOB_DEPARTED = "job_departed"
class MPacket(object): # Monitor
  def __init__(self, _id, event_str):
    self._id = _id
    self.event_str = event_str
  
  def deep_copy(self):
   return MPacket(self._id, self.event_str)
  
  def __repr__(self):
    return "MPacket[_id= {}, event_str= {}]".format(self._id, self.event_str)

class PacketGenerator(object):
  def __init__(self, env, _id, adist, sdist, initial_delay=0, finish=float("inf"), flow_id=0):
    self._id = _id
    self.env = env
    self.adist = adist
    self.sdist = sdist
    self.initial_delay = initial_delay
    self.finish = finish
    self.n_sent = 0
    self.flow_id = flow_id
    
    self.out = None
    self.action = env.process(self.run())  # starts the run() method as a SimPy process

  def run(self):
    yield self.env.timeout(self.initial_delay)
    while self.env.now < self.finish:
      # wait for next transmission
      yield self.env.timeout(self.adist() )
      self.n_sent += 1
      p = Packet(time=self.env.now, size=self.sdist(), _id=self.n_sent, src=self._id, flow_id=self.flow_id)
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
    
    self.p_list = []
    self.p_in_serv = None
    self.cancel = None
    self.n_recved = 0
    self.n_dropped = 0
    self.size_n = 0  # Current size of the queue in n
    self.size_B = 0  # Current size of the queue in bytes
    self.busy = 0  # Used to track if a packet is currently being sent
    self.wt_list = []
    self.qt_list = []
    
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
    return len(self.p_list) + self.busy
  
  def state(self, job_id_to_exclude=[]):
    if not job_id_to_exclude:
      return [self.length() ]
    
    p_list = []
    for p in self.p_list:
      if p.job_id not in job_id_to_exclude:
        p_list.append(p)
    state = len(p_list)
    if self.p_in_serv != None and (self.p_in_serv.job_id not in job_id_to_exclude):
      state += 1
    return [state]
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      self.p_list.append(p)
      self.syncer.put(1)
  
  def run_helper(self): # To implement del from self.store
    while True:
      (yield self.syncer.get() )
      if len(self.p_list) == 0:
        continue # log(ERROR, "self.p_list is empty!") # May happen because of task cancellations
      self.p_in_serv = self.p_list.pop(0)
      self.wt_list.append(self.env.now - self.p_in_serv.ref_time)
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
        self.qt_list.append(self.env.now - self.p_in_serv.ref_time)
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
      
      for p in self.p_list:
        if p.job_id == cp._id:
          self.p_list.remove(p)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)
    
class QMonitor(object):
  def __init__(self, env, q, dist):
    self.q = q
    self.env = env
    self.dist = dist
    
    self.t_list = [] # time steps that the numbers polled from the q
    self.n_list = [] # num of packets in the q
    self.action = env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.dist() )
      
      self.t_list.append(self.env.now)
      self.n_list.append(len(self.q.store.items) + self.q.busy)

# *******************************************  MDSQ  ********************************************* #
class JSink(object): # Join
  def __init__(self, _id, env):
    self._id = _id
    self.env = env
    
    self.st_list = [] # total time in system
    self.qid__num_win_map = {}
    
    self.store = simpy.Store(env)
    self.out = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
  
  def __repr__(self):
    return "JSink[_id= {}]".format(self._id)
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      self.st_list.append(self.env.now - p.entrance_time)
      if p.prev_hop_id not in self.qid__num_win_map:
        self.qid__num_win_map[p.prev_hop_id] = 0
      self.qid__num_win_map[p.prev_hop_id] += 1
      
      if self.out is not None:
        self.out.put(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    return self.store.put(p)

class JQ(object): # JoinQ for MDS; completion of any k tasks out of n means job completion
  def __init__(self, _id, env, k, input_qid_list):
    self._id = _id
    self.env = env
    self.k = k
    self.input_qid_list = input_qid_list
    
    self.fj = False
    if self.k == len(input_qid_list):
      self.fj = True
    self.input_id__pq_map = {i: [] for i in input_qid_list}
    self.qt_list = []
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
    return "JQ[_id= {}, k= {}, input_qid_list= [{}]]".format(self._id, self.k, ",".join(self.input_qid_list) )
  
  def state(self):
    return [len(pq) for i, pq in self.input_id__pq_map.items() ]
  
  def check_for_job_completion(self, possible_winner_id):
    now = self.env.now
    len_list = [len(pq) for i, pq in self.input_id__pq_map.items() ]
    num_non_zero = len([l for l in len_list if l > 0] )
    if num_non_zero > self.k:
      log(ERROR, "num_non_zero= {} > k= {}".format(num_non_zero, self.k) )
      return 1
    elif num_non_zero < self.k:
      return 0
    
    ref_p = None
    for j, pq in self.input_id__pq_map.items():
      if len(pq) == 0:
        continue
      p = pq.pop(0)
      if ref_p is None:
        ref_p = p
      else:
        if (p.prev_hop_id == ref_p.prev_hop_id) or \
           (p.entrance_time != ref_p.entrance_time) or \
           (p.job_id != ref_p.job_id):
          log(ERROR, "supposed to be tasks of the same job;\n\tp= {}\n\tref_p= {}".format(p, ref_p) )
          return 1
      self.qt_list.append(now - p.ref_time)
    if not self.fj:
      self.out_c.put_c(CPacket(_id=ref_p.job_id, prev_hop_id=self._id) )
    
    ref_p.prev_hop_id = possible_winner_id
    self.out.put(ref_p)
    self.length = max([len(pq) for i, pq in self.input_id__pq_map.items() ] )
    if self.out_m is not None:
      self.out_m.put_m(MPacket(_id=ref_p.job_id, event_str=MPACKET_JOB_DEPARTED) )
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      if p.prev_hop_id not in self.input_qid_list:
        log(ERROR, "packet can NOT continue {}; packet= {}".format(self, p) )
        return 1
      # state = list_to_str(state() )
      # if state not in self.state__num_found_map:
      #   self.state__num_found_map[state] = 0
      # self.state__num_found_map[state] += 1
      
      self.input_id__pq_map[p.prev_hop_id].append(p)
      self.check_for_job_completion(p.prev_hop_id)
  
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

"""
Fork incoming job to all sub-q's, wait for any k task to complete, cancel the remainings.
"""
class MDSQ(object):
  def __init__(self, _id, env, k, qid_list, qserv_rate_list, out=None):
    self._id = _id
    self.env = env
    self.k = k
    self.qid_list = qid_list
    self.qserv_rate_list = qserv_rate_list
    # self.out = out
    
    self.num_q = len(qid_list)
    self.join_sink = JSink(_id, env)
    self.join_sink.out = out
    self.join_q = JQ(_id, env, k, qid_list)
    self.join_q.out = self.join_sink
    self.join_q.out_c = self
    self.id_q_map = {}
    for i in range(self.num_q):
      qid = qid_list[i]
      m1_q = S1_Q(_id=qid, env=env, serv_rate=qserv_rate_list[i] )
      m1_q.out = self.join_q
      self.id_q_map[qid] = m1_q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.action = env.process(self.run() ) # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
  
  def __repr__(self):
    return "MDSQ[k= {}, qid_list= [{}] ]".format(self.k, ",".join(self.qid_list) )
  
  def length(self):
    return max([q.length() for i, q in self.id_q_map.items() ] )
  
  def state(self, job_id_to_exclude=[]):
    # return [q.length() for i, q in self.id_q_map.items() ]
    state = []
    for i, q in self.id_q_map.items():
      state += q.state(job_id_to_exclude)
    return state
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      for i, q in self.id_q_map.items():
        q.put(p.deep_copy() )
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    if p.entrance_time is None:
      p.entrance_time = self.env.now
    if p.job_id is None:
      self.job_id_counter += 1
      p.job_id = self.job_id_counter
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      for i, q in self.id_q_map.items():
        q.put_c(cp)
      if cp.prev_hop_id != self._id: # To avoid inifinite loops in case cp came from a different jq then self.join_q
        self.join_q.put_c(cp)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

class MDSQMonitor(object):
  def __init__(self, env, q, poll_dist):
    self.q = q
    self.env = env
    self.poll_dist = poll_dist
    
    self.t_list = [] # Time steps that the numbers polled from the aq
    # self.n_list = [] # Num of packets in the aq
    self.polled_state__counter_map = {}
    
    self.action = env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_dist() )
      
      self.t_list.append(self.env.now)
      state = self.q.state()
      # print("state= {}".format(list_to_str(state) ) )
      num_job = max(state)
      # self.n_list.append(num_job)
      # rel_state = ",".join("%d" % (num_job-l) for l in state)
      rel_state = list_to_str(state)
      if rel_state not in self.polled_state__counter_map:
        self.polled_state__counter_map[rel_state] = 0
      self.polled_state__counter_map[rel_state] += 1
# **************************************  Availability Q  **************************************** #
"""
  Implements (n, k, r, t)-LRC Q
  A file consists of k-pieces and mapped to n code-words with an MDS code.
  The (r, t)-availability ensures that each systematic node can be regenerated using one of the t 
  disjoint repair groups of other storage nodes, each of size at most r (typically << k)
  Ref: When do the Availability Codes Make the Stored Data More Available?
"""
class AVQ(object): # Availability
  def __init__(self, _id, env, k, r, t, qid_list, qserv_rate_list, out=None):
    self._id = _id
    self.env = env
    self.k = k
    self.r = r
    self.t = t
    self.qid_list = qid_list
    # self.out = out
    
    self.num_q = len(qid_list)
    if self.num_q != (1 + t*r):
      log(ERROR, "self.num_q= {} != (1 + t*r)= {}".format(self.num_q, (1 + t*r) ) )
      return 1
    self.join_sink = JSink(_id, env)
    self.join_sink.out = out
    self.join_q = JQ(_id=_id, env=env, k=1, input_qid_list=qid_list)
    self.join_q.out = self.join_sink
    self.join_q.out_c = self
    self.join_q.out_m = None # can be set by the caller if desired
    self.group_id__q_map = {}
    for g in range(1, t + 1 + 1):
      q = None
      if g == 1:
        q = S1_Q(_id=qid_list[0], env=env, serv_rate=qserv_rate_list[0] )
        # plot_exp_dist("q1", lambda: random.expovariate(qserv_rate_list[0] ) )
        log(DEBUG, "g= {}, q= {}, qserv_rate_list[0]= {}".format(g, q, qserv_rate_list[0]) )
        q.out = self.join_q
      else:
        li = 1+(g-2)*r
        ri = li + r
        q = MDSQ(_id="".join(["%s," % i for i in qid_list[li:ri] ] ),
                 env=env, k=k, qid_list=qid_list[li:ri], qserv_rate_list=qserv_rate_list[li:ri], out=self.join_q)
        log(DEBUG, "g= {}, q= {}, qserv_rate_list[li:ri]= {}".format(g, q, qserv_rate_list[li:ri] ) )
      self.group_id__q_map[g] = q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.out = None
    self.action = env.process(self.run() ) # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
  
  def __repr__(self):
    return "AVQ[k= {}, r= {}, t= {}]".format(self.k, self.r, self.t)
    # return "AVQ[k= {}, r= {}, t= {}, qid_list= [{}] ]".format(self.k, self.r, ",".join(self.qid_list) )
  
  def state(self, job_id_to_exclude=[]):
    state = []
    for g, q in self.group_id__q_map.items():
      state += q.state(job_id_to_exclude)
    
    return state
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      for g, q in self.group_id__q_map.items():
        if g != -1:
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
      for g, q in self.group_id__q_map.items():
        if g != -1:
          q.put_c(cp.deep_copy() )
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

class AVQMonitor(object):
  def __init__(self, env, aq, poll_dist):
    self.aq = aq
    self.env = env
    self.poll_dist = poll_dist
    
    self.t_list = [] # Time steps that the numbers polled from the aq
    # self.n_list = [] # Num of jobs in the aq
    self.polled_state__counter_map = {}
    self.state__num_found_by_job_departed_map = {}
    self.start_setup__num_found_by_job_departed_map = {}
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    self.action = env.process(self.run_m() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_dist() )
      # self.t_list.append(self.env.now)
      aq_state = self.aq.state()
      # print("aq_state= {}".format(pprint.pformat(aq_state) ) )
      sub_qs_state = [max(aq_state) for i in range(len(aq_state) ) ]
      num_jobs = max(sub_qs_state)
      for i, s in enumerate(sub_qs_state):
        sub_qs_state[i] -= aq_state[i]
      
      # num_job = max(state)
      # self.n_list.append(num_job)
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
  