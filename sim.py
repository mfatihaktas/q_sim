import simpy, random, copy, pprint
from rvs import *
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
    # return "Packet[flow_id: {}, _id: {}]".format(self.flow_id, self._id)
    return "Packet[_id: {}, sym: {}]".format(self.job_id, self.sym)

class CPacket(object): # Control
  def __init__(self, _id, sym=None, prev_hop_id=None, departed_qid_l=[]):
    self._id = _id
    self.sym = sym
    self.prev_hop_id = prev_hop_id
    self.departed_qid_l = departed_qid_l
  
  def deep_copy(self):
   return CPacket(self._id, self.sym, self.prev_hop_id, list(self.departed_qid_l) )
  
  def __repr__(self):
    # return "CPacket[_id= {}, prev_hop_id= {}, departed_qid_l= {}]".format(self._id, self.prev_hop_id, self.departed_qid_l)
    return "CPacket[_id= {}, sym= {}, departed_qid_l= {}]".format(self._id, self.sym, self.departed_qid_l)

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
  def __init__(self, env, _id, ar, flow_id=0, psize=None, psize_dist_m=None):
    self._id = _id
    self.env = env
    self.ar = ar
    self.flow_id = flow_id
    
    if psize == "Pareto":
      self.psize_dist = Pareto(psize_dist_m['l'], psize_dist_m['a'] )
    else:
      self.psize_dist = None
    
    self.n_sent = 0
    self.out = None
  
  def init(self):
    self.env.process(self.run() )
  
  def run(self):
    while 1:
      yield self.env.timeout(random.expovariate(self.ar) )
      self.n_sent += 1
      
      psize = 1 if self.psize_dist is None else self.psize_dist.gen_sample()
      p = Packet(time=self.env.now, _id=self.n_sent, flow_id=self.flow_id, size=psize)
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
    self.action = env.process(self.run() )
  
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
  def __init__(self, _id, env):
    self._id = _id
    self.env = env

class FCFS(Q): # First Come First Serve
  def __init__(self, _id, env, sdist_m, out=None):
    super().__init__(_id, env)
    self.servdist = sdist_m['dist']
    self.stime = rv_from_m(sdist_m)
    self.out = out
    
    self.p_inserv = None
    self.p_l = []
    
    self.cancel, self.preempt = None, None
    self.cancel_flag, self.preempt_flag = False, False
    # self.wt_l, self.qt_l = [], []
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.syncer = simpy.Store(env) # simpy.Resource(env, capacity=1)
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
    self.action = env.process(self.run_helper() )
  
  def __repr__(self):
    return "FCFS[_id= {}, stime= {}]".format(self._id, self.stime)
  
  def busy(self):
    return self.length() > 0
  
  def length(self):
    return len(self.p_l) + (self.p_inserv is not None)
  
  # def pstate_l(self):
  #   l = ["{}, ji= {}".format(self.p_inserv.sym, self.p_inserv.job_id) ]
  #   for p in self.p_l:
  #     l.append("{}, ji= {}".format(p.sym, p.job_id) )
  #   return l
  
  def num_sym_in(self, sym):
    num = 0
    if (self.p_inserv is not None) and self.p_inserv.sym == sym:
      num += 1
    for p in self.p_l:
      if p.sym == sym:
        num += 1
    return num
  
  def state(self, job_id_to_exclude=[]):
    if not job_id_to_exclude:
      return [self.length() ]
    
    p_l = []
    for p in self.p_l:
      if p.job_id not in job_id_to_exclude:
        p_l.append(p)
    state = len(p_l)
    if self.p_inserv != None and (self.p_inserv.job_id not in job_id_to_exclude):
      state += 1
    return [state]
  
  def _in(self, job_id):
    if self.p_inserv is None:
      return False
    elif self.p_inserv.job_id == job_id:
      return True
    for p in self.p_l:
      if job_id == p.job_id:
        return True
    return False
  
  def run(self):
    while True:
      yield self.store.get()
      self.syncer.put(1)
  
  def run_helper(self): # To implement del from self.store
    while True:
      yield self.syncer.get()
      if len(self.p_l) == 0:
        continue # log(ERROR, "self.p_l is empty!") # May happen because of task cancellations
      self.p_inserv = self.p_l.pop(0)
      # self.wt_l.append(self.env.now - self.p_inserv.ref_time)
      
      self.cancel = self.env.event()
      self.preempt = self.env.event()
      
      clk_start_time = self.env.now
      t = self.p_inserv.size*self.stime.gen_sample()
      sim_log(DEBUG, self.env, self, "starting {}-clock! on ".format(self.servdist), self.p_inserv)
      yield (self.cancel | self.preempt | self.env.timeout(t) )
      
      if self.cancel_flag: # task got cancelled
        sim_log(DEBUG, self.env, self, "cancelled clock on ", self.p_inserv)
        self.cancel_flag = False
      # elif self.preempt_flag: # task got preempted
      #   self.p_l.insert(1, self.p_inserv)
      #   sim_log(DEBUG, self.env, self, "preempted ", self.p_inserv)
      #   self.preempt_flag = False
      else:
        sim_log(DEBUG, self.env, self, "serv done in {}s on ".format(self.env.now-clk_start_time), self.p_inserv)
        # self.qt_l.append(self.env.now - self.p_inserv.ref_time)
        if self.out is not None and self.p_inserv.sym != SYS_TRAFF_SYM:
          sim_log(DEBUG, self.env, self, "forwarding", self.p_inserv)
          self.p_inserv.prev_hop_id = self._id
          self.out.put(self.p_inserv)
        else:
          sim_log(DEBUG, self.env, self, "finished serv", self.p_inserv)
      self.p_inserv = None
  
  def put(self, p, preempt=False):
    p.ref_time = self.env.now
    self.p_l.append(p.deep_copy() )
    sim_log(DEBUG, self.env, self, "recved", p)
    
    # if preempt and (self.p_inserv is not None):
    #   self.preempt_flag = True
    #   self.preempt.succeed()
    #   self.p_l.insert(0, p)
    #   self.syncer.put(1)
    #   return
    return self.store.put(p)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    
    if cp._id == -1 and self.p_inserv is not None and self.p_inserv.sym == cp.sym: # for fairnessfirst
      self.p_inserv = None
      self.cancel_flag = True
    elif self.p_inserv is not None and self.p_inserv.job_id == cp._id:
      self.p_inserv = None
      self.cancel_flag = True
      
    # This has to be run always for mixed-traffic since requests can depart out of order
    for p in self.p_l:
      if p.job_id == cp._id:
        self.p_l.remove(p)
    
    if self.cancel_flag:
      self.cancel.succeed()

class QMonitor(object):
  def __init__(self, env, q, poll_interval):
    self.q = q
    self.env = env
    self.poll_interval = poll_interval
    
    self.pollt_l = []
    self.qlength_l = []
    self.action = env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_interval)
      
      self.pollt_l.append(self.env.now)
      self.qlength_l.append(self.q.length() )
