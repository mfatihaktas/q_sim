import simpy, random, copy, pprint
from simpy.core import BoundClass
from simpy.resources import base
from heapq import heappush, heappop

from simplex_sim_components import *
from patch import *
from deprecated import *

class MDS_PacketGenerator(PacketGenerator):
  def __init__(self, env, _id, mds_arr_dist, size_dist, initial_delay=0, finish=float("inf"), flow_id=0, sys_arr_dist=None):
    super().__init__(env, _id, mds_arr_dist, size_dist, initial_delay, finish, flow_id)
    self.sys_arr_dist = sys_arr_dist
    
    self.mds_n_sent = 0
    self.sys_n_sent = 0
    self.out = None
    
    self.mds_action = env.process(self.mds_run() )  # starts the run() method as a SimPy process
    if sys_arr_dist is not None:
      self.sys_action = env.process(self.sys_run() )
  
  def __repr__(self):
    return "MDS_PacketGenerator[mds_n_sent= {}, sys_n_sent= {}]".format(self.mds_n_sent, self.sys_n_sent)
  
  def mds_run(self):
    yield self.env.timeout(self.initial_delay)
    while self.env.now < self.finish:
      yield self.env.timeout(self.arr_dist() )
      self.mds_n_sent += 1
      self.out.put(Packet(time=self.env.now, size=self.size_dist(), _id=self.mds_n_sent, sym=MDS_TRAFF_SYM) )
  
  def sys_run(self):
    yield self.env.timeout(self.initial_delay)
    while self.env.now < self.finish:
      yield self.env.timeout(self.sys_arr_dist() )
      self.sys_n_sent += 1
      self.out.put(Packet(time=self.env.now, size=self.size_dist(), _id=self.sys_n_sent, sym=REGULAR_TRAFF_SYM) )
  
# *******************************************  MDSQ  ********************************************* #
# Fork incoming job to all sub-q's, wait for any k task to complete, cancel the remainings.
class MDSQ(object):
  def __init__(self, _id, env, k, qid_l, qserv_rate_l, r=None, preempt=False, out=None):
    self._id = _id
    self.env = env
    self.n = len(qid_l)
    self.r = r
    self.k = k
    self.qid_l = qid_l
    self.qserv_rate_l = qserv_rate_l
    self.preempt = preempt
    # self.out = out
    
    self.join_sink = JSink(_id, env)
    self.join_sink.out = out
    self.join_q = JQ(_id, env, k, qid_l)
    self.join_q.out = self.join_sink
    self.join_q.out_c = self
    self.id_q_map = {}
    for i in range(self.n):
      qid = qid_l[i]
      m1_q = S1_Q(_id=qid, env=env, serv_rate=qserv_rate_l[i] )
      m1_q.out = self.join_q
      self.id_q_map[qid] = m1_q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.action = env.process(self.run() ) # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
  
  def __repr__(self):
    return "MDSQ[k= {}, qid_l= [{}] ]".format(self.k, ",".join(self.qid_l) )
  
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
      self.qid_to_split_l = []
      if p.sym == MDS_TRAFF_SYM:
        if self.r is None:
          self.qid_to_split_l = self.qid_l
        else:
          while len(self.qid_to_split_l) < self.r:
            qid = self.qid_l[random.randint(0, self.n-1) ]
            if qid not in self.qid_to_split_l:
              self.qid_to_split_l.append(qid)
      elif p.sym == REGULAR_TRAFF_SYM:
        self.qid_to_split_l = self.qid_l
      else:
        log(ERROR, "Unexpected p.sym= {}".format(p.sym) )
        return 1
      for qid in self.qid_to_split_l:
        self.id_q_map[qid].put(p, preempt=self.preempt)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    if p.entrance_time is None:
      p.entrance_time = self.env.now
    if p.job_id is None:
      self.job_id_counter += 1
      p.job_id = self.job_id_counter
    return self.store.put(p)
  
  # def put(self, p):
  #   sim_log(DEBUG, self.env, self, "recved", p)
  #   if p.entrance_time is None:
  #     p.entrance_time = self.env.now
  #   if p.job_id is None:
  #     self.job_id_counter += 1
  #     p.job_id = self.job_id_counter
  #   # return self.store.put(p)
  #   self.qid_to_split_l = []
  #   if p.sym == MDS_TRAFF_SYM:
  #     if self.r is None:
  #       self.qid_to_split_l = self.qid_l
  #     else:
  #       while len(self.qid_to_split_l) < self.r:
  #         qid = self.qid_l[random.randint(0, self.n-1) ]
  #         if qid not in self.qid_to_split_l:
  #           self.qid_to_split_l.append(qid)
  #   elif p.sym == REGULAR_TRAFF_SYM:
  #     # self.qid_to_split_l = self.qid_l
  #     self.qid_to_split_l.append(self.qid_l[random.randint(0, self.n-1) ] )
  #   else:
  #     log(ERROR, "Unexpected p.sym= {}".format(p.sym) )
  #     return 1
  #   for qid in self.qid_to_split_l:
  #     self.id_q_map[qid].put(p, preempt=self.preempt)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      for i, q in self.id_q_map.items():
        if q._id not in cp.departed_qid_l:
          q.put_c(cp)
      if cp.prev_hop_id != self._id: # To avoid inifinite loops forward only when cp comes from a different JQ than self
        self.join_q.put_c(cp)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)
  
  # def put_c(self, cp):
  #   sim_log(DEBUG, self.env, self, "recved", cp)
  #   # return self.store_c.put(cp)
  #   for i, q in self.id_q_map.items():
  #     if q._id not in cp.departed_qid_l:
  #       q.put_c(cp)
  #   if cp.prev_hop_id != self._id: # To avoid inifinite loops forward only when cp comes from a different JQ than self
  #     self.join_q.put_c(cp)
  
class MDSQMonitor(object):
  def __init__(self, env, q, poll_dist):
    self.q = q
    self.env = env
    self.poll_dist = poll_dist
    
    self.t_l = [] # Time steps that the numbers polled from the aq
    # self.n_l = [] # Num of packets in the aq
    self.state__counter_map = {}
    
    self.action = env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_dist() )
      
      self.t_l.append(self.env.now)
      state = self.q.state()
      # print("state= {}".format(list_to_str(state) ) )
      num_job = max(state)
      # self.n_l.append(num_job)
      # rel_state = ",".join("%d" % (num_job-l) for l in state)
      rel_state = list_to_str(state)
      if rel_state not in self.state__counter_map:
        self.state__counter_map[rel_state] = 0
      self.state__counter_map[rel_state] += 1
