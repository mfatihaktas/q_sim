import simpy, random, copy, pprint

from sim_components import *
from patch import *

class MDS_PG(PG):
  def __init__(self, env, _id, inter_arr_dist):
    super().__init__(env, _id, inter_arr_dist, flow_id=0)
    
    self.n_sent = 0
    self.out = None
    
    env.process(self.run() )
  
  def __repr__(self):
    return "MDS_PG[n_sent= {}]".format(self.n_sent)
  
  def run(self):
    while 1:
      yield self.env.timeout(self.inter_arr_dist() )
      self.n_sent += 1
      self.out.put(Packet(time=self.env.now, size=1, _id=self.n_sent, sym=MDS_TRAFF_SYM) )
  
# *******************************************  MDSQ  ********************************************* #
# Split arrivals to all, wait for any k for download, cancel the outstanding remainings.
class MDSQ(object):
  # r: split each arrival to randomly any r servers
  def __init__(self, _id, env, k, qid_l, qmu_l, serv="Exp", r=None, preempt=False, out=None):
    self._id = _id
    self.env = env
    self.n = len(qid_l)
    self.r = r
    self.k = k
    self.qid_l = qid_l
    self.qmu_l = qmu_l
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
      q = FCFS(_id=qid, env=env, rate=qmu_l[i], serv=serv)
      q.out = self.join_q
      self.id_q_map[qid] = q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.action = env.process(self.run() ) # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
    self.servtype__num_m = k*[0]
  
  def __repr__(self):
    return "MDSQ[k= {}, qid_l= {} ]".format(self.k, self.qid_l)
  
  def length(self):
    return max([q.length() for i, q in self.id_q_map.items() ] )
  
  def state(self, job_id_to_exclude=[]):
    # return [q.length() for i, q in self.id_q_map.items() ]
    state = []
    for i, q in self.id_q_map.items():
      state += q.state(job_id_to_exclude)
    return state
  
  def n_servers_in(self, _id):
    n = 0
    for i,q in self.id_q_map.items():
      if q._in(_id): n += 1
    return n
  
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
      else:
        self.qid_to_split_l = self.qid_l
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
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      # 
      next_job_id = cp._id + 1
      type_ = self.n - self.n_servers_in(next_job_id)
      if type_ == self.n: # next job is not in
        type_ = 0
      self.servtype__num_m[type_] += 1
      # 
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
