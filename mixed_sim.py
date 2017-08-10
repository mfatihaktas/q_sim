from sim import *
from patch import *

# *************************************  Mixed Packet Generator  ********************************* #
class MixedPG(object):
  def __init__(self, env, _id, qlambda_l):
    self._id = _id
    self.env = env
    self.qlambda_l = qlambda_l
    
    self.i__n_sent = len(qlambda_l)*[0]
    for i,l in enumerate(qlambda_l):
      self.env.process(self.run(i, l) )
    self.out = None
  
  def __repr__(self):
    return "MixedPG[_id={}, qlambda_l={}]".format(self._id, self.qlambda_l)
  
  def run(self, i, l):
    while 1:
      yield self.env.timeout(random.expovariate(l) )
      self.i__n_sent[i] += 1
      self.out.put(Packet(time=self.env.now, _id=self.i__n_sent[i], flow_id=i) )

# ********************************************  Slave Q  ***************************************** #
class SlaveQ(Q): # Release HoL at command
  def __init__(self, _id, env):
    super().__init__(_id, env, 1)
    
    self.p_l = []
    self.n_recved = 0
    self.n_released = 0
    self.qt_l = []
    
    # self.store = simpy.Store(env)
    # env.process(self.run() )

  def __repr__(self):
    return "SlaveQ[_id={}]".format(self._id)
  
  def length(self):
    return len(self.p_l)
  
  def avg_qtime(self):
    # return sum(self.qt_l)/len(self.qt_l)
    nonzero_qt_l = [t for t in self.qt_l if t > 0.000001]
    return sum(nonzero_qt_l)/len(nonzero_qt_l)
  
  # def run(self):
  #   while True:
  #     p = (yield self.store.get() )
  #     self.p_l.append(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    self.n_recved += 1
    p.ref_time = self.env.now
    # return self.store.put(p)
    
    self.p_l.append(p)
  
  def release(self):
    if len(self.p_l):
      p = self.p_l.pop(0)
      self.qt_l.append(self.env.now - p.ref_time)
      sim_log(DEBUG, self.env, self, "released", p)
      self.n_released += 1

# *************************************  Mixed Net  ****************************************** #
# n servers with Poisson arrivals, once any k servers are busy, Hol is immediately released
class MixedNet(object): # Network
  def __init__(self, env, n, k):
    self.env = env
    self.n = n
    self.k = k
    
    self.id_q_map = []
    for i in range(self.n):
      self.id_q_map.append(SlaveQ(_id=i, env=env) )
    
    self.start_time = env.now
    # self.store = simpy.Store(env)
    # env.process(self.run() )
  
  def __repr__(self):
    return "MixedNet[n={}, k={}]".format(self.n, self.k)
  
  def E_T(self):
    # Assuming each q is identical
    E_T_sum = 0
    for q in self.id_q_map:
      E_T_sum += q.avg_qtime()
    return E_T_sum/self.n
  
  def throughput(self):
    n_released = sum([q.n_released for i,q in enumerate(self.id_q_map) ] )
    return n_released/(self.env.now - self.start_time)
  
  # def run(self):
  #   while True:
  #     p = (yield self.store.get() )
  #     self.id_q_map[p.flow_id].put(p)
      
  #     n_busy = 0
  #     for i,q in enumerate(self.id_q_map):
  #       if q.length():
  #         n_busy += 1
  #     if n_busy >= self.k:
  #       for i,q in enumerate(self.id_q_map):
  #         q.release()
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    # return self.store.put(p)
    
    self.id_q_map[p.flow_id].put(p)
    
    n_busy = 0
    for i,q in enumerate(self.id_q_map):
      if q.length():
        n_busy += 1
    if n_busy >= self.k:
      for i,q in enumerate(self.id_q_map):
        q.release()

class MixedNetMonitor(object):
  def __init__(self, env, mixednet, poll_interval):
    self.env = env
    self.mixednet = mixednet
    self.poll_interval = poll_interval
    
    self.qid__state_counter_map_map = {}
    for i in range(self.mixednet.n):
      self.qid__state_counter_map_map[i] = {}
    env.process(self.run() )
  
  def __repr__(self):
    return "Monitor:{}".format(self.mixednet)
  
  def steadystate_prob_map(self):
    # Assuming each q is identical
    state_counter_map = {}
    for i in range(self.mixednet.n):
      for s, c in self.qid__state_counter_map_map[i].items():
        if s not in state_counter_map:
          state_counter_map[s] = 0
        state_counter_map[s] += c
    
    total_c = sum([c for s,c in state_counter_map.items() ] )
    return {s:c/total_c for s,c in state_counter_map.items() }
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_interval)
      
      for i in range(self.mixednet.n):
        state_counter_map = self.qid__state_counter_map_map[i]
        s = self.mixednet.id_q_map[i].length()
        # print("polled s= {}".format(s) )
        if s not in state_counter_map:
          state_counter_map[s] = 0
        state_counter_map[s] += 1
  