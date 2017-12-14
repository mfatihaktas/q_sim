from sim import *
from patch import *
from rvs import *

# *************************************  Mixed Packet Generator  ********************************* #
class MixedPG(object):
  def __init__(self, env, _id, qarrdist_m_l):
    self._id = _id
    self.env = env
    self.qarrdist_m_l = qarrdist_m_l
    
    self.i__n_sent = len(qarrdist_m_l)*[0]
    for i, dist_m in enumerate(qarrdist_m_l):
      self.env.process(self.run(i, dist_m) )
    self.out = None
  
  def __repr__(self):
    return "MixedPG[_id={}, qarrdist_m_l=\n {}]".format(self._id, self.qarrdist_m_l)
  
  def run(self, i, dist_m):
    while 1:
      dist = dist_m['dist']
      if dist == 'Exp':
        rv = Exp(dist_m['mu'] )
      elif dist == 'Pareto':
        rv = Pareto(dist_m['loc'], dist_m['a'])
      yield self.env.timeout(rv.gen_sample() )
      self.i__n_sent[i] += 1
      self.out.put(Packet(time=self.env.now, _id=self.i__n_sent[i], flow_id=i) )

# ********************************************  Slave Q  ***************************************** #
class SlaveQ(Q): # Release HoL at command
  def __init__(self, _id, env):
    super().__init__(_id, env)
    
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
    return sum(self.qt_l)/len(self.qt_l)
    # nonzero_qt_l = [t for t in self.qt_l if t > 0.000001]
    # return sum(nonzero_qt_l)/len(nonzero_qt_l)
  
  def avg_qtime2(self):
    return sum([t**2 for t in self.qt_l] )/len(self.qt_l)
    
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
  def __init__(self, env, n, k, attacker=None):
    self.env = env
    self.n = n
    self.k = k
    self.attacker = attacker
    
    self.i_q_map = []
    for i in range(self.n):
      self.i_q_map.append(SlaveQ(_id=i, env=env) )
    
    self.start_time = env.now
    # self.store = simpy.Store(env)
    # env.process(self.run() )
    
    self.attack_on = False
    self.attacker_startt = None
    self.attack_success = None
    self.attackt_l = []
    env.process(self.run_attacker() )
  
  def __repr__(self):
    return "MixedNet[n={}, k={}]".format(self.n, self.k)
  
  def state(self):
    ql_l = [len(self.i_q_map[i] ) for i in range(self.n) ]
    return ','.join(map(str, ql_l) )
  
  def ET_ET2(self):
    # Assuming each q is identical
    ET_sum, ET2_sum = 0, 0
    for q in self.i_q_map:
      ET_sum += q.avg_qtime()
      ET2_sum += q.avg_qtime2()
    return ET_sum/self.n, ET2_sum/self.n
  
  def throughput(self):
    n_released = sum([q.n_released for i,q in enumerate(self.i_q_map) ] )
    return n_released/(self.env.now - self.start_time)
  
  # def run(self):
  #   while True:
  #     p = (yield self.store.get() )
  #     self.i_q_map[p.flow_id].put(p)
      
  #     n_busy = 0
  #     for i,q in enumerate(self.i_q_map):
  #       if q.length():
  #         n_busy += 1
  #     if n_busy >= self.k:
  #       for i,q in enumerate(self.i_q_map):
  #         q.release()
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    # return self.store.put(p)
    
    self.i_q_map[p.flow_id].put(p)
    
    if self.attack_on:
      self.attacker.in_packet(p.flow_id)
    
    busy_qid_l = []
    for i, q in enumerate(self.i_q_map):
      if q.length():
        busy_qid_l.append(q.length() )
    
    if len(busy_qid_l) >= self.k:
      for i, q in enumerate(self.i_q_map):
        q.release()
      
      if self.attack_on:
        self.attacker.out_frame(busy_qid_l)
        
        if self.attacker.state() == self.state():
          self.attackt_l.append(self.env.now - self.attacker_startt)
          self.attack_on = False
          self.attack_success.succeed()
  
  def run_attacker(self):
    while True:
      yield self.env.timeout(5000)
      self.attack_on = True
      self.attacker_startt = self.env.now
      self.attacker.reset()
      
      self.attack_success = self.env.event()
      yield (self.attack_success)

class Attacker(object):
  def __init__(self, env, n, k):
    self.env = env
    self.n = n
    self.k = k
    
    self.i_ql_map = [0]*self.n
  
  def state(self):
    return ','.join(map(str, self.i_ql_map) )
  
  def reset(self):
    self.i_ql_map.clear()
  
  def in_packet(self, i):
    self.i_ql_map[i] += 1
  
  def out_frame(self, o_l):
    for o in o_l:
      _l = self.i_ql_map[o]
      self.i_ql_map[o] = max(_l-1, 0)

'''
class Deanonymizer(object):
  def __init__(self, env, n, k):
    self.env = env
    self.n = n
    self.k = k
    
    self.deanont_l = []
    
    self.attack_on = False
    self.empty = None
    self.success = None
    env.process(self.run() )
    
    self.deanon_startt = None
    self.i__o_l_m = {i:None for i in range(self.n) }
    self.in_frame = []
  
  def run(self):
    while True:
      yield self.env.timeout(5000)
      
      self.attack_on = True
      self.deanon_startt = self.env.now
      for i in self.i__o_l_m:
        self.i__o_l_m[i] = [j for j in range(self.n) ]
      
      self.empty = self.env.event()
      yield (self.empty)
      
      self.success = self.env.event()
      yield (self.success)
  
  def in_packet(self, i):
    if self.attack_on:
      self.in_frame.append(i)
  
  def out_frame(self, o_l):
    if not self.attack_on:
      return
    
    if len(self.in_frame) < self.k: # not empty yet
      self.in_frame.clear()
      return
    elif self.empty is not None and len(self.in_frame) == self.k: # now empty
      self.empty.succeed()
      self.empty = None
    # else:
      for i in self.in_frame:
        for o in self.i__o_l_m[i]:
          if o not in o_l:
            self.i__o_l_m[i].remove(o)
      self.in_frame.clear()
      
      for i, o_l in self.i__o_l_m.items():
        if len(o_l) > 1:
          return
      self.deanont_l.append(self.env.now - self.deanon_startt)
      self.attack_on = False
      self.success.succeed()
'''

class MNMonitor(object):
  def __init__(self, env, mn, poll_interval):
    self.env = env
    self.mn = mn
    self.poll_interval = poll_interval
    
    self.qid__scounter_map_map = {}
    for i in range(self.mn.n):
      self.qid__scounter_map_map[i] = {}
    env.process(self.run() )
  
  def __repr__(self):
    return "Monitor:{}".format(self.mn)
  
  def ss_prob_map(self):
    # Assuming each q is identical
    scounter_map = {}
    for i in range(self.mn.n):
      for s, c in self.qid__scounter_map_map[i].items():
        if s not in scounter_map:
          scounter_map[s] = 0
        scounter_map[s] += c
    
    sum_c = sum([c for s,c in scounter_map.items() ] )
    return {s:c/sum_c for s,c in scounter_map.items() }
  
  def EL_EL2(self):
    # Assuming each q is identical
    scounter_map = {}
    for i in range(self.mn.n):
      for s, c in self.qid__scounter_map_map[i].items():
        if s not in scounter_map:
          scounter_map[s] = 0
        scounter_map[s] += c
    
    sum_c, sum_s, sum_s2 = 0, 0, 0
    for s, c in scounter_map.items():
      sum_c += c
      sum_s += c*s
      sum_s2 += c*s**2
    return sum_s/sum_c, sum_s2/sum_c
  
  def run(self):
    while True:
      yield self.env.timeout(self.poll_interval)
      
      for i in range(self.mn.n):
        scounter_map = self.qid__scounter_map_map[i]
        s = self.mn.id_q_map[i].length()
        # print("polled s= {}".format(s) )
        if s not in scounter_map:
          scounter_map[s] = 0
        scounter_map[s] += 1
  