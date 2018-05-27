import numpy as np

from sim import *
from patch import *
from rvs import *

# ********************************************  Slave Q  ***************************************** #
class SlaveQ(Q): # Release HoL at command
  def __init__(self, _id, env):
    super().__init__(_id, env)
    self.p_l = []
    
    self.n_recved = 0
    self.n_released = 0
    self.qt_l = []
  
  def __repr__(self):
    return "SlaveQ[id={}]".format(self._id)
  
  def length(self):
    return len(self.p_l)
  
  def avg_qtime(self):
    return np.mean(self.qt_l)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    self.n_recved += 1
    p.ref_time = self.env.now
    self.p_l.append(p)
  
  def release(self):
    if len(self.p_l):
      p = self.p_l.pop(0)
      self.qt_l.append(self.env.now - p.ref_time)
      sim_log(DEBUG, self.env, self, "released", p)
      self.n_released += 1

# *************************************  SamplekMix  ****************************************** #
# n servers with Poisson arrivals, every time a message arrives, k queues are selected at uniformly random and released.
class SamplekMix(object):
  def __init__(self, env, n, k, pd=None):
    self.env = env
    self.n = n
    self.k = k
    self.pd = pd if pd is not None else 1/n
    
    self.i_q_m = []
    for i in range(self.n):
      self.i_q_m.append(SlaveQ(_id=i, env=env) )
    
  def __repr__(self):
    return "SamplekMix[n={}, k={}]".format(self.n, self.k)
  
  def qt_l(self):
    l = []
    for q in self.i_q_m:
      l.extend(q.qt_l)
    return l
  
  def ET(self):
    avg_qtime_l = [q.avg_qtime() for q in self.i_q_m]
    log(WARNING, "avg_qtime_l= {}".format(avg_qtime_l) )
    return np.mean([q.avg_qtime() for q in self.i_q_m] )
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    self.i_q_m[p.flow_id].put(p)
    
    l = list(range(self.n) )
    l.remove(p.flow_id)
    if np.random.uniform(0, 1) <= self.pd:
      i_l = [p.flow_id] + [l[i] for i in np.random.choice(self.n-1, self.k-1, replace=False) ]
    else:
      i_l = [l[i] for i in np.random.choice(self.n-1, self.k, replace=False) ]
    # print("i_l= {}".format(i_l) )
    for i in i_l:
      self.i_q_m[i].release()
    
  