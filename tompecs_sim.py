import simpy, random

from log_utils import *

# *******************************  JobGen  ****************************** #
class Job(object):
  def __init__(self, _id, k, u):
    self._id = _id
    self.k = k
    self.u = u
  
  def __repr__(self):
    return "Job[id= {}, k= {}, u= {}]".format(self._id, self.k, self.u)
  
  def deep_copy(self):
    return Job(self._id, self.k, self.u)

class JobGen(object):
  def __init__(self, env, X, K, u=1):
    self.env = env
    self.X = X # inter-generation time
    self.K = K # number of tasks within a job
    self.u = u
    
    self.ngened = 0
    self.out = None
  
  def init(self):
    self.env.process(self.run() )
  
  def run(self):
    while 1:
      yield self.env.timeout(self.X.sample() )
      self.ngened += 1
      
      j = Job(_id=self.ngened, k=self.K.sample(), u=self.u)
      self.out.put(j)

class FCFS(object):
  def __init__(self, _id, env, V, out_m):
    self._id = _id
    self.env = env
    self.V = V # service time
    self.out_m = out_m
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
  
  def __repr__(self):
    return "FCFS[id= {}, V= {}]".format(self._id, self.V)
  
  def run(self):
    while True:
      j = yield self.store.get()
      
      slog(DEBUG, self.env, self, "starting clock on ", j)
      t = self.V.sample()*j.u
      yield self.env.timeout(t)
      slog(DEBUG, self.env, self, "serv done in {}s on ".format(t), j)
      
      self.out_m.put_m({
        'type': 'jdone',
        'jid': j._id} )
  
  def put(self, j):
    slog(DEBUG, self.env, self, "recved", j)
    return self.store.put(j)

class MixedFJQ(object):
  def __init__(self, _id, env, n, V):
    self._id = _id
    self.env = env
    self.n = n
  
    self.q_l = []
    for i in range(n):
      q = FCFS('FCFS-{}'.format(i), env, V, self)
      q.out = self
      self.q_l.append(q)
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    
    self.jid_info_m = {}
  
  def __repr__(self):
    return "MixedFJQ[id= {}, n= {}]".format(self._id, self.n)
  
  def run(self):
    while True:
      j = yield self.store.get()
      for q in random.sample(self.q_l, j.k):
        q.put(j.deep_copy() )
  
  def put(self, j):
    slog(DEBUG, self.env, self, "recved", j)
    self.jid_info_m[j._id] = {
      'entry_t': self.env.now,
      'k': j.k, 'finish_t_l': [] }
    return self.store.put(j)
  
  def put_m(self, m):
    slog(DEBUG, self.env, self, "recved", m)
    
    if m['type'] == 'jdone':
      jid = m['jid']
      self.jid_info_m[jid]['finish_t_l'].append(self.env.now)
    else:
      log(ERROR, "Unrecognized; m['type']= {}".format(m['type'] ) )
  
      