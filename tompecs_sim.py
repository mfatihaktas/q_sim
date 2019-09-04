import simpy, random

from sim import FCFS
from mds_sim import MDSQ

"""
  Implements a FJ Q over (n, k, r, t)-LRC.
  k objects are encoded into n codewords.
  The (r, t)-availability ensures that each systematic node can be regenerated using one of the t
  disjoint repair groups of other servers, each of size at most r (typically << k)
  [When do the Availability Codes Make the Stored Data More Available?]
"""
class AVQ(object): # Availability
  def __init__(self, _id, env, t, r, k, servdist_m, sching, out=None):
    self._id = _id
    self.env = env
    self.t = t
    self.r = r
    self.k = k
    self.sching = sching
    self.out = out
    
    self.jsink = JSink(_id, env)
    self.jsink.out = out
    self.id__q_map = {}
    
    num_q = int(1 + t*r)
    qid_l = [i for i in range(num_q) ]
    self.qid_l = []
    for g in range(1, t + 1 + 1):
      if g == 1:
        self.qid_l.append(qid_l[0] )
      else:
        li = 1+(g-2)*r
        ri = li + r
        self.qid_l.append("FJ{}".format(qid_l[li:ri] ) )
    
    self.join_q = JQ(_id, env, k=1, input_qid_l=self.qid_l)
    self.join_q.out = self.jsink
    self.join_q.out_c = self
    self.join_q.out_m = None # can be set by the caller if desired
    for g in range(1, t + 1 + 1):
      q = None
      if g == 1:
        q = FCFS(self.qid_l[g-1], env, servdist_m)
        log(DEBUG, "g= {}, q= {}".format(g, q) )
        q.out = self.join_q
      else:
        li = 1+(g-2)*r
        ri = li + r
        q = MDSQ("FJ{}".format(list(qid_l[li:ri] ) ), env, k, qid_l[li:ri], servdist_m, out=self.join_q)
        log(DEBUG, "g= {}, q= {}".format(g, q) )
      self.id__q_map[g] = q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.action = env.process(self.run() )
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
    self.servtype__num_m = (t+1)*[0]
  
  def __repr__(self):
    return "AVQ[k= {}, r= {}, t= {}]".format(self.k, self.r, self.t)
    # return "AVQ[k= {}, r= {}, t= {}, qid_l= [{}] ]".format(self.k, self.r, ",".join(self.qid_l) )
  
  def state(self, job_id_to_exclude=[]):
    state = []
    for g, q in self.id__q_map.items():
      state += q.state(job_id_to_exclude)
    
    return state
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      if self.sching == "rep-to-all":
        for i, q in self.id__q_map.items():
          q.put(p.deep_copy() )
      else: # select-one at random
        id_l = [i for i,q in self.id__q_map.items() ]
        i = id_l[random.randint(0, len(id_l)-1) ]
        self.id__q_map[i].put(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.job_id_counter += 1
    p.job_id = self.job_id_counter
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      #
      next_job_id = cp._id + 1
      type_ = 0
      for g, q in self.id__q_map.items():
        if g == 1:
          if not q._in(next_job_id): break
        else:
          if q.n_servers_in(next_job_id) < self.r: type_ += 1
      self.servtype__num_m[type_] += 1
      #
      for g, q in self.id__q_map.items():
        if q._id not in cp.departed_qid_l:
          q.put_c(cp.deep_copy() )
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)


