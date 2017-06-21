import numpy, pprint

# from sim import *
from simplex_sim import *
from fairness_model import *
from patch import *

POP_SYM = 'a'

class FF_JSink(object): # Join
  def __init__(self, _id, env):
    self._id = _id
    self.env = env
    
    self.st_pop_l = []
    self.st_unpop_l = [] # system time
    self.qid__num_win_map = {}
    
    self.store = simpy.Store(env)
    self.out = None
  
  def __repr__(self):
    return "JSink[_id= {}]".format(self._id)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    
    st = self.env.now - p.entrance_time
    if p.sym == POP_SYM:
      self.st_pop_l.append(st)
      
      if p.winner_id not in self.qid__num_win_map:
        self.qid__num_win_map[p.winner_id] = 0
      self.qid__num_win_map[p.winner_id] += 1
    else:
      self.st_unpop_l.append(st)
    
    if self.out is not None: 
      self.out.put(p)

# **********************************  Fairness First AVQ  ******************************* #
class FF_AVQ(object): # Fairness First
  def __init__(self, _id, env, t, qmu_l, serv, sym__rgroup_l_map, sym__sysqid_map, out=None):
    self._id = _id
    self.env = env
    self.sym__rgroup_l_map = sym__rgroup_l_map
    self.sym__sysqid_map = sym__sysqid_map
    self.out = out
    
    self.num_q = len(qmu_l)
    self.qid_l = [i for i in range(self.num_q) ]
    
    self.jsink = FF_JSink(_id, env)
    self.jsink.out = out
    self.id_q_map = {}
    
    self.jq = MT_JQ(_id=_id, env=env, input_qid_l=self.qid_l, sym__rgroup_l_map=sym__rgroup_l_map)
    self.jq.out = self.jsink # data outlet
    self.jq.out_c = self # control outlet
    for i, qid in enumerate(self.qid_l):
      q = FCFS(_id=qid, env=env, rate=qmu_l[i], serv=serv)
      q.out = self.jq
      self.id_q_map[qid] = q
    # 
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    env.process(self.run() )
    env.process(self.run_c() )
    
    self.job_id_counter = 0
    self.pop_l = []
    self.starttype__num_map = {i:0 for i in range(t+1) }
    
    self.pop_store = simpy.Store(env)
    env.process(self.send_pop() )
    self.release_pop = None
  
  def __repr__(self):
    return "MT_AVQ[qid_l= {}]".format(self.qid_l)
  
  def state(self):
    return {i:q.length() for i, q in self.id_q_map.items() }
  
  def send_pop(self):
    while True:
      p = (yield self.pop_store.get() )
      
      sys_q = self.id_q_map[self.sym__sysqid_map[p.sym] ]
      if sys_q.busy:
        self.release_pop = self.env.event()
        yield self.release_pop
      
      sys_q.put(p.deep_copy() )
      st = 0
      for r_l in self.sym__rgroup_l_map[p.sym]:
        if len(r_l) == 1: continue
        
        rep = True
        for r in r_l:
          q = self.id_q_map[r]
          if q.busy and q.p_in_serv.sym != POP_SYM: # q may be busy with pop_sym because simpy did not register cancellation yet
            print("server-{} busy with sym= {}, while sending p= {}".format(r, q.p_in_serv.sym, p) )
            rep = False
            break
        if rep:
          for r in r_l:
            self.id_q_map[r].put(p.deep_copy() )
          st += 1
      self.starttype__num_map[st] += 1
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      
      if p.sym == POP_SYM:
        self.pop_store.put(p)
      else:
        sys_q = self.id_q_map[self.sym__sysqid_map[p.sym] ]
        if sys_q.p_in_serv is not None and sys_q.p_in_serv.sym == POP_SYM:
          sys_qid = self.sym__sysqid_map[p.sym]
          for r_l in self.sym__rgroup_l_map[POP_SYM]:
            if sys_qid in r_l:
              for r in r_l:
                self.id_q_map[r].put_c(CPacket(_id=sys_q.p_in_serv.job_id, prev_hop_id=self._id) )
                # print("sending cancel to q= {}".format(self.id_q_map[r] ) )
              
        sys_q.put(p.deep_copy() )
  
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
      # log(WARNING, "recved cp= {}".format(cp) )
      # print("state= {}".format(pprint.pformat(self.state() ) ) )
      
      for g, q in self.id_q_map.items():
        if q._id not in cp.departed_qid_l:
          q.put_c(cp.deep_copy() )
      
      if cp.sym == POP_SYM:
        if self.release_pop is not None:
          self.release_pop.succeed()
          self.release_pop = None
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

# *********************************  Exp  *********************************** #
def test_ff_avq(num_f_run, arr_rate, mu, t, r, k, serv, pop_arr_rate):
  num_q = int(1 + t*r)
  qmu_l = num_q*[mu]
  
  E_T_pop_f_sum, E_T_unpop_f_sum = 0, 0
  for f in range(num_f_run):
    log(WARNING, "arr_rate= {}, pop_arr_rate= {}, mu= {}, t= {}, r= {}, k= {}, qmu_l= {}".format(arr_rate, pop_arr_rate, mu, t, r, k, qmu_l) )
    env = simpy.Environment()
    
    sym__sysqid_map, sym__rgroup_l_map = {}, {}
    if t == 1: sym_l = ['a', 'b']
    elif t == 3: sym_l = ['a', 'b', 'c']
    for sym in sym_l:
      rgroup_l = []
      if t == 1:
        if sym == 'a':
          sym__sysqid_map[sym] = 0
          rgroup_l.append([0] )
          rgroup_l.append([1, 2] )
        elif sym == 'b':
          sym__sysqid_map[sym] = 1
          rgroup_l.append([1] )
          rgroup_l.append([0, 2] )
      elif t == 3:
        if sym == 'a':
          sym__sysqid_map[sym] = 0
          rgroup_l.append([0] )
          rgroup_l.append([1, 2] )
          rgroup_l.append([3, 4] )
          rgroup_l.append([5, 6] )
        elif sym == 'b':
          sym__sysqid_map[sym] = 1
          rgroup_l.append([1] )
          rgroup_l.append([0, 2] )
          rgroup_l.append([3, 5] )
          rgroup_l.append([4, 6] )
        elif sym == 'c':
          sym__sysqid_map[sym] = 3
          rgroup_l.append([3] )
          rgroup_l.append([0, 4] )
          rgroup_l.append([1, 5] )
          rgroup_l.append([2, 6] )
      sym__rgroup_l_map[sym] = rgroup_l
    pg = MT_PG(env, "pg", arr_rate, sym_l, pop_sym=POP_SYM, pop_arr_rate=pop_arr_rate)
    log(WARNING, "sym__rgroup_l_map=\n {}\n sym__sysqid_map=\n {}".format(pprint.pformat(sym__rgroup_l_map), pprint.pformat(sym__sysqid_map) ) )
    avq = FF_AVQ("ff_avq", env, t, qmu_l, serv, sym__rgroup_l_map, sym__sysqid_map)
    pg.out = avq
    pg.init()
    env.run(until=50000)
    
    l = avq.jsink.st_pop_l
    # print("avq.jsink.st_pop_l= {}".format(l) )
    if len(l): E_T_pop_f_sum += float(sum(l) )/len(l)
    l = avq.jsink.st_unpop_l
    if len(l): E_T_unpop_f_sum += float(sum(l) )/len(l)
    
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    print("pg.sym__n_sent= {}".format(pg.sym__n_sent) )
    
    print("avq.starttype__num_map= {}".format(pprint.pformat(avq.starttype__num_map) ) )
    
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  E_T_pop = E_T_pop_f_sum/num_f_run
  E_T_unpop = E_T_unpop_f_sum/num_f_run
  print(">> E_T_pop= {}, E_T_unpop= {}".format(E_T_pop, E_T_unpop) )
  if E_T_unpop > 100: return (None, None)
  return (E_T_pop, E_T_unpop)
  
def plot_ff_simplex():
  t, r, k = 1, 2, 2
  serv = "Exp" # "Dolly" # "Exp"
  mu = 1
  unpop_arr_rate = 0.0001 # 0.05 # 0.5
  if serv == "Exp":
    if t == 1: arr_rate_ub = 0.95
    elif t == 3: arr_rate_ub = 2.4
    elif t == 7: arr_rate_ub = 3
  elif serv == "Dolly":
    if t == 1: arr_rate_ub = 0.28
    elif t == 3: arr_rate_ub = 0.4
  log(WARNING, "t= {}, r= {}, k= {}, mu= {}, serv= {}, arr_rate_ub= {}".format(t, r, k, mu, serv, arr_rate_ub) )
  
  E_T_unpop_sim_simplex_l, E_T_pop_sim_simplex_l, E_T_simplex_l = [], [], []
  E_T_pop_ub_simplex_l, E_T_pop_lb_simplex_l = [], []
  
  sim_simplex = False
  if serv == "Exp":
    if t == 11:
      if unpop_arr_rate == 0.0001:
        E_T_pop_sim_simplex_l= [
          0.987915280539363,
          1.304227687592849,
          1.5861329900584953,
          2.0338480009189297,
          2.6568538834829756,
          4.25223189197621,
          4.572270849279754,
          4.772813399725301,
          5.745868507381323,
          6.806214731538167,
          7.509779490220204,
          8.512967502760787,
          9.711561399133826,
          11.721099073290514,
          23.72900999922603]
      elif unpop_arr_rate == 0.05:
        E_T_pop_sim_simplex_l= [
          0.8848064222998808,
          1.2293433967049092,
          1.4965561257946725,
          1.8706965341214936,
          2.6273876667369085,
          4.058122745636425,
          4.475293417500701,
          5.169425997589906,
          5.155439787006647,
          6.424170233975005,
          8.091666619742423,
          9.059079920569532,
          8.922207674393974,
          14.597795898918362,
          16.089683306249366]
      elif unpop_arr_rate == 0.5:
        E_T_pop_sim_simplex_l= [
          1.0845274986990392,
          1.18656625530774,
          1.483262548844013,
          1.861180534019326,
          2.674900768108865,
          4.052899648396583,
          4.6692017115889,
          5.33978762187163,
          5.883951181517289,
          6.259235022219047,
          7.214618010931689,
          9.075137367876843,
          10.538302731496527,
          12.159546579037903,
          21.76150746101844]
    elif t == 33:
      pass
    elif t == 77:
      pass
    else: sim_simplex = True
  if serv == "Dolly":
    if t == 11:
      pass
    elif t == 33:
      pass
    else: sim_simplex = True
  
  num_f_run = 2
  arr_rate_l = []
  # for arr_rate in [*numpy.linspace(0.05, 0.8*arr_rate_ub, 5, endpoint=False), *numpy.linspace(0.8*arr_rate_ub, arr_rate_ub, 10) ]:
  for arr_rate in numpy.linspace(0.05, arr_rate_ub, 2):
    arr_rate_l.append(arr_rate)
    if sim_simplex:
      (E_T_pop, E_T_unpop) = test_ff_avq(num_f_run, unpop_arr_rate, mu, t, r, k, serv, pop_arr_rate=arr_rate)
      E_T_pop_sim_simplex_l.append(E_T_pop)
      E_T_unpop_sim_simplex_l.append(E_T_unpop)
    
    # E_T_pop_ub_ff_simplex(arr_rate, t, mu)
    E_T_pop_ub_simplex_l.append(1/(mu - arr_rate) )
    E_T_pop_lb_simplex_l.append(E_T_pop_lb_ff_simplex(arr_rate, t, mu) )
  
  mew, ms = 3, 8
  log(WARNING, "E_T_pop_sim_simplex_l= {}".format(pprint.pformat(E_T_pop_sim_simplex_l) ) )
  plot.plot(arr_rate_l, E_T_pop_sim_simplex_l, label="Simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # log(WARNING, "E_T_unpop_sim_simplex_l= {}".format(pprint.pformat(E_T_unpop_sim_simplex_l) ) )
  # plot.plot(arr_rate_l, E_T_unpop_sim_simplex_l, label="Simulation, unpopular", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.plot(arr_rate_l, E_T_pop_ub_simplex_l, label="Upper bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  print("E_T_pop_lb_simplex_l= {}".format(pprint.pformat(E_T_pop_lb_simplex_l) ) )
  plot.plot(arr_rate_l, E_T_pop_lb_simplex_l, label="Lower bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'Arrival rate $\lambda$ (Request/sec)', fontsize=12)
  plot.ylabel(r'Average download time (sec)', fontsize=13)
  plot.title(r'Servers $\sim Exp(\mu={})$, availability $t={}$'.format(mu, t) )
  # plot.title(r'Servers $\sim$ Dolly, availability $t={}$'.format(t) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_ff_simplex_t_{}.pdf".format(t) )
  log(WARNING, "done; t= {}, r= {}, k= {}".format(t, r, k) )

if __name__ == "__main__":
  plot_ff_simplex()
  