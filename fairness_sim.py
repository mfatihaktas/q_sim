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
    self.qid__num_pop_win_map = {}
    
    self.store = simpy.Store(env)
    self.out = None
  
  def __repr__(self):
    return "JSink[_id= {}]".format(self._id)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    
    st = self.env.now - p.entrance_time
    if p.sym == POP_SYM:
      self.st_pop_l.append(st)
      
      if p.winner_id not in self.qid__num_pop_win_map:
        self.qid__num_pop_win_map[p.winner_id] = 0
      self.qid__num_pop_win_map[p.winner_id] += 1
    else:
      self.st_unpop_l.append(st)
    
    if self.out is not None:
      self.out.put(p)

# **********************************  Fairness First AVQ  ******************************* #
class FF_AVQ(object): # Fairness First
  def __init__(self, _id, env, t, qmu_l, serv, sym__rgroup_l_map, sym_sysqid_map, out=None):
    self._id = _id
    self.env = env
    self.sym__rgroup_l_map = sym__rgroup_l_map
    self.sym_sysqid_map = sym_sysqid_map
    self.out = out
    
    self.num_q = len(qmu_l)
    self.qid_l = [i for i in range(self.num_q) ]
    
    self.jsink = FF_JSink(_id, env)
    self.jsink.out = out
    
    self.jq = MT_AV_JQ(_id, env, self.qid_l, sym__rgroup_l_map)
    self.jq.out = self.jsink # data outlet
    self.jq.out_c = self # control outlet
    self.id_q_map = {}
    for i, qid in enumerate(self.qid_l):
      q = FCFS(qid, env, qmu_l[i], serv)
      q.out = self.jq
      self.id_q_map[qid] = q
    #
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    env.process(self.run() )
    env.process(self.run_c() )
    
    self.job_id_counter = 0
    self.starttype__num_map = {i:0 for i in range(t+1) }
    
    self.pop_store = simpy.Store(env)
    env.process(self.send_pop() )
    self.release_pop = None
    
    self.qid__job_id_l_map = {i:[] for i in range(self.num_q) }
    
    self.m = FF_AVQMonitor(env, self.id_q_map[1], poll_rate=1)
  
  def __repr__(self):
    return "FF_AVQ[qid_l= {}]".format(self.qid_l)
  
  def state(self):
    return {i:q.length() for i,q in self.id_q_map.items() }
  
  def send_pop(self):
    while True:
      p = (yield self.pop_store.get() )
      
      sys_qid = self.sym_sysqid_map[p.sym]
      if len(self.qid__job_id_l_map[sys_qid] ) > 0:
        self.release_pop = self.env.event()
        yield self.release_pop
      
      self.id_q_map[sys_qid].put(p.deep_copy() )
      self.qid__job_id_l_map[sys_qid].append(p.job_id)
      st = 0
      for r_l in self.sym__rgroup_l_map[POP_SYM]:
        if len(r_l) == 1: continue
        
        rep = True
        for r in r_l:
          # q = self.id_q_map[r]
          # if q.length() == 1:
          #   # q.p_in_serv is not None and
          #   if q.p_in_serv.sym != POP_SYM: # q may be busy with pop_sym because simpy did not register cancellation yet
          #     print("server-{} busy with sym= {}, while sending p= {}, pstate_l= {}".format(r, q.p_in_serv.sym, p, q.pstate_l() ) )
          #     rep = False
          #     break
          # elif q.busy():
          #   print("server-{} busy with sym= {}, while sending p= {}, q.length= {}, pstate_l= {}".format(r, q.p_in_serv.sym, p, q.length(), q.pstate_l() ) )
          #   rep = False
          #   break
          if len(self.qid__job_id_l_map[r] ) > 0:
            rep = False
            break
        if rep:
          for r in r_l:
            self.id_q_map[r].put(p.deep_copy() )
            self.qid__job_id_l_map[r].append(p.job_id)
          st += 1
      self.starttype__num_map[st] += 1
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      
      if p.sym == POP_SYM:
        self.pop_store.put(p)
      else:
        sys_qid = self.sym_sysqid_map[p.sym]
        sys_q = self.id_q_map[sys_qid]
        
        for r_l in self.sym__rgroup_l_map[POP_SYM]:
          if sys_qid in r_l:
            for r in r_l:
              self.id_q_map[r].put_c(CPacket(_id=p.job_id, sym=POP_SYM, prev_hop_id=self._id) )
              # print("sending cancel to server-{}".format(r) )
            break
        sys_q.put(p.deep_copy() )
        self.qid__job_id_l_map[sys_qid].append(p.job_id)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.job_id_counter += 1
    p.job_id = self.job_id_counter
    return self.store.put(p.deep_copy() )
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      # sim_log(WARNING, self.env, self, "recved", cp)
      
      for g, q in self.id_q_map.items():
        if q._id not in cp.departed_qid_l:
          q.put_c(cp.deep_copy() )
      
      if cp.sym == POP_SYM:
        if self.release_pop is not None:
          self.release_pop.succeed()
          self.release_pop = None
      
      for i, job_id_l in self.qid__job_id_l_map.items():
        try: job_id_l.remove(cp._id)
        except: pass
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp.deep_copy() )

class FF_AVQMonitor(object):
  def __init__(self, env, q, poll_rate):
    self.q = q
    self.env = env
    self.poll_rate = poll_rate
    
    # self.t_l = [] # Time steps that the numbers polled from the q
    self.num_pop_in_l = []
    
    env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(1/self.poll_rate)
      self.num_pop_in_l.append(self.q.num_sym_in(POP_SYM) )

# *********************************  Experiment  *********************************** #
def test_ff_avq(num_f_run, unpop_ar, mu, t, r, k, serv, pop_ar, mds=False):
  num_q = int(1 + t*r)
  qmu_l = num_q*[mu]
  
  E_T_pop_f_sum, E_T_unpop_f_sum = 0, 0
  for f in range(num_f_run):
    log(WARNING, "unpop_ar= {}, pop_ar= {}, mu= {}, t= {}, r= {}, k= {}, qmu_l= {}".format(unpop_ar, pop_ar, mu, t, r, k, qmu_l) )
    env = simpy.Environment()
    
    sym_sysqid_map, sym__rgroup_l_map = {}, {}
    if t == 1: sym_l = ['a', 'b']
    elif t == 3: sym_l = ['a', 'b', 'c']
    for sym in sym_l:
      rgroup_l = []
      if t == 1:
        if sym == 'a':
          sym_sysqid_map[sym] = 0
          rgroup_l.append([0] )
          rgroup_l.append([1, 2] )
        elif sym == 'b':
          sym_sysqid_map[sym] = 1
          rgroup_l.append([1] )
          rgroup_l.append([0, 2] )
      elif t == 3:
        if not mds:
          if sym == 'a':
            sym_sysqid_map[sym] = 0
            rgroup_l.append([0] )
            rgroup_l.append([1, 2] )
            rgroup_l.append([3, 4] )
            rgroup_l.append([5, 6] )
          elif sym == 'b':
            sym_sysqid_map[sym] = 1
            rgroup_l.append([1] )
            rgroup_l.append([0, 2] )
            rgroup_l.append([3, 5] )
            rgroup_l.append([4, 6] )
          elif sym == 'c':
            sym_sysqid_map[sym] = 3
            rgroup_l.append([3] )
            rgroup_l.append([0, 4] )
            rgroup_l.append([1, 5] )
            rgroup_l.append([2, 6] )
        elif mds:
          if sym == 'a':
            sym_sysqid_map[sym] = 0
            rgroup_l.append([0] )
            rgroup_l.append([1, 2, 3] )
            rgroup_l.append([4, 5, 6] )
          elif sym == 'b':
            sym_sysqid_map[sym] = 1
            rgroup_l.append([1] )
          elif sym == 'c':
            sym_sysqid_map[sym] = 2
            rgroup_l.append([2] )
      sym__rgroup_l_map[sym] = rgroup_l
    pg = MT_PG(env, "pg", (len(sym_l)-1)*unpop_ar, sym_l, pop_sym=POP_SYM, pop_ar=pop_ar)
    log(WARNING, "sym__rgroup_l_map=\n {}\n sym_sysqid_map=\n {}".format(pprint.pformat(sym__rgroup_l_map), pprint.pformat(sym_sysqid_map) ) )
    avq = FF_AVQ("ff_avq", env, t, qmu_l, serv, sym__rgroup_l_map, sym_sysqid_map)
    pg.out = avq
    pg.init()
    env.run(until=2*50000) # 50000
    
    # for n in avq.m.num_pop_in_l:
    #   if n > 1:
    #     log(ERROR, "# of pop= {} > 1!".format(n) )
    
    l = avq.jsink.st_pop_l
    # print("avq.jsink.st_pop_l= {}".format(l) )
    if len(l): E_T_pop_f_sum += float(sum(l) )/len(l)
    l = avq.jsink.st_unpop_l
    if len(l): E_T_unpop_f_sum += float(sum(l) )/len(l)
    
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_pop_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    print("pg.sym__n_sent= {}".format(pg.sym__n_sent) )
    
    print("avq.starttype__num_map= {}".format(pprint.pformat(avq.starttype__num_map) ) )
    
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_pop_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  E_T_pop = E_T_pop_f_sum/num_f_run
  E_T_unpop = E_T_unpop_f_sum/num_f_run
  print(">> E_T_pop= {}, E_T_unpop= {}".format(E_T_pop, E_T_unpop) )
  # if E_T_unpop > E_T_UB: return (None, None)
  return (E_T_pop, E_T_unpop)
  
def plot_ff_simplex():
  t, r, k = 3, 2, 2
  serv = "Exp" # "Dolly"
  mu = 1
  unpop_ar = 0.5 # 0.05 # 0.5 # 0.9
  if serv == "Exp":
    ar_ub = pop_ar_ub_ff_simplex_approx(unpop_ar, t, mu)
  elif serv == "Dolly":
    if t == 1: ar_ub = 0.28
    elif t == 3: ar_ub = 0.4
  log(WARNING, "t= {}, r= {}, k= {}, mu= {}, serv= {}, unpop_ar= {}, ar_ub= {}".format(t, r, k, mu, serv, unpop_ar, ar_ub) )
  
  E_T_unpop_sim_simplex_l, E_T_pop_sim_simplex_l = [], []
  E_T_pop_ub_simplex_l, E_T_pop_lb_simplex_l, E_T_pop_approx_simplex_l = [], [], []
  
  sim_simplex = False
  if serv == "Exp":
    if t == 1:
      if unpop_ar == 0.0001:
        E_T_pop_sim_simplex_l= [
          0.6860751158403763,
          0.8005599202866529,
          0.9723723138044678,
          1.2166101108790657,
          1.7216640731447055,
          2.8736382051836835,
          3.1840946825045027,
          3.7291881704421277,
          4.261254826933591,
          4.91211849239902,
          6.471718165339034,
          9.07663506621845,
          11.229855186330713,
          16.227913935225132,
          55.58059597074051]
      elif unpop_ar == 0.05:
        E_T_pop_sim_simplex_l= [
          0.7327683934087881,
          0.8560955329404556,
          1.0270311319092937,
          1.2910150057416987,
          1.8762614517674574,
          3.2613156631457647,
          3.5816483537239643,
          4.111518589016356,
          4.83902479158454,
          6.057062700103652,
          8.045091092797048,
          11.214589660849484,
          14.586167287114673,
          33.44738505096981,
          None]
      elif unpop_ar == 0.5:
        E_T_pop_sim_simplex_l= [
          0.9355006832487013,
          1.1208049761868801,
          1.3893483102228363,
          1.8971749729934722,
          2.7524297546978698,
          5.5850641343205005,
          6.067310795156397,
          7.208867107864213,
          9.197031614327727,
          12.038706806221821,
          18.5928197729436,
          33.37017843702431,
          None, # 73.25485711218015,
          None,
          None]
      elif unpop_ar == 0.9:
        E_T_pop_sim_simplex_l= [
          1.0359952211542578,
          1.2442409356984336,
          1.4932059436942178,
          1.9499350937409239,
          2.874305826627953,
          5.0719054300074085,
          5.89966254880032,
          6.705821426225249,
          7.594024097153369,
          9.117411167967962,
          11.701326304349212,
          14.160534580060528,
          23.040827931310787,
          34.05609250262459,
          None]
    elif t == 3:
      if unpop_ar == 0.05:
        E_T_pop_sim_simplex_l= [
          0.474682115366536,
          0.5513283357182526,
          0.6616703388858086,
          0.8347162899985489,
          1.1783099870860145,
          2.000525389700597,
          2.272541259951746,
          2.552886411822349,
          2.9449171918650148,
          3.2821169148860987,
          3.932213628674842,
          5.823877483207954,
          10.17432639498029,
          15.792913977390553,
          None]
          # 115.26928517746343]
      elif unpop_ar == 0.5:
        E_T_pop_sim_simplex_l= [
          0.589749465116383,
          0.6946720057825672,
          0.8451953891228225,
          1.0929876952490964,
          1.5946791649136038,
          3.13546222142155,
          3.337230268177806,
          3.879738208621262,
          5.091169628769583,
          5.804356046565282,
          8.503231612031973,
          15.717938719789103,
          22.48950553473519,
          None,
          None]
          # 97.33546499091221,
          # 1274.0523614601243]
      elif unpop_ar == 0.9:
        E_T_pop_sim_simplex_l= [
          0.6854469411718157,
          0.7775158117765747,
          0.9530942589412483,
          1.2234353036129992,
          1.6963500873668853,
          3.0641066301181565,
          3.5851219129343748,
          3.8900424463235637,
          4.724835883838208,
          5.850094918231982,
          7.09750223966705,
          10.454638502172712,
          13.540757905874694,
          26.460116374413456,
          None]
          # 216.13373742455795]
    elif t == 77:
      pass
    else: sim_simplex = True
  if serv == "Dolly":
    if t == 11:
      pass
    elif t == 33:
      pass
    else: sim_simplex = True
  
  sim_mds = False
  if serv == "Exp":
    if t == 11:
      pass
    elif t == 33:
      pass
    else: sim_mds = True
  
  mew, ms = 2, 8
  num_f_run = 1
  ar_l = []
  for pop_ar in [*numpy.linspace(0.05, 0.8*ar_ub, 5, endpoint=False), *numpy.linspace(0.8*ar_ub, ar_ub, 10) ]:
  # for pop_ar in numpy.linspace(0.0001, 0.0001, 1):
    ar_l.append(pop_ar)
    if sim_simplex:
      (E_T_pop, E_T_unpop) = test_ff_avq(num_f_run, unpop_ar, mu, t, r, k, serv, pop_ar=pop_ar)
      E_T_pop_sim_simplex_l.append(E_T_pop)
      E_T_unpop_sim_simplex_l.append(E_T_unpop)
    E_T_pop_approx_simplex_l.append(E_T_pop_approx_ff_simplex(pop_ar, unpop_ar, t, mu) )
  log(WARNING, "E_T_pop_sim_simplex_l= {}".format(pprint.pformat(E_T_pop_sim_simplex_l) ) )
  plot.plot(ar_l, E_T_pop_sim_simplex_l, label="Simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(ar_l, E_T_pop_approx_simplex_l, label="Approximation", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  ar_ub_min = pop_ar_ub_min_ff_simplex(t, mu)
  for pop_ar in [*numpy.linspace(0.05, 0.8*ar_ub_min, 5, endpoint=False), *numpy.linspace(0.8*ar_ub_min, ar_ub_min, 10) ]:
    ar_l.append(pop_ar)
    E_T_pop_ub_simplex_l.append(E_T_pop_ub_ff_simplex(pop_ar, t, mu) )
  plot.plot(ar_l, E_T_pop_ub_simplex_l, label="Upper bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  ar_ub_max = pop_ar_ub_max_ff_simplex(t, mu)
  for pop_ar in [*numpy.linspace(0.05, 0.8*ar_ub_max, 5, endpoint=False), *numpy.linspace(0.8*ar_ub_max, ar_ub_max, 10) ]:
    ar_l.append(pop_ar)
    E_T_pop_lb_simplex_l.append(E_T_pop_lb_ff_simplex(pop_ar, t, mu) )
  plot.plot(ar_l, E_T_pop_lb_simplex_l, label="Lower bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  # log(WARNING, "E_T_unpop_sim_simplex_l= {}".format(pprint.pformat(E_T_unpop_sim_simplex_l) ) )
  # plot.plot(ar_l, E_T_unpop_sim_simplex_l, label="Simulation, unpopular", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  mds_ar_ub = ar_ub
  for pop_ar in [*numpy.linspace(0.05, 0.8*mds_ar_ub, 5, endpoint=False), *numpy.linspace(0.8*mds_ar_ub, mds_ar_ub, 10) ]:
    ar_l.append(pop_ar)
    if sim_mds:
      (E_T_pop, E_T_unpop) = test_ff_avq(num_f_run, unpop_ar, mu, t, r, k, serv, pop_ar=pop_ar, mds=True)
      E_T_pop_sim_mds_l.append(E_T_pop)
      E_T_unpop_sim_mds_l.append(E_T_unpop)
    # E_T_pop_approx_mds_l.append(E_T_pop_approx_ff_mds(pop_ar, unpop_ar, t, mu) )
  log(WARNING, "E_T_pop_sim_mds_l= {}".format(pprint.pformat(E_T_pop_sim_mds_l) ) )
  plot.plot(ar_l, E_T_pop_sim_mds_l, label="Simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  # plot.plot(ar_l, E_T_pop_approx_simplex_l, label="Approximation", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'$\lambda_p$ (Request/sec)', fontsize=12)
  plot.ylabel(r'Average download time (sec)', fontsize=13)
  plot.title(r'$X \sim Exp(\mu={})$, $t={}$, $\lambda_u={}$'.format(mu, t, unpop_ar) )
  # plot.title(r'Servers $\sim$ Dolly, availability $t={}$'.format(t) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_ff_simplex_t_{}_ls_{}.pdf".format(t, unpop_ar) )
  log(WARNING, "done; t= {}, r= {}, k= {}".format(t, r, k) )

if __name__ == "__main__":
  plot_ff_simplex()
  