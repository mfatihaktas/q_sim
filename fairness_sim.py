import numpy, pprint

from simplex_models import simplex_sym_l__sym__rgroup_l_m
from simplex_sim import *
from fairness_model import *
from patch import *

HOT_SYM = 'a'

class FF_JSink(object): # Join
  def __init__(self, _id, env):
    self._id = _id
    self.env = env
    
    self.st_hot_l = []
    self.st_cold_l = [] # system time
    self.qid__num_hot_win_map = {}
    
    self.store = simpy.Store(env)
    self.out = None
  
  def __repr__(self):
    return "JSink[_id= {}]".format(self._id)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    
    st = self.env.now - p.entrance_time
    if p.sym == HOT_SYM:
      self.st_hot_l.append(st)
      
      if p.winner_id not in self.qid__num_hot_win_map:
        self.qid__num_hot_win_map[p.winner_id] = 0
      self.qid__num_hot_win_map[p.winner_id] += 1
    else:
      self.st_cold_l.append(st)
    
    if self.out is not None:
      self.out.put(p)

# **********************************  Fairness First AVQ  ******************************* #
class FF_AVQ(object): # Fairness First
  def __init__(self, _id, env, t, sym__rgroup_l_map, serv, serv_dist_m, out=None):
    self._id = _id
    self.env = env
    self.sym__rgroup_l_map = sym__rgroup_l_map
    self.out = out
    
    self.sym_sysqid_map = {}
    for s, rg_l in sym__rgroup_l_map.items():
      sys_qid = None
      for g in rg_l:
        if len(g) == 1:
          sys_qid = g[0]
      self.sym_sysqid_map[s] = sys_qid
    
    self.num_q = int(1 + t*2)
    self.qid_l = [i for i in range(self.num_q) ]
    
    self.jsink = FF_JSink(_id, env)
    self.jsink.out = out
    
    self.jq = MT_AV_JQ(_id, env, self.qid_l, sym__rgroup_l_map)
    self.jq.out = self.jsink # data outlet
    self.jq.out_c = self # control outlet
    self.id_q_map = {}
    for i in self.qid_l:
      q = FCFS(i, env, serv, serv_dist_m)
      q.out = self.jq
      self.id_q_map[i] = q
    #
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    env.process(self.run() )
    env.process(self.run_c() )
    
    self.job_id_counter = 0
    self.starttype__num_map = {i:0 for i in range(t+1) }
    
    self.hot_store = simpy.Store(env)
    env.process(self.send_hot() )
    self.release_hot = None
    
    self.qid__job_id_l_map = {i:[] for i in range(self.num_q) }
    
    self.m = FF_AVQMonitor(env, self.id_q_map[1], poll_rate=1)
  
  def __repr__(self):
    return "FF_AVQ[qid_l= {}]".format(self.qid_l)
  
  def state(self):
    return {i:q.length() for i,q in self.id_q_map.items() }
  
  def send_hot(self):
    while True:
      p = (yield self.hot_store.get() )
      
      sys_qid = self.sym_sysqid_map[p.sym]
      if len(self.qid__job_id_l_map[sys_qid] ) > 0:
        self.release_hot = self.env.event()
        yield self.release_hot
      
      self.qid__job_id_l_map[sys_qid].append(p.job_id)
      self.id_q_map[sys_qid].put(p.deep_copy() )
      st = 0
      for r_l in self.sym__rgroup_l_map[HOT_SYM]:
        if len(r_l) == 1: continue
        
        rep = True
        for r in r_l:
          if len(self.qid__job_id_l_map[r] ) > 0:
            rep = False
            break
        if rep:
          for r in r_l:
            self.qid__job_id_l_map[r].append(p.job_id)
            self.id_q_map[r].put(p.deep_copy() )
          st += 1
      self.starttype__num_map[st] += 1
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      
      if p.sym == HOT_SYM:
        self.hot_store.put(p)
      else:
        sys_qid = self.sym_sysqid_map[p.sym]
        sys_q = self.id_q_map[sys_qid]
        
        for r_l in self.sym__rgroup_l_map[HOT_SYM]:
          if sys_qid in r_l:
            for r in r_l:
              self.id_q_map[r].put_c(CPacket(_id=-1, sym=HOT_SYM, prev_hop_id=self._id) )
              # print("sending cancel to server-{}".format(r) )
            break
        # To check if any hot data left
        # for p in sys_q.p_l:
        #   if p.sym == HOT_SYM:
        #     log(ERROR, "Hot data still in p= {}".format(p) )
        
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
      
      if cp.sym == HOT_SYM:
        for g, q in self.id_q_map.items():
          if q._id not in cp.departed_qid_l:
            q.put_c(cp.deep_copy() )
        
        if self.release_hot is not None:
          self.release_hot.succeed()
          self.release_hot = None
      
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
    self.num_hot_in_l = []
    
    env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(1/self.poll_rate)
      self.num_hot_in_l.append(self.q.num_sym_in(HOT_SYM) )

# *********************************  Fairness-First  *********************************** #
def test_ff_avq(num_f_run, t, cold_ar, hot_ar, serv, serv_dist_m, mds=False):
  E_T_hot_f_sum, E_T_cold_f_sum = 0, 0
  for f in range(num_f_run):
    log(WARNING, "t= {}, cold_ar= {}, hot_ar= {}, serv= {}, serv_dist_m= {}".format(t, cold_ar, hot_ar, serv, serv_dist_m) )
    
    (sym_l, sym__rgroup_l_map) = simplex_sym_l__sym__rgroup_l_m(t)
    log(WARNING, "sym__rgroup_l_map=\n {}".format(pprint.pformat(sym__rgroup_l_map) ) )
    
    env = simpy.Environment()
    pg = MT_PG(env, "pg", (len(sym_l)-1)*cold_ar, sym_l, HOT_SYM, hot_ar)
    avq = FF_AVQ("ff_avq", env, t, sym__rgroup_l_map, serv, serv_dist_m)
    pg.out = avq
    pg.init()
    env.run(until=2*50000)
    
    # for n in avq.m.num_hot_in_l:
    #   if n > 1:
    #     log(ERROR, "# of hot= {} > 1!".format(n) )
    
    l = avq.jsink.st_hot_l
    if len(l): E_T_hot_f_sum += float(sum(l) )/len(l)
    l = avq.jsink.st_cold_l
    if len(l): E_T_cold_f_sum += float(sum(l) )/len(l)
    
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_hot_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    print("pg.sym__n_sent= {}".format(pg.sym__n_sent) )
    
    print("avq.starttype__num_map= {}".format(pprint.pformat(avq.starttype__num_map) ) )
    
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_hot_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  E_T_hot = E_T_hot_f_sum/num_f_run
  E_T_cold = E_T_cold_f_sum/num_f_run
  print(">> E_T_hot= {}, E_T_cold= {}".format(E_T_hot, E_T_cold) )
  return (E_T_cold, E_T_hot)

def plot_ff_simplex():
  t = 1
  serv = "Exp" # "Dolly"
  mu = 1
  cold_ar = 0.1 # 0.5
  if serv == "Exp":
    serv_dist_m = {'mu': mu}
    ar_ub = hot_ar_ub_ff_simplex_approx(cold_ar, t, mu)
  elif serv == "Dolly":
    if t == 1: ar_ub = 0.28
    elif t == 3: ar_ub = 0.4
  log(WARNING, "t= {}, cold_ar= {}, serv= {}, serv_dist_m= {}, ar_ub= {}".format(t, cold_ar, serv, serv_dist_m, ar_ub) )
  
  E_T_cold_sim_simplex_l, E_T_hot_sim_simplex_l = [], []
  E_T_hot_ub_simplex_l, E_T_hot_lb_simplex_l, E_T_hot_approx_simplex_l = [], [], []
  
  num_f_run = 3
  sim_simplex = False
  if serv == "Exp":
    if t == 1:
      if cold_ar == 0.1:
        # num_f_run = 3
        E_T_hot_sim_simplex_l= [
          0.7508814818497344,
          0.8877425592134958,
          1.0935571693945558,
          1.4111299366247418,
          2.040702365903609,
          3.703025182692372,
          4.273110212627999,
          4.628558659826815,
          5.7066698290831654,
          7.160420826964484,
          9.640831852038255,
          13.791015734978886,
          21.9774519416767,
          None, # 72.98658763976105
          None] # 707.4689410004945
      elif cold_ar == 0.5:
        E_T_hot_sim_simplex_l= [
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
    elif t == 3:
      if cold_ar == 0.1:
        pass
      elif cold_ar == 0.5:
        E_T_hot_sim_simplex_l= [
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
          None, # 22.48950553473519,
          None, # 97.33546499091221,
          None] # 1274.0523614601243
    elif t == 77:
      pass
    else: sim_simplex = True
  elif serv == "Dolly":
    if t == 11:
      pass
    elif t == 33:
      pass
    else: sim_simplex = True
  else: sim_simplex = True
  
  sim_mds = False
  if serv == "Exp":
    if t == 11:
      pass
    elif t == 33:
      pass
    else: sim_mds = True
  
  mew, ms = 2, 8
  ar_l = []
  for hot_ar in [*numpy.linspace(0.05, 0.8*ar_ub, 5, endpoint=False), *numpy.linspace(0.8*ar_ub, ar_ub, 10) ]:
  # for hot_ar in numpy.linspace(0.1, 0.1, 1):
    ar_l.append(hot_ar)
    if sim_simplex:
      (E_T_cold, E_T_hot) = test_ff_avq(num_f_run, t, cold_ar, hot_ar, serv, serv_dist_m)
      E_T_cold_sim_simplex_l.append(E_T_cold)
      E_T_hot_sim_simplex_l.append(E_T_hot)
    E_T_hot_approx_simplex_l.append(E_T_hot_approx_ff_simplex(hot_ar, cold_ar, t, serv, serv_dist_m) )
  log(WARNING, "E_T_hot_sim_simplex_l= {}".format(pprint.pformat(E_T_hot_sim_simplex_l) ) )
  plot.plot(ar_l, E_T_hot_sim_simplex_l, label="Simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(ar_l, E_T_hot_approx_simplex_l, label="Approximation", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  ar_ub_min = hot_ar_ub_min_ff_simplex(t, mu)
  for hot_ar in [*numpy.linspace(0.05, 0.8*ar_ub_min, 5, endpoint=False), *numpy.linspace(0.8*ar_ub_min, ar_ub_min, 10) ]:
    ar_l.append(hot_ar)
    E_T_hot_ub_simplex_l.append(E_T_hot_ub_ff_simplex(hot_ar, t, mu) )
  plot.plot(ar_l, E_T_hot_ub_simplex_l, label="Upper bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  ar_ub_max = hot_ar_ub_max_ff_simplex(t, mu)
  for hot_ar in [*numpy.linspace(0.05, 0.8*ar_ub_max, 5, endpoint=False), *numpy.linspace(0.8*ar_ub_max, ar_ub_max, 10) ]:
    ar_l.append(hot_ar)
    E_T_hot_lb_simplex_l.append(E_T_hot_lb_ff_simplex(hot_ar, t, mu) )
  plot.plot(ar_l, E_T_hot_lb_simplex_l, label="Lower bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  # log(WARNING, "E_T_cold_sim_simplex_l= {}".format(pprint.pformat(E_T_cold_sim_simplex_l) ) )
  # plot.plot(ar_l, E_T_cold_sim_simplex_l, label="Simulation, cold data", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  # ar_l = []
  # mds_ar_ub = ar_ub
  # for hot_ar in [*numpy.linspace(0.05, 0.8*mds_ar_ub, 5, endpoint=False), *numpy.linspace(0.8*mds_ar_ub, mds_ar_ub, 10) ]:
  #   ar_l.append(hot_ar)
  #   if sim_mds:
  #     (E_T_cold, E_T_hot) = test_ff_avq(num_f_run, t, cold_ar, hot_ar, serv, serv_dist_m, mds=True)
  #     E_T_cold_sim_mds_l.append(E_T_cold)
  #     E_T_hot_sim_mds_l.append(E_T_hot)
  #   # E_T_hot_approx_mds_l.append(E_T_hot_approx_ff_mds(hot_ar, cold_ar, t, mu) )
  # log(WARNING, "E_T_hot_sim_mds_l= {}".format(pprint.pformat(E_T_hot_sim_mds_l) ) )
  # plot.plot(ar_l, E_T_hot_sim_mds_l, label="Simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'Hot data arrival rate $\lambda$ (Request/s)', fontsize=12)
  plot.ylabel(r'Average hot data download time (s)', fontsize=12)
  plot.title(r'$S \sim Exp(\mu={})$, $t={}$, $\lambda_c={}$'.format(mu, t, cold_ar) )
  # plot.title(r'Servers $\sim$ Dolly, availability $t={}$'.format(t) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.3)
  fig.tight_layout()
  plot.savefig("plot_ff_simplex_t_{}_lc_{}.pdf".format(t, cold_ar) )
  log(WARNING, "done; t= {}".format(t) )

# *********************************  Rep-to-all  *********************************** #
def test_simplex_reptoall(num_f_run, t, cold_ar, hot_ar, serv, serv_dist_m):
  E_T_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "t= {}, cold_ar= {}, hot_ar= {}, serv= {}, serv_dist_m= {}, ".format(t, cold_ar, hot_ar, serv, serv_dist_m) )
    
    (sym_l, sym__rgroup_l_map) = simplex_sym_l__sym__rgroup_l_m(t)
    log(WARNING, "sym__rgroup_l_map=\n {}".format(pprint.pformat(sym__rgroup_l_map) ) )
    
    env = simpy.Environment()
    pg = MT_PG(env, "pg", (len(sym_l)-1)*cold_ar, sym_l, HOT_SYM, hot_ar)
    avq = MT_AVQ("mt_avq", env, t, sym__rgroup_l_map, serv, serv_dist_m)
    pg.out = avq
    pg.init()
    c = 4 if serv == "Pareto" else 1
    env.run(until=c*50000) # 20
    
    # print("pg.sym__n_sent= {}".format(pprint.pformat(pg.sym__n_sent) ) )
    st_l = avq.jsink.st_l
    if len(st_l) > 0:
      E_T_f_sum += float(sum(st_l) )/len(st_l)
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  E_T = E_T_f_sum/num_f_run
  print(">> E_T= {}".format(E_T) )
  if E_T > 100: return None
  return E_T

def plot_reptoall_over_ff():
  t = 3
  serv = "Exp" # "Dolly"
  mu = 1
  cold_ar = 0.5 # 0.1 # 0.5
  if serv == "Exp":
    serv_dist_m = {'mu': mu}
    ar_ub = hot_ar_ub_ff_simplex_approx(cold_ar, t, mu)
  elif serv == "Dolly":
    if t == 1: ar_ub = 0.28
    elif t == 3: ar_ub = 0.4
  log(WARNING, "t= {}, cold_ar= {}, serv= {}, serv_dist_m= {}, ar_ub= {}".format(t, cold_ar, serv, serv_dist_m, ar_ub) )
  
  E_T_reptoall_cold_sim_l, E_T_reptoall_hotcold_sim_l = [], []
  E_T_ff_cold_sim_l, E_T_ff_hot_sim_l = [], []
  
  num_f_run = 3
  sim = False
  if serv == "Exp":
    if t == 1:
      if cold_ar == 0.1:
        E_T_reptoall_hotcold_sim_l= [
          0.7201505037958582,
          0.8137041264569737,
          0.9576988819341777,
          1.1556965248077378,
          1.4620610157759621,
          2.1137845402305997,
          2.244746996515651,
          2.371611126513781,
          2.5237271493291082,
          2.8241953240042132,
          2.9999307033411413,
          3.435526854246021,
          3.8673792894924,
          4.502426067143943,
          5.060620345144394]
        E_T_ff_hot_sim_l= [
          0.7528640774127929,
          0.8883769671299149,
          1.0910725127844978,
          1.411467908695325,
          2.0765024736783175,
          3.7413272910643443,
          4.429733209072531,
          4.700459068818758,
          5.766598589325757,
          7.562666174028496,
          9.230175145500544,
          13.539093729563461,
          20.50708165354279,
          64.65613140956314,
          None] # 632.4531468836764
        E_T_ff_cold_sim_l= [
          1.1130973761654523,
          1.1163667464165845,
          1.1251563242508384,
          1.1088422654185905,
          1.1131072736103387,
          1.1221882386517144,
          1.10861950300541,
          1.1162110193221944,
          1.105966514755852,
          1.1230458893264734,
          1.1077928002418798,
          1.1069175571746173,
          1.1132450899382687,
          1.1291230722102326,
          1.1138810147686706]
      elif cold_ar == 0.5:
        # num_f_run = 3
        E_T_reptoall_hotcold_sim_l= [
          0.9363139990245205,
          1.0592291690482762,
          1.231488381151278,
          1.5228553548169674,
          2.0170475230291993,
          3.1620874832406725,
          3.473205375661949,
          3.866809382647841,
          4.381243387049081,
          4.842042025556865,
          5.664916818649885,
          6.9026595925007515,
          7.719995100957113,
          10.64296466417875,
          15.515612762889816]
        E_T_ff_hot_sim_l= [
          0.9319581938424295,
          1.115623991352767,
          1.4227089127566297,
          1.8722202284409732,
          2.795986722212484,
          5.505853558891492,
          6.139470869379135,
          7.536927829171991,
          9.47190831179973,
          12.762001358056489,
          16.063831672136484,
          24.901178025970523,
          181.76685991791166,
          None, # 872.4424105424245
          None] # 1647.3312422694653
        E_T_ff_cold_sim_l= [
          1.977366720016086,
          1.9922420452600413,
          1.9636128875179872,
          2.0094453457891093,
          2.007427338320136,
          2.022709999075777,
          1.9970721345210727,
          1.9747366413959602,
          2.006676286948779,
          2.0017823160561843,
          2.0122851716339354,
          2.0141071336728356,
          1.9805536042871632,
          1.9915303364254016,
          1.9846564963461162]
    elif t == 33:
      if cold_ar == 0.1:
        pass
      elif cold_ar == 0.5:
        pass
    elif t == 77:
      pass
    else: sim = True
  else: sim = True
  
  mew, ms = 2, 8
  ar_l = []
  for hot_ar in [*numpy.linspace(0.05, 0.8*ar_ub, 5, endpoint=False), *numpy.linspace(0.8*ar_ub, ar_ub, 10) ]:
    ar_l.append(hot_ar)
    if sim:
      E_T = test_simplex_reptoall(num_f_run, t, cold_ar, hot_ar, serv, serv_dist_m)
      E_T_reptoall_cold_sim_l.append(E_T)
      E_T_reptoall_hotcold_sim_l.append(E_T)
      
      (E_T_cold, E_T_hot) = test_ff_avq(num_f_run, t, cold_ar, hot_ar, serv, serv_dist_m)
      E_T_ff_cold_sim_l.append(E_T_cold)
      E_T_ff_hot_sim_l.append(E_T_hot)
  log(WARNING, "E_T_reptoall_hotcold_sim_l= {}".format(pprint.pformat(E_T_reptoall_hotcold_sim_l) ) )
  log(WARNING, "E_T_ff_hot_sim_l= {}".format(pprint.pformat(E_T_ff_hot_sim_l) ) )
  log(WARNING, "E_T_ff_cold_sim_l= {}".format(pprint.pformat(E_T_ff_cold_sim_l) ) )
  
  gain_l, pain_l = [], []
  for i, E_T_reptoall in enumerate(E_T_reptoall_hotcold_sim_l):
    if E_T_ff_hot_sim_l[i] is None:
      gain_in_hot = None
      pain_in_cold = None
    else:
      gain_in_hot = (E_T_ff_hot_sim_l[i] - E_T_reptoall)/E_T_ff_hot_sim_l[i] * 100
      pain_in_cold = (E_T_reptoall - E_T_ff_cold_sim_l[i])/E_T_ff_cold_sim_l[i] * 100
    gain_l.append(gain_in_hot)
    pain_l.append(pain_in_cold)
  plot.plot(ar_l, gain_l, label="Gain in hot", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(ar_l, pain_l, label="Pain in cold", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'Hot data arrival rate $\lambda$ (Request/s)', fontsize=12)
  plot.ylabel(r'Pain or gain in percentange', fontsize=12)
  plot.title('Pain and gain of Rep-to-all over Fairness-First\n' + r'$S \sim Exp(\mu={})$, $t={}$, $\lambda_c={}$'.format(mu, t, cold_ar) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_reptoall_over_ff_t_{}_lc_{}.pdf".format(t, cold_ar) )
  log(WARNING, "done; t= {}".format(t) )

if __name__ == "__main__":
  # plot_ff_simplex()
  plot_reptoall_over_ff()
