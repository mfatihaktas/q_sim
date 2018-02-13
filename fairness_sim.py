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
  def __init__(self, _id, env, t, sym__rgroup_l_m, sdist_m, out=None):
    self._id = _id
    self.env = env
    self.sym__rgroup_l_m = sym__rgroup_l_m
    self.out = out
    
    self.sym_sysqid_map = {}
    for s, rg_l in sym__rgroup_l_m.items():
      sys_qid = None
      for g in rg_l:
        if len(g) == 1:
          sys_qid = g[0]
      self.sym_sysqid_map[s] = sys_qid
    
    self.num_q = int(1 + t*2)
    self.qid_l = [i for i in range(self.num_q) ]
    
    self.jsink = FF_JSink(_id, env)
    self.jsink.out = out
    
    self.jq = MT_AV_JQ(_id, env, self.qid_l, sym__rgroup_l_m)
    self.jq.out = self.jsink # data outlet
    self.jq.out_c = self # control outlet
    self.id_q_map = {}
    for i in self.qid_l:
      q = FCFS(i, env, sdist_m)
      q.out = self.jq
      self.id_q_map[i] = q
    #
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    env.process(self.run() )
    # env.process(self.run_c() )
    
    self.jid_counter = 0
    self.starttype_numj_map = {i: 0 for i in range(t+1) }
    
    self.hot_store = simpy.Store(env)
    env.process(self.send_hot() )
    self.release_hot = None
    
    # self.m = FF_AVQMonitor(env, self.id_q_map[1], poll_rate=1)
  
  def __repr__(self):
    return "FF_AVQ[qid_l= {}]".format(self.qid_l)
  
  def state(self):
    return {i:q.length() for i,q in self.id_q_map.items() }
  
  def send_hot(self):
    while True:
      p = (yield self.hot_store.get() )
      
      sys_qid = self.sym_sysqid_map[p.sym]
      sys_q = self.id_q_map[sys_qid]
      if sys_q.busy():
        self.release_hot = self.env.event()
        yield self.release_hot
      self.release_hot = None
      sys_q.put(p.deep_copy() )
      # yield self.env.timeout(0.0001)
      
      st = 0
      for r_l in self.sym__rgroup_l_m[HOT_SYM]:
        if len(r_l) == 1: continue
        
        rep = True
        for r in r_l:
          if self.id_q_map[r].busy():
            rep = False
            q = self.id_q_map[r]
            sim_log(WARNING, self.env, self, "busy q!", p)
            print("q= {}\n q.p_inserv= {}, q.p_l= {}".format(q, q.p_inserv, q.p_l) )
            # print("q.p_inserv.job_id in q.recved_canceljid_l= {}".format(q.p_inserv.job_id in q.recved_canceljid_l) )
            # print("q.recved_canceljid_l= {}".format(q.recved_canceljid_l) )
            break
        if rep:
          for r in r_l:
            self.id_q_map[r].put(p.deep_copy() )
          st += 1
      self.starttype_numj_map[st] += 1
  
  def run(self):
    while True:
      p = yield self.store.get()
      
      if p.sym == HOT_SYM:
        self.hot_store.put(p)
      else:
        sys_qid = self.sym_sysqid_map[p.sym]
        sys_q = self.id_q_map[sys_qid]
        if sys_q.busy():
          sys_q.put_c(CPacket(_id=-1, sym=HOT_SYM, prev_hop_id=self._id) )
        
        # for r_l in self.sym__rgroup_l_m[HOT_SYM]:
        #   if sys_qid in r_l:
        #     for r in r_l:
        #       self.id_q_map[r].put_c(CPacket(_id=-1, sym=HOT_SYM, prev_hop_id=self._id) )
        #       # print("sending cancel to server-{}".format(r) )
        #     break
        
        # To check if any hot data left
        # for p in sys_q.p_l:
        #   if p.sym == HOT_SYM:
        #     log(ERROR, "Hot data still in p= {}".format(p) )
        
        sys_q.put(p.deep_copy() )
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.jid_counter += 1
    p.job_id = self.jid_counter
    return self.store.put(p.deep_copy() )
  
  # def run_c(self):
  #   while True:
  #     cp = (yield self.store_c.get() )
  #     # sim_log(WARNING, self.env, self, "recved", cp)
      
  #     if cp.sym == HOT_SYM:
  #       for _, q in self.id_q_map.items():
  #         if q._id not in cp.departed_qid_l:
  #           q.put_c(cp.deep_copy() )
        
  #       if self.release_hot is not None:
  #         self.release_hot.succeed()
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    # return self.store_c.put(cp.deep_copy() )
    
    if cp.sym == HOT_SYM:
      for _, q in self.id_q_map.items():
        if q._id not in cp.departed_qid_l:
          # q.put_c(cp.deep_copy() )
          q.put_c(cp)
      
      if self.release_hot is not None:
        self.release_hot.succeed()

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
def sim_ff(num_srun, t, car, har, sdist_m):
  ETh_sum, ETc_sum = 0, 0
  for _ in range(num_srun):
    log(WARNING, "t= {}, car= {}, har= {}, sdist_m= {}".format(t, car, har, sdist_m) )
    
    (sym_l, sym__rgroup_l_m) = simplex_sym_l__sym__rgroup_l_m(t)
    log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
    
    env = simpy.Environment()
    pg = MT_PG(env, "pg", (len(sym_l)-1)*car, sym_l, HOT_SYM, har)
    avq = FF_AVQ("ff_avq", env, t, sym__rgroup_l_m, sdist_m)
    pg.out = avq
    pg.init()
    env.run(until=2*50000)
    
    # for n in avq.m.num_hot_in_l:
    #   if n > 1:
    #     log(ERROR, "# of hot= {} > 1!".format(n) )
    
    l = avq.jsink.st_hot_l
    if len(l): ETh_sum += float(sum(l) )/len(l)
    l = avq.jsink.st_cold_l
    if len(l): ETc_sum += float(sum(l) )/len(l)
    
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_hot_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    print("pg.sym__n_sent= {}".format(pg.sym__n_sent) )
    
    print("avq.starttype_numj_map= {}".format(pprint.pformat(avq.starttype_numj_map) ) )
    
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_hot_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  ETh = ETh_sum/num_srun
  ETc = ETc_sum/num_srun
  print(">> ETh= {}, ETc= {}".format(ETh, ETc) )
  return (ETc, ETh)

def plot_fairnessfirst():
  t = 1
  serv = "Exp" # "Dolly"
  car = 0.000001 # 0.1 # 0.5
  if serv == "Exp":
    mu = 1
    sdist_m = {'dist': 'Exp', 'mu': mu}
    har_ub = ff_har_ub(car, t, sdist_m)
  elif serv == "Pareto":
    loc, a = 1, 3
    sdist_m = {'dist': 'Pareto', 'loc': loc, 'a': a}
    har_ub = ff_har_ub(car, t, sdist_m)
  elif serv == "Dolly":
    sdist_m = {'dist': 'Dolly'}
    if t == 1: har_ub = 0.28
    elif t == 3: har_ub = 0.4
  log(WARNING, "t= {}, car= {}, sdist_m= {}, har_ub= {}".format(t, car, serv, sdist_m, har_ub) )
  
  ETc_sim_l, ETh_sim_l = [], []
  ETh_ub_l, ETh_lb_l, ETh_approx_l = [], [], []
  
  num_srun = 1
  sim_simplex = False
  if serv == "Exp":
    if t == 11:
      if car == 0.1:
        # num_srun = 3
        ETh_sim_l= [
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
      elif car == 0.5:
        ETh_sim_l= [
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
      if car == 0.1:
        ETh_sim_l= [
          0.4939003286095989,
          0.5715098498251754,
          0.6839367636377401,
          0.8648712234306651,
          1.2162541140583187,
          2.1526032250152802,
          2.3332214430320124,
          2.77347691353428,
          3.1892878047378503,
          3.8146350092352512,
          5.053509482106505,
          6.571664953851065,
          10.8035629393868,
          25.898741689274527,
          None] # 484.0153086904779
      elif car == 0.5:
        ETh_sim_l= [
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
    else: sim_simplex = True
  elif serv == "Dolly":
    if t == 11:
      pass
    elif t == 33:
      pass
    else: sim_simplex = True
  else: sim_simplex = True
  
  mew, ms = 2, 6
  ar_l = []
  # for har in [*numpy.linspace(0.05, 0.8*har_ub, 5, endpoint=False), *numpy.linspace(0.8*har_ub, har_ub, 5) ]:
  for har in numpy.linspace(0.5, 0.5, 1):
    ar_l.append(har)
    if sim_simplex:
      (ETc, ETh) = sim_ff(num_srun, t, car, har, sdist_m)
      ETc_sim_l.append(ETc)
      ETh_sim_l.append(ETh)
    ETh_approx = ff_ETh_approx_(har, car, t, sdist_m)
    print("ETh_approx= {}".format(ETh_approx) )
    ETh_approx_l.append(ETh_approx)
  log(WARNING, "ETh_sim_l= {}".format(pprint.pformat(ETh_sim_l) ) )
  plot.plot(ar_l, ETh_sim_l, label="Simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(ar_l, ETh_approx_l, label="Approximation", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  har_ub_min = ff_har_ub_min(t, sdist_m)
  for har in [*numpy.linspace(0.05, 0.8*har_ub_min, 5, endpoint=False), *numpy.linspace(0.8*har_ub_min, har_ub_min, 10) ]:
    ar_l.append(har)
    ETh_ub_l.append(ff_ETh_ub(har, t, sdist_m) )
  plot.plot(ar_l, ETh_ub_l, label="Upper bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  har_ub_max = ff_har_ub_max(t, sdist_m)
  for har in [*numpy.linspace(0.05, 0.8*har_ub_max, 5, endpoint=False), *numpy.linspace(0.8*har_ub_max, har_ub_max, 10) ]:
    ar_l.append(har)
    ETh_lb_l.append(ff_ETh_lb(har, t, sdist_m) )
  plot.plot(ar_l, ETh_lb_l, label="Lower bound", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'Hot data arrival rate $\lambda$', fontsize=12)
  plot.ylabel(r'Avg hot data download time', fontsize=12)
  plot.title(r'$V \sim Exp(\mu={})$, $t={}$, $\lambda_c={}$'.format(mu, t, car) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.3)
  fig.tight_layout()
  plot.savefig("plot_fairnessfirst_t_{}_lc_{}.pdf".format(t, car) )
  log(WARNING, "done; t= {}".format(t) )

# *********************************  Rep-to-all  *********************************** #
def test_simplex_reptoall(num_srun, t, car, har, serv, sdist_m):
  ET_sum = 0
  for f in range(num_srun):
    log(WARNING, "t= {}, car= {}, har= {}, serv= {}, sdist_m= {}, ".format(t, car, har, serv, sdist_m) )
    
    (sym_l, sym__rgroup_l_m) = simplex_sym_l__sym__rgroup_l_m(t)
    log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
    
    env = simpy.Environment()
    pg = MT_PG(env, "pg", (len(sym_l)-1)*car, sym_l, HOT_SYM, har)
    avq = MT_AVQ("mt_avq", env, t, sym__rgroup_l_m, serv, sdist_m)
    pg.out = avq
    pg.init()
    c = 4 if serv == "Pareto" else 1
    env.run(until=c*50000)
    
    # print("pg.sym__n_sent= {}".format(pprint.pformat(pg.sym__n_sent) ) )
    st_l = avq.jsink.st_l
    if len(st_l) > 0:
      ET_sum += float(sum(st_l) )/len(st_l)
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  ET = ET_sum/num_srun
  print(">> ET= {}".format(ET) )
  if ET > 100: return None
  return ET

def plot_reptoall_over_ff():
  t = 3
  serv = "Dolly" # "Exp"
  mu = 1
  sdist_m = {'mu': mu}
  car = 0.5 # 0.1 # 0.5
  if serv == "Exp":
    har_ub = 0.9 * ff_har_ub(car, t, mu)
  elif serv == "Dolly":
    if t == 1: har_ub = 0.28
    elif t == 3: har_ub = 0.4
  log(WARNING, "t= {}, car= {}, serv= {}, sdist_m= {}, har_ub= {}".format(t, car, serv, sdist_m, har_ub) )
  
  ETc_reptoall_sim_l, EThc_reptoall_sim_l = [], []
  ETc_ff_sim_l, ETh_ff_sim_l = [], []
  
  num_srun = 3
  sim = False
  if serv == "Exp":
    if t == 1:
      if car == 0.1:
        EThc_reptoall_sim_l= [
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
        ETh_ff_sim_l= [
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
        ETc_ff_sim_l= [
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
      elif car == 0.5:
        # num_srun = 3
        EThc_reptoall_sim_l= [
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
        ETh_ff_sim_l= [
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
        ETc_ff_sim_l= [
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
    elif t == 3:
      if car == 0.1:
        pass
      elif car == 0.5:
        EThc_reptoall_sim_l= [
          0.6840010904077812,
          0.7859533196791282,
          0.9274142482662505,
          1.1707063798084094,
          1.6891848540864387,
          3.055811156306783,
          3.587534635804969,
          4.143524780469451,
          5.352753879415308,
          7.174386591939361,
          11.264606827298268,
          21.55016485251362,
          None, # 157.58532507782218
          None, # 380.1302925876542
          None] # 699.7387982678907
        ETh_ff_sim_l= [
          0.5907714676944101,
          0.6710283396413771,
          0.8068957630374675,
          1.000439845878967,
          1.3513351803595934,
          2.0536752009713415,
          2.2300187670792453,
          2.429466341104507,
          2.6903484553485266,
          3.000872233306323,
          3.4410421781408966,
          3.935681394156022,
          None, # 4.643556010602577
          None, # 5.59704091169721
          None] # 7.344123932313496
        ETc_ff_sim_l= [
          2.009345254639259,
          1.9861399350607394,
          1.9959564824533977,
          1.9931726443065612,
          1.999767945808216,
          1.9806110151269938,
          2.0030805971338523,
          2.007778635142899,
          1.9774619376748654,
          1.9949836103359992,
          1.993188037134782,
          1.9826785520434822,
          1.99915258644733,
          2.001614284033611,
          2.0032270043889446]
    elif t == 77:
      pass
    else: sim = True
  else: sim = True
  
  mew, ms = 2, 8
  ar_l = []
  for har in [*numpy.linspace(0.05, 0.8*har_ub, 5, endpoint=False), *numpy.linspace(0.8*har_ub, har_ub, 10) ]:
    ar_l.append(har)
    if sim:
      ET = test_simplex_reptoall(num_srun, t, car, har, serv, sdist_m)
      # ETc_reptoall_sim_l.append(ET)
      EThc_reptoall_sim_l.append(ET)
      
      (ETc, ETh) = sim_ff(num_srun, t, car, har, serv, sdist_m)
      ETc_ff_sim_l.append(ETc)
      ETh_ff_sim_l.append(ETh)
  log(WARNING, "EThc_reptoall_sim_l= {}".format(pprint.pformat(EThc_reptoall_sim_l) ) )
  log(WARNING, "ETh_ff_sim_l= {}".format(pprint.pformat(ETh_ff_sim_l) ) )
  log(WARNING, "ETc_ff_sim_l= {}".format(pprint.pformat(ETc_ff_sim_l) ) )
  
  gain_l, pain_l = [], []
  for i, ET_reptoall in enumerate(EThc_reptoall_sim_l):
    if ETh_ff_sim_l[i] is None:
      gain_in_hot = None
      pain_in_cold = None
    else:
      gain_in_hot = (ET_reptoall - ETh_ff_sim_l[i])/ETh_ff_sim_l[i] * 100
      pain_in_cold = (ET_reptoall - ETc_ff_sim_l[i])/ETc_ff_sim_l[i] * 100
    gain_l.append(gain_in_hot)
    pain_l.append(pain_in_cold)
  plot.plot(ar_l, gain_l, label="Gain in hot data", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(ar_l, pain_l, label="Pain in cold data", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'Hot data arrival rate $\lambda$ (Request/s)', fontsize=12)
  plot.ylabel(r'Pain or gain in percentage', fontsize=12)
  plot.title('Pain and gain of Rep-to-all over Fairness-First\n' + r'$S \sim Exp(\mu={})$, $t={}$, $\lambda_c={}$'.format(mu, t, car) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_reptoall_over_ff_t_{}_lc_{}.pdf".format(t, car) )
  log(WARNING, "done; t= {}".format(t) )

if __name__ == "__main__":
  plot_fairnessfirst()
  # plot_reptoall_over_ff()
