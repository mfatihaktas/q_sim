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
            # q = self.id_q_map[r]
            # sim_log(WARNING, self.env, self, "busy q!", p)
            # print("q= {}\n q.p_inserv= {}, q.p_l= {}".format(q, q.p_inserv, q.p_l) )
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
  serv = "Pareto" # "Exp" # "Pareto" # "Dolly"
  car = 0.3 # 0.1 # 0.5
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
  V = rv_from_m(sdist_m)
  log(WARNING, "t= {}, car= {}, sdist_m= {}, har_ub= {}".format(t, car, sdist_m, har_ub) )
  
  ETc_sim_l, ETh_sim_l = [], []
  ETh_ub_l, ETh_lb_l, ETh_approx_l = [], [], []
  
  num_srun = 3
  sim_simplex = False
  if serv == "Exp":
    if t == 1:
      if car == 0.1:
        ETh_sim_l= [
          0.7409446245913069,
          0.8691944553872436,
          1.0596513900016429,
          1.3430765997093623,
          1.9143424134519744,
          3.08430647071154,
          3.4120498780040727,
          3.8258463740663924,
          4.283456962750768,
          4.989571474137262,
          6.152026410860675,
          7.440060692819363,
          9.505361521577766,
          None, # 15.067634355415697,
          None] # 22.767168502143562]
      elif car == 0.5:
        ETh_sim_l= [
          0.9213906136285814,
          1.079391764141003,
          1.316370892694153,
          1.6957597966631406,
          2.3821791307448024,
          3.98475560438631,
          4.342188135067258,
          4.824047459304102,
          5.590832654385245,
          6.258986176974797,
          7.381370440746494,
          8.733970147990197,
          11.30012985193308,
          None, # 13.067793614094803,
          None] # 26.007985252827424]
    elif t == 3:
      if car == 0.1:
        ETh_sim_l= [
          0.4884404089708351,
          0.5599196256213207,
          0.6680821792758612,
          0.8409653112226637,
          1.1496162745785325,
          1.8904759704555996,
          2.1399639333981635,
          2.3053084133910713,
          2.695316742054405,
          3.091031654333508,
          3.8030094936237986,
          4.540371522911487,
          5.882876885214188,
          8.876245687732691,
          18.939562973061]
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
  elif serv == "Pareto":
    if t == 1:
      if car == 0.1:
        # loc, a = 1, 3
        ETh_sim_l= [
          1.3846238320956463,
          1.5462852274928522,
          1.767452898642065,
          2.1366627984341435,
          2.7224216797464003,
          4.307322557250082,
          4.653124052893923,
          5.565779098536901,
          6.010815631470543,
          6.795093510347267,
          8.834250495532023,
          11.070140430188268,
          16.40864497329819,
          None, # 34.009466918567675,
          None] # 86.0404055347701]
      elif car == 0.3:
        ETh_sim_l= [
          1.462211187084237,
          1.6486059593502425,
          1.9146010529140192,
          2.388782884574221,
          3.041516781416121,
          4.888687307178,
          5.641152074380211,
          6.053718044645767,
          6.992545442656206,
          8.117053085690948,
          10.082084657522154,
          13.245017062284292,
          None, # 18.62550055240398,
          None, # 35.50299305129379,
          None] # 99.77979372603328]
      elif car == 0.5:
        ETh_sim_l= [
          1.551260625733348,
          1.7272026851183273,
          2.016358819849658,
          2.4707991007399293,
          3.4357747815241146,
          5.144106330516206,
          5.83749870316689,
          6.786917446173142,
          7.963602903118772,
          8.686250854689352,
          10.323510959844379,
          14.193429330450575,
          None, # 21.435936778092756,
          None, # 46.936136975195645,
          None] # 91.4627632437775]
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
  for har in [*numpy.linspace(0.05, 0.8*har_ub, 5, endpoint=False), *numpy.linspace(0.8*har_ub, har_ub, 10) ]:
  # for har in numpy.linspace(0.5, 0.5, 1):
    ar_l.append(har)
    if sim_simplex:
      (ETc, ETh) = sim_ff(num_srun, t, car, har, sdist_m)
      ETc_sim_l.append(ETc)
      ETh_sim_l.append(ETh)
    # ETh_approx = ff_ETh_approx_(har, car, t, sdist_m)
    # print("ETh_approx= {}".format(ETh_approx) )
    ETh_approx = ff_ETh_newapprox(har, car, t, sdist_m)
    print("ETh_newapprox= {}".format(ETh_approx) )
    ETh_approx_l.append(ETh_approx)
  log(WARNING, "ETh_sim_l= {}".format(pprint.pformat(ETh_sim_l) ) )
  plot.plot(ar_l, ETh_sim_l, label="Simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(ar_l, ETh_approx_l, label="Approximation", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  ar_l = []
  har_ub_min = ff_har_ub_min(t, sdist_m)
  # for har in [*numpy.linspace(0.05, 0.8*har_ub_min, 5, endpoint=False), *numpy.linspace(0.8*har_ub_min, har_ub_min, 10) ]:
  for har in numpy.linspace(0.05, har_ub_min, 100):
    ar_l.append(har)
    ETh_ub_l.append(ff_ETh_ub(har, t, sdist_m) )
  plot.plot(ar_l, ETh_ub_l, label="Upper bound", color=next(dark_color), linestyle='--', lw=2)
  
  ar_l = []
  har_ub_max = ff_har_ub_max(t, sdist_m)
  # for har in [*numpy.linspace(0.05, 0.8*har_ub_max, 5, endpoint=False), *numpy.linspace(0.8*har_ub_max, har_ub_max, 10) ]:
  for har in numpy.linspace(0.05, har_ub_max, 100):
    ar_l.append(har)
    ETh_lb_l.append(ff_ETh_lb(har, t, sdist_m) )
  plot.plot(ar_l, ETh_lb_l, label="Lower bound", color=next(dark_color), linestyle='-.', lw=2)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'Hot data arrival rate $\lambda$', fontsize=14)
  plot.ylabel(r'Avg hot data download time', fontsize=14)
  # plot.title(r'Fairness-first $t={}$, Servers $\sim {}$, $\lambda_c={}$'.format(t, V, car) )
  plot.title(r'Fairness-first $t= {}$, Servers$\sim {}$'.format(t, V) + "\n" + r'Cold data arrival rate= {}'.format(car) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.3)
  fig.tight_layout()
  plot.savefig("plot_ff_t{}_{}_car_{}.pdf".format(t, serv, car) )
  log(WARNING, "done; t= {}".format(t) )

# *******************************************  Rep-to-all  *************************************** #
def sim_reptoall(num_srun, t, car, har, sdist_m):
  ET_sum = 0
  for f in range(num_srun):
    log(WARNING, "t= {}, car= {}, har= {}, sdist_m= {}, ".format(t, car, har, sdist_m) )
    
    (sym_l, sym__rgroup_l_m) = simplex_sym_l__sym__rgroup_l_m(t)
    log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
    
    env = simpy.Environment()
    pg = MT_PG(env, "pg", (len(sym_l)-1)*car, sym_l, HOT_SYM, har)
    avq = MT_AVQ("mt_avq", env, t, sym__rgroup_l_m, sdist_m)
    pg.out = avq
    pg.init()
    env.run(until=2*50000)
    
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
  return ET

def plot_reptoall_vs_ff():
  t = 1
  serv = "Exp" # "Dolly"
  car = 0.5 # 0.1 # 0.5
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
  V = rv_from_m(sdist_m)
  har_ub *= 0.8
  log(WARNING, "t= {}, car= {}, sdist_m= {}, har_ub= {}".format(t, car, sdist_m, har_ub) )
  
  ETc_reptoall_l, EThc_reptoall_l = [], []
  ETc_ff_l, ETh_ff_l = [], []
  
  num_srun = 5
  sim = False # True
  if serv == "Exp":
    if t == 1:
      if car == 0.1:
        EThc_reptoall_l= [
          0.7177278908809522,
          0.787854715660468,
          0.8843059766403686,
          1.01115883628257,
          1.1645103388284417,
          1.4153819150532145,
          1.4485019610292265,
          1.5051838624607468,
          1.5422094438301317,
          1.6076951807154118,
          1.6606679793218717,
          1.7352947059527633,
          1.8289551708906704,
          1.8888738930821816,
          1.9708155101648706]
        ETh_ff_l= [
          0.7290809989976729,
          0.8397711006160605,
          0.963281224987861,
          1.1413494209470774,
          1.4069706607282797,
          1.837523214196881,
          1.906570583309542,
          2.029332360218048,
          2.1045127158381036,
          2.2501605597167234,
          2.358511980918369,
          2.5965759411384863,
          2.6646908485377234,
          2.891551534492633,
          3.1821415153859913]
        ETc_ff_l= [
          1.111910834380286,
          1.1073307992556616,
          1.137540659437823,
          1.1132737587908357,
          1.092300419429767,
          1.1064880608649308,
          1.0959154457973779,
          1.1142249152952592,
          1.1091535187111279,
          1.1103199379895863,
          1.10947072241001,
          1.1126121300773992,
          1.1114864059929912,
          1.1093576986128773,
          1.1307572696780974]
      elif car == 0.5:
        # EThc_reptoall_l= [
        #   0.9357724942120637,
        #   1.0195073997857127,
        #   1.1390600110812168,
        #   1.2839997790832902,
        #   1.508433291917468,
        #   1.8353683210888068,
        #   1.9182499166648572,
        #   2.005499474648738,
        #   2.1262070566662343,
        #   2.158185158260862,
        #   2.267768633383888,
        #   2.331300182993372,
        #   2.4693773038135114,
        #   2.558164399620566,
        #   2.6993829524203203]
        EThc_reptoall_l= [
          0.9412838787757932,
          1.0226867661452232,
          1.1271541586439433,
          1.2852656403912046,
          1.5213904086219823,
          1.8484614880081793,
          1.9329176879280179,
          1.9662625487946195,
          2.080318851000592,
          2.1613783654833854,
          2.2397408521746343,
          2.333284621648475,
          2.457727501060236,
          2.590198586725196,
          2.730893670865941]
        # ETh_ff_l= [
        #   0.9235816250143966,
        #   1.0432936019456611,
        #   1.210933971711239,
        #   1.4425091636665603,
        #   1.7765407222538554,
        #   2.3481727775719268,
        #   2.431752034302232,
        #   2.5534547595530372,
        #   2.6230174844081,
        #   2.86747562230979,
        #   3.08040820848498,
        #   3.2263201208846035,
        #   3.3525433256976203,
        #   3.6285838718416534,
        #   3.840151807003022]
        ETh_ff_l= [
          0.8988449705708461,
          1.0422460107787392,
          1.2072340729101128,
          1.4348460941999397,
          1.7785372934364996,
          2.32529339487659,
          2.4361372855169408,
          2.5667379736338223,
          2.7055939669080917,
          2.8449083674648623,
          3.0538179996365695,
          3.1869798429664042,
          3.35925912433934,
          3.6086801263382773,
          3.9278158322372647]
        # ETc_ff_l= [
        #   2.002169040112291,
        #   1.9865201856659145,
        #   2.0029129814973348,
        #   2.0026859424489474,
        #   1.9921846936850471,
        #   2.00073164130422,
        #   1.9992582177117828,
        #   2.027578688717941,
        #   1.984036470630879,
        #   1.997604087674219,
        #   1.9904973362315517,
        #   1.9850721470152282,
        #   2.0006263019829813,
        #   1.9936189860324502,
        #   2.01215295024685]
        ETc_ff_l= [
          1.9799401626936597,
          2.0128526627309293,
          1.993172187101203,
          2.0059887312972533,
          1.9950806941435169,
          2.002054238086192,
          2.0116933638566,
          2.0080279888803956,
          1.9877198454729186,
          2.011780919800691,
          2.0153969494113317,
          2.011096360545667,
          1.9827236077353896,
          2.016087449883659,
          1.997929666975192]
    elif t == 33:
      if car == 0.1:
        pass
      elif car == 0.5:
        pass
    else: sim = True
  
  mew, ms = 2, 6
  ar_l = []
  for har in [*numpy.linspace(0.05, 0.8*har_ub, 5, endpoint=False), *numpy.linspace(0.8*har_ub, har_ub, 10) ]:
    ar_l.append(har)
    if sim:
      ET = sim_reptoall(num_srun, t, car, har, sdist_m)
      # ETc_reptoall_l.append(ET)
      (ETc, ETh) = sim_ff(num_srun, t, car, har, sdist_m)
      
      EThc_reptoall_l.append(ET)
      ETc_ff_l.append(ETc)
      ETh_ff_l.append(ETh)
  
  log(WARNING, "EThc_reptoall_l= {}".format(pprint.pformat(EThc_reptoall_l) ) )
  log(WARNING, "ETh_ff_l= {}".format(pprint.pformat(ETh_ff_l) ) )
  log(WARNING, "ETc_ff_l= {}".format(pprint.pformat(ETc_ff_l) ) )
  
  gain_in_hot_l, gain_in_cold_l = [], []
  for i, ET_reptoall in enumerate(EThc_reptoall_l):
    if ETh_ff_l[i] is None:
      gain_in_hot = None
      gain_in_cold = None
    else:
      gain_in_hot = (ETh_ff_l[i] - ET_reptoall)/ETh_ff_l[i] * 100
      gain_in_cold = (ETc_ff_l[i] - ET_reptoall)/ETc_ff_l[i] * 100
    gain_in_hot_l.append(gain_in_hot)
    gain_in_cold_l.append(gain_in_cold)
  plot.plot(ar_l, gain_in_hot_l, label="In hot data", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  plot.plot(ar_l, gain_in_cold_l, label="In cold data", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'Hot data arrival rate', fontsize=14)
  plot.ylabel(r'Gain in percentage', fontsize=14)
  # plot.title('Gain of Rep-to-all over Fairness-First\n' + r'Servers $\sim {}$, $t={}$, $\lambda_c={}$'.format(V, t, car) )
  plot.title(r'Servers$\sim {}$, $t={}$'.format(V, t) + '\nCold data arrival rate= {}'.format(car) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_reptoall_vs_ff_{}_t{}_car_{}.pdf".format(serv, t, car) )
  log(WARNING, "done; t= {}".format(t) )

if __name__ == "__main__":
  # plot_fairnessfirst()
  plot_reptoall_vs_ff()
