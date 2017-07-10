import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import numpy, math, collections
from simplex_sim import *

def eg_pgen_psink():
  env = simpy.Environment()
  ps = PacketSink(env, debug=True)  # debugging enable for simple output
  pg1 = PG(env, _id="pg1", adist=lambda: 1.5, sdist=lambda: expovariate(0.01) )
  pg2 = PG(env, _id="pg2", adist=lambda: 2, sdist=lambda: expovariate(0.01) )
  
  pg1.out = ps
  pg2.out = ps
  #
  env.run(until=20)

def eg_overloaded_m1_q():
  env = simpy.Environment()
  
  ps = PacketSink(env, debug=True)
  pg = PG(env, _id="pg", adist=lambda: 1.5, sdist=lambda: 100)
  m1_q = S1_Q(_id=0, env=env, rate=200.0, qlimit_B=300, debug=True)
  
  pg.out = m1_q
  m1_q.out = ps
  #
  env.run(until=20)
  print("waits: {}".format(ps.waits) )
  print("pg.n_sent= {}\nps.n_recved= {}\nm1_q.n_dropped= {}"
        .format(pg.n_sent, ps.n_recved, m1_q.n_dropped) )

def test_m_m_1():
  env = simpy.Environment()
  
  # ps = PacketSink(env, debug=False, rec_arrivals=True)
  arr_rate = 0.5
  mu = 1.0
  # pg = PG(env, "pg", arr_rate)
  pg = MT_PG(env, "pg", arr_rate, ['a', 'b'], pop_sym='a', pop_arr_rate=0.25)
  # pg = MT_PG(env, "pg", arr_rate, ['a', 'b'], pop_sym=None)
  q = FCFS(_id=0, env=env, rate=mu)
  pg.out = q
  
  qm = QMonitor(env, q=q, dist=lambda: 1)
  pg.init()
  env.run(until=50000)
  # 
  print("received= {}, dropped= {}, sent= {}".format(q.n_recved, q.n_dropped, pg.n_sent) )
  print("arr_rate= {}, mu= {}".format(arr_rate, mu) )
  E_N = arr_rate/(mu - arr_rate)
  print("E_N= {}, E_N_sim= {:.3f}".format(E_N, float(sum(qm.n_l) )/len(qm.n_l) ) )
  E_T = E_N/arr_rate
  print("E_T= {:.3f}, E_T_sim= {:.3f}".format(E_T, float(sum(q.qt_l) )/len(q.qt_l) ) )
  E_W = E_T - 1/mu
  print("E_W= {:.3f}, E_W_sim= {:.3f}".format(E_W, float(sum(q.wt_l) )/len(q.wt_l) ) )
  
  # plot.plot(qm.t_l, qm.n_l, 'r-')
  plot.step(qm.t_l, qm.n_l, 'r-', where='mid', label='mid')
  # plot.axis([0, 6, 0, 20] )
  plot.savefig("m_m_1_q.png")

def test_fj():
  def H_k(k):
    return float(sum([1/float(i) for i in range(1, k + 1) ] ) )
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    arr_rate = 0.5
    mu = 1.0
    pg = PG(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = 2
    # qid_list = ["m1_q_%s" % i for i in range(1, num_q) ]
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_l = [lambda: random.expovariate(mu) for i in range(num_q) ]
    fj_q = MDSQ("fj_q", env, num_q, qid_list, qserv_dist_l)
    # mds_1_q = MDSQ("fj_q", env, 1, qid_list, qserv_dist_l)
    # qm = QMonitor(env, q=m1_q, dist=lambda: 500)
    
    pg.out = fj_q
    
    # env.run(until=5)
    # env.run(until=100)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    E_T = None
    ro = arr_rate/mu
    E_T_2 = (12 - ro)/(8*(mu - arr_rate) )
    if num_q == 2:
      E_T = E_T_2
    else:
      x = H_k(num_q)/H_k(2)
      # Approximate Analysis of Fork/Join Synchronization in Parallel Queues
      E_T = (x + (4.0/11.0)*(1-x)*ro)*E_T_2
    print("num_q= {}, arr_rate= {}, mu= {}".format(num_q, arr_rate, mu) )
    st_l = fj_q.join_sink.st_l
    if len(st_l) > 0:
      sim_E_T = float(sum(st_l) )/len(st_l)
      print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
      diff_list.append(abs(E_T - sim_E_T) )
  print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )
  
  print("fj_q.join_sink.qid__num_win_map= {}".format(pprint.pformat(fj_q.join_sink.qid__num_win_map) ) )
  total_num_wins = sum([n for i, n in fj_q.join_sink.qid__num_win_map.items() ] )
  print("pg.n_sent= {}, total_num_wins= {}".format(pg.n_sent, total_num_wins) )
  qid__win_freq_map = {i:float(n)/total_num_wins for i, n in fj_q.join_sink.qid__num_win_map.items() }
  print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )

def test_mds_n_1(arr_rate, mu, num_q):
  print("test_mds_n_1:: arr_rate= {}, mu= {}, num_q= {}"
        .format(arr_rate, mu, num_q) )
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    pg = PG(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_l = [lambda: random.expovariate(mu) for i in range(num_q) ]
    mdsq = MDSQ("mdsq", env, 1, qid_list, qserv_dist_l)
    # qm = QMonitor(env, q=m1_q, dist=lambda: 500)
    
    pg.out = mdsq
    
    # env.run(until=5)
    # env.run(until=100)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    print("num_q= {}, arr_rate= {}, mu= {}".format(num_q, arr_rate, mu) )
    # ro = arr_rate/mu
    E_T = 1.0/(num_q*mu - arr_rate)
    st_l = mdsq.join_sink.st_l
    if len(st_l) > 0:
      sim_E_T = float(sum(st_l) )/len(st_l)
      # print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
      diff_list.append(abs(E_T - sim_E_T) )
      return sim_E_T
  print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )

def plot_split_merge_approx_error(arr_rate, mu, n):
  gamma = arr_rate
  k__i_list_map, k__p_i_list_map = {}, {}
  k_list, k__avg_active_num_server_list = [], []
  for k in range(2, n+1):
    if k not in k__p_i_list_map:
      k__i_list_map[k], k__p_i_list_map[k] = [], []
    
    p_k_1 = 1/((n-k+1)*mu/gamma * (harmonic_sum(n) - harmonic_sum(n-k+1) ) + 1)
    def p_i(i):
      return (n-k+1)*mu/((n-i)*gamma) * p_k_1 if i < k-1 else p_k_1
    for i in range(0, k):
      k__i_list_map[k].append(i)
      k__p_i_list_map[k].append(p_i(i) )
    #
    k_list.append(k)
    k__avg_active_num_server_list.append((n-k+1)*(k-1)*mu/gamma * p_k_1)
  # print("k__i_list_map=\n{}".format(pprint.pformat(k__i_list_map) ) )
  # print("k__p_i_list_map=\n{}".format(pprint.pformat(k__p_i_list_map) ) )
  plot.figure(1)
  plot.subplot(211)
  for k in k__i_list_map:
    plot.plot(k__i_list_map[k], k__p_i_list_map[k], '-o', label="k={}".format(k), color=plot_color_list[k%len(plot_color_list) ] )
    # plot.legend()
  plot.xlabel('i')
  plot.ylabel(r'$p_i$')
  plot.title(r'n= {}, $\lambda$= {}, $\mu$= {}'.format(n, arr_rate, mu) )
  plot.subplot(212)
  plot.plot(k_list, k__avg_active_num_server_list, 'b-o')
  # plot.legend()
  plot.xlabel('k')
  plot.ylabel('Average active number of servers')
  # plot.title(r'n= {}, $\lambda$= {}, $\mu$= {}'.format(n, arr_rate, mu) )
  plot.tight_layout()
  # ax = plot.gca()
  # ax.grid(True)
  plot.savefig("plot_split_merge_approx_error__n_{}.png".format(n) )

def test_simplex(arr_rate=None):
  k, r, t = 2, 2, 1
  if arr_rate == None:
    arr_rate = 1.1
  num_q = int(1 + t*r)
  qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
  qmu_list = [0.5, 1, 1]
  # qmu_list = [1, 1, 1]
  log(WARNING, "arr_rate= {}, k= {}, r= {}, t= {}, qmu_list= {}".format(arr_rate, k, r, t, pprint.pformat(qmu_list) ) )
  
  env = simpy.Environment()
  pg = PG(env, _id="p_gen",
                       adist=lambda: random.expovariate(arr_rate),
                       sdist=lambda: 1)
  a_q = AVQ("a_q", env, k, r, t, qid_list, qmu_list=qmu_list)
  aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 0.1)
  a_q.join_q.out_m = aq_monitor
  pg.out = a_q
  # env.run(until=50000) # env.run(until=500)
  env.run(until=50000)
  
  total_num_wins = sum([n for i, n in a_q.join_sink.qid__num_win_map.items() ] )
  qid__win_freq_map = {i:float(n)/total_num_wins for i, n in a_q.join_sink.qid__num_win_map.items() }
  print("arr_rate= {}, qid__win_freq_map= {}".format(arr_rate, pprint.pformat(qid__win_freq_map) ) )
  
  st_l = a_q.join_sink.st_l
  if len(st_l) > 0:
    print("sim_E_T= {}".format(float(sum(st_l) )/len(st_l)) )
  # """
  print("\n")
  # print("aq_monitor.polled_state__counter_map= {}".format(pprint.pformat(aq_monitor.polled_state__counter_map) ) )
  total_counter = sum([c for rs, c in aq_monitor.polled_state__counter_map.items() ] )
  polled_state__counter_map = {rs:float(c)/total_counter for rs, c in aq_monitor.polled_state__counter_map.items() }
  print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
  return polled_state__counter_map['0,(0,0)']
  
  print("\n")
  # print("aq_monitor.state__num_found_by_job_departed_map= {}".format(pprint.pformat(aq_monitor.state__num_found_by_job_departed_map) ) )
  total_counter = sum([c for rs, c in aq_monitor.state__num_found_by_job_departed_map.items() ] )
  state__freq_found_by_job_departed_map = {rs:float(c)/total_counter for rs, c in aq_monitor.state__num_found_by_job_departed_map.items() }
  print("state__freq_found_by_job_departed_map= {}".format(pprint.pformat(state__freq_found_by_job_departed_map) ) )
  
  print("\n")
  # print("aq_monitor.start_setup__num_found_by_job_departed_map= {}".format(pprint.pformat(aq_monitor.start_setup__num_found_by_job_departed_map) ) )
  total_counter = sum([c for rs, c in aq_monitor.start_setup__num_found_by_job_departed_map.items() ] )
  start_setup__freq_found_by_job_departed_map = {rs:float(c)/total_counter for rs, c in aq_monitor.start_setup__num_found_by_job_departed_map.items() }
  print("start_setup__freq_found_by_job_departed_map= {}".format(pprint.pformat(start_setup__freq_found_by_job_departed_map) ) )
  # """

def plot_simplex_w_varying_mu_alloc(num_q):
  k, r, t = 2, 2, 1 # 2
  mu = 1.0
  arr_rate_ub = simplex_inner_bound_on_arr_rate(r, t, mu)
  arr_rate = mu
  log(WARNING, "k= {}, r= {}, t= {}, mu= {}, arr_rate_ub={}".format(k, r, t, mu, arr_rate_ub) )
  
  Cap = num_q*mu
  def qmu_l(c):
    mu_ = Cap/(c+2)
    gamma = c*mu_
    return [gamma, mu_, mu_]
  color = iter(cm.rainbow(numpy.linspace(0, 1, 10) ) )
  # for c in numpy.arange(0.05, 1, 0.1):
  for arr_rate in numpy.linspace(0.2, mu, 3):
    c_l = []
    E_T_simplex_sm_l, sim_hetero_E_T_simplex_l = [], []
    hetero_E_T_simplex_l = []
    for c in numpy.linspace(0.25, 5, 10):
      c_l.append(c)
      # sim
      E_T_simplex_sm_l.append(E_T_simplex_sn(t, arr_rate, mu, c) )
      
      num_f_run = 3
      sim_hetero_simplex_E_T = test_avq(num_f_run, arr_rate, mu, k, r, t, qmu_l(c) )
      sim_hetero_E_T_simplex_l.append(sim_hetero_simplex_E_T)
      
      hetero_E_T_simplex_l.append(simplex_w_one_repair__E_T(arr_rate, mu, c) )
    plot.plot(c_l, E_T_simplex_sm_l, 'o', color=next(color), label="SM-Simplex $\lambda={0:0.2f}$".format(arr_rate) )
    plot.plot(c_l, sim_hetero_E_T_simplex_l, 'o', color=next(color), label=r'Simplex $\lambda={0:0.2f}$'.format(arr_rate) )
    plot.plot(c_l, hetero_E_T_simplex_l, 'o', color=next(color), label=r'LB-Simplex $\lambda={0:0.2f}$'.format(arr_rate) )
    
  plot.xlabel("c")
  plot.ylabel("E[T] (s)")
  plot.title(r'Simplex(t:{}), heterogeneous servers; $c=\gamma/\mu$'.format(t) )
  plot.savefig("plot_simplex_w_varying_mu_alloc__t_{}.png".format(t) )
  log(WARNING, "done; n= {} k= {} t= {}".format(num_q, k, t) )

def simplex_t_1__zero_state():
  arr_rate__zero_state_prob_map = {}
  for arr_rate in numpy.arange(0.1, 1.3, 0.1):
    arr_rate__zero_state_prob_map[arr_rate] = test_simplex(arr_rate)
  log(WARNING, "arr_rate__zero_state_prob_map= {}".format(pprint.pformat(arr_rate__zero_state_prob_map) ) )

def plot_dist(dist=None):
  N = 10000 * 100
  if dist == None:
    data = [random.expovariate(1) for i in range(N) ]
  else:
    data = [dist() for i in range(N) ]
  plot.clf()
  
  # plot.hist(data, bins='auto', normed=True)
  
  # w_l = numpy.ones_like(data)/float(len(data) )
  # plot.hist(data, weights=w_l)
  
  plot.hist(data, bins=numpy.arange(min(data)-0.5, max(data)-0.5, 1), normed=True)
  d_f_m = {d:f/len(data) for d,f in collections.Counter(data).items() }
  print("d_f_m= {}".format(d_f_m) )
  
  plot.savefig("plot_dist.png")

def plot_binomial_dist__approx():
  n = 15
  p = 0.4
  def comp_dist(k):
    if k > n:
      return 0
    sum_ = 0
    for i in range(k, n+1):
      sum_ += binomial(n, i) * p**i * (1-p)**(n-i)
    return sum_
  def dist(k):
    if k > n:
      return 0
    sum_ = 0
    for i in range(0, k+1):
      sum_ += binomial(n, i) * p**i * (1-p)**(n-i)
    return sum_
  def chernoff_bound_on_upper_tail(k):
    p_ = k/n
    print("p_= {}".format(p_) )
    return math.exp(n*((p_*math.log(p/p_) ) + (1-p_)*math.log((1-p)/(1-p_) ) ) )
  
  k_l, dist_l, approx_dist_l = [], [], []
  # for k in range(0, n+2):
  for k in range(int(p*n), n):
    k_l.append(k)
    dist_l.append(comp_dist(k) )
    approx_dist_l.append(chernoff_bound_on_upper_tail(k) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  # print("k_l= {}".format(pprint.pformat(k_l) ) )
  plot.plot(k_l, approx_dist_l, color='red', label=r'$Pr\{N \leq k\}_{UB}$', marker=next(marker), linestyle='', mew=2)
  plot.plot(k_l, dist_l, color='black', label=r'$Pr\{N \leq k\}$', marker=next(marker), linestyle='', mew=2)
  plot.xlabel(r'$k$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, p= {}'.format(n, p) )
  plot.savefig("plot_binomial_dist__approx_n_{}.png".format(n) )
  plot.gcf().clear()
  log(WARNING, "done; n= {}, p= {}".format(n, p) )

def sum_of_harmonics():
  k = 10
  mu = 1
  
  def E_C(n):
    sum_ = sum([H(n-i) for i in range(1, k+1) ] )
    return 1/mu*(k*H(n) - sum_ + (n-k)*(H(n)-H(n-k) ) )
  
  for n in range(k, 4*k):
    print("n= {}, E_C= {}".format(n, E_C(n) ) )

# **********************************  Fairness First MDSQ  ******************************* #
class MT_MDS_JQ(object):
  def __init__(self, _id, env, k, input_qid_l, sym_sysqid_map):
    self._id = _id
    self.env = env
    self.input_qid_l = input_qid_l
    self.sym_sysqid_map = sym_sysqid_map
    
    self.job_id__p_l_map = {}
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.out = None
    self.out_c = None
    self.out_m = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
  
  def __repr__(self):
    return "MT_MDS_JQ[_id= {}, input_qid_l= {}]".format(self._id, self.input_qid_l)
  
  def check_for_job_completion(self, p):
    p_l = self.job_id__p_l_map[p.job_id]
    recved_from_qid_l = [p.prev_hop_id for p in p_l]
    sys_qid = self.sym_sysqid_map[p.sym]
    success = False
    
    if sys_qid in p_l or len(p_l) >= k:
      self.out_c.put_c(CPacket(_id=p.job_id, sym=p.sym, prev_hop_id=self._id, departed_qid_l=recved_from_qid_l) )
      
      p.winner_id = p.prev_hop_id
      p.prev_hop_id = self._id
      self.out.put(p)
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      if p.job_id not in self.job_id__p_l_map:
        self.job_id__p_l_map[p.job_id] = []
      self.job_id__p_l_map[p.job_id].append(p)
      self.check_for_job_completion(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.ref_time = self.env.now
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      self.job_id__p_l_map.pop(p.job_id, None)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)

class FF_MDSQ(object): # Fairness First
  def __init__(self, _id, env, n, k, sym_l, serv, out=None):
    self._id = _id
    self.env = env
    self.n = n
    self.k = k
    self.sym_l = sym_l
    self.out = out
    
    self.qid_l = [i for i in range(self.num_q) ]
    self.sym_sysqid_map = {s:i for i,s in enumerate(self.sym_l) }
    
    self.jsink = FF_JSink(_id, env)
    self.jsink.out = out
    self.id_q_map = {}
    
    self.jq = MT_MDS_JQ(_id, env, self.qid_l, sym_sysqid_map)
    self.jq.out = self.jsink
    self.jq.out_c = self
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
  
  def __repr__(self):
    return "MT_MDSQ[n= {}, k= {}]".format(self.n, self.k)
  
  def state(self):
    return {i:q.length() for i,q in self.id_q_map.items() }
  
  def send_pop(self):
    while True:
      p = (yield self.pop_store.get() )
      
      sys_qid = self.sym_sysqid_map[p.sym]
      sys_q = self.id_q_map[sys_qid]
      if sys_q.busy:
        self.release_pop = self.env.event()
        yield self.release_pop
      
      serving_qid_l = []
      sys_q.put(p.deep_copy() )
      serving_qid_l.append(sys_qid)
      
      qid_l_ = []
      for qid in list(self.qid_l).remove(sys_qid):
        if not (q.busy and q.p_in_serv.sym != POP_SYM): # q may be busy with pop_sym because simpy did not register cancellation yet
          qid_l_.append(qid)
      
      if len(qid_l_) >= k:
        for qid in qid_l_:
          self.id_q_map[qid].put(p.deep_copy() )
          serving_qid_l.append(qid)
      self.starttype__num_map[len(serving_qid_l)-1] += 1
      self.job_id__serving_qid_l_map[p.job_id] = serving_qid_l
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      
      if p.sym == POP_SYM:
        self.pop_store.put(p)
      else:
        sys_q = self.id_q_map[self.sym_sysqid_map[p.sym] ]
        if sys_q.p_in_serv is not None and sys_q.p_in_serv.sym == POP_SYM:
          ji = sys_q.p_in_serv.job_id
          serving_qid_l = self.job_id__serving_qid_l_map[ji]
          serving_qid_l.remove(self.sym_sysqid_map[p.sym] )
          
          sys_q.put_c(CPacket(_id=ji, prev_hop_id=self._id) )
          if len(serving_qid_l) < k:
            for qid in serving_qid_l:
              self.id_q_map[r].put_c(CPacket(_id=ji, prev_hop_id=self._id) )
        sys_q.put(p.deep_copy() )
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.job_id_counter += 1
    p.job_id = self.job_id_counter
    return self.store.put(p)
  
  # def run_c(self):
  #   while True:
  #     cp = (yield self.store_c.get() )
  #     # 
  #     # log(WARNING, "recved cp= {}".format(cp) )
  #     # print("state= {}".format(pprint.pformat(self.state() ) ) )
      
  #     for g, q in self.id_q_map.items():
  #       if q._id not in cp.departed_qid_l:
  #         q.put_c(cp.deep_copy() )
      
  #     if cp.sym == POP_SYM:
  #       if self.release_pop is not None:
  #         self.release_pop.succeed()
  #         self.release_pop = None
  
  # def put_c(self, cp):
  #   sim_log(DEBUG, self.env, self, "recved", cp)
  #   return self.store_c.put(cp)

if __name__ == "__main__":
  test_m_m_1()
  # test_simplex()
  
  # plot_dist()
  # plot_dist(dist=dolly_slowdown_dist)
  
  # simplex_t_1__zero_state()
  # sum_of_harmonics()
  