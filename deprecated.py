import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import numpy, math, collections
from simplex_sim_components import *

def eg_pgen_psink():
  env = simpy.Environment()
  ps = PacketSink(env, debug=True)  # debugging enable for simple output
  pg1 = PacketGenerator(env, _id="pg1", adist=lambda: 1.5, sdist=lambda: expovariate(0.01) )
  pg2 = PacketGenerator(env, _id="pg2", adist=lambda: 2, sdist=lambda: expovariate(0.01) )
  
  pg1.out = ps
  pg2.out = ps
  #
  env.run(until=20)

def eg_overloaded_m1_q():
  env = simpy.Environment()
  
  ps = PacketSink(env, debug=True)
  pg = PacketGenerator(env, _id="pg", adist=lambda: 1.5, sdist=lambda: 100)
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
  serv_rate = 1.0
  pg = PacketGenerator(env, _id="pg",
                       adist=lambda: random.expovariate(arr_rate),
                       sdist=lambda: 1.0)
                      # sdist=lambda: random.expovariate(serv_rate) )
  # m1_q = S1_Q(_id=0, env=env, rate=1.0) # rate (Bps)
  m1_q = S1_Q(_id=0, env=env, serv_dist=lambda: random.expovariate(serv_rate) ) # rate (Bps)
  # qm = QMonitor(env, q=m1_q, dist=lambda: 0.5)
  qm = QMonitor(env, q=m1_q, dist=lambda: 500)
  
  pg.out = m1_q
  # m1_q.out = ps
  #
  # env.run(until=8000)
  env.run(until=50000)
  # env.run(until=500000)
  # print("Last 10 waits: "  + ", ".join(["{:.3f}".format(x) for x in ps.waits[-10:] ] ) )
  print("Last 10 in qm.n_list: {}".format(qm.n_list[-10:] ) )
  # print("Last 10 sink arrival times: " + ", ".join(["{:.3f}".format(x) for x in ps.arrivals[-10:] ] ) )
  # print("average wait = {:.3f}".format(sum(ps.waits)/len(ps.waits) ) )
  print("received: {}, dropped {}, sent {}".format(m1_q.n_recved, m1_q.n_dropped, pg.n_sent) )
  print("loss rate: {}".format(float(m1_q.n_dropped)/m1_q.n_recved) )
  
  print("arr_rate= {}, serv_rate= {}".format(arr_rate, serv_rate) )
  E_N = arr_rate/(serv_rate - arr_rate)
  print("E[N]= {}, sim_E[N]= {:.3f}".format(E_N, float(sum(qm.n_list) )/len(qm.n_list) ) )
  E_T = E_N/arr_rate
  print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, float(sum(m1_q.qt_list) )/len(m1_q.qt_list) ) )
  E_W = E_T - 1/serv_rate
  print("E[W]= {:.3f}, sim_E[W]= {:.3f}".format(E_W, float(sum(m1_q.wt_list) )/len(m1_q.wt_list) ) )
  
  # plot.plot(qm.t_list, qm.n_list, 'r-')
  plot.step(qm.t_list, qm.n_list, 'r-', where='mid', label='mid')
  # plot.axis([0, 6, 0, 20] )
  # plot.show()
  plot.savefig("m_m_1_q.png")

def test_fj():
  def H_k(k):
    return float(sum([1/float(i) for i in range(1, k + 1) ] ) )
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    arr_rate = 0.5
    serv_rate = 1.0
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = 2
    # qid_list = ["m1_q_%s" % i for i in range(1, num_q) ]
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
    fj_q = MDSQ("fj_q", env, num_q, qid_list, qserv_dist_list)
    # mds_1_q = MDSQ("fj_q", env, 1, qid_list, qserv_dist_list)
    # qm = QMonitor(env, q=m1_q, dist=lambda: 500)
    
    pg.out = fj_q
    
    # env.run(until=5)
    # env.run(until=100)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    E_T = None
    ro = arr_rate/serv_rate
    E_T_2 = (12 - ro)/(8*(serv_rate - arr_rate) )
    if num_q == 2:
      E_T = E_T_2
    else:
      x = H_k(num_q)/H_k(2)
      # Approximate Analysis of Fork/Join Synchronization in Parallel Queues
      E_T = (x + (4.0/11.0)*(1-x)*ro)*E_T_2
    print("num_q= {}, arr_rate= {}, serv_rate= {}".format(num_q, arr_rate, serv_rate) )
    st_list = fj_q.join_sink.st_list
    if len(st_list) > 0:
      sim_E_T = float(sum(st_list) )/len(st_list)
      print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
      diff_list.append(abs(E_T - sim_E_T) )
  print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )
  
  print("fj_q.join_sink.qid__num_win_map= {}".format(pprint.pformat(fj_q.join_sink.qid__num_win_map) ) )
  total_num_wins = sum([n for i, n in fj_q.join_sink.qid__num_win_map.items() ] )
  print("pg.n_sent= {}, total_num_wins= {}".format(pg.n_sent, total_num_wins) )
  qid__win_freq_map = {i:float(n)/total_num_wins for i, n in fj_q.join_sink.qid__num_win_map.items() }
  print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )

def test_mds_n_1(arr_rate, serv_rate, num_q):
  print("test_mds_n_1:: arr_rate= {}, serv_rate= {}, num_q= {}"
        .format(arr_rate, serv_rate, num_q) )
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
    mdsq = MDSQ("mdsq", env, 1, qid_list, qserv_dist_list)
    # qm = QMonitor(env, q=m1_q, dist=lambda: 500)
    
    pg.out = mdsq
    
    # env.run(until=5)
    # env.run(until=100)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    print("num_q= {}, arr_rate= {}, serv_rate= {}".format(num_q, arr_rate, serv_rate) )
    # ro = arr_rate/serv_rate
    E_T = 1.0/(num_q*serv_rate - arr_rate)
    st_list = mdsq.join_sink.st_list
    if len(st_list) > 0:
      sim_E_T = float(sum(st_list) )/len(st_list)
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
  pg = PacketGenerator(env, _id="p_gen",
                       adist=lambda: random.expovariate(arr_rate),
                       sdist=lambda: 1)
  a_q = AVQ("a_q", env, k, r, t, qid_list, qserv_rate_list=qmu_list)
  aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 0.1)
  a_q.join_q.out_m = aq_monitor
  pg.out = a_q
  # env.run(until=50000) # env.run(until=500)
  env.run(until=50000)
  
  total_num_wins = sum([n for i, n in a_q.join_sink.qid__num_win_map.items() ] )
  qid__win_freq_map = {i:float(n)/total_num_wins for i, n in a_q.join_sink.qid__num_win_map.items() }
  print("arr_rate= {}, qid__win_freq_map= {}".format(arr_rate, pprint.pformat(qid__win_freq_map) ) )
  
  st_list = a_q.join_sink.st_list
  if len(st_list) > 0:
    print("sim_E_T= {}".format(float(sum(st_list) )/len(st_list)) )
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

def plot_simplex_w_varying_serv_rate_alloc(num_q):
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
  plot.savefig("plot_simplex_w_varying_serv_rate_alloc__t_{}.png".format(t) )
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

if __name__ == "__main__":
  # random.seed(33)
  # test_simplex()
  
  # plot_dist()
  plot_dist(dist=dolly_slowdown_dist)
  
  # simplex_t_1__zero_state()
  # sum_of_harmonics()
  