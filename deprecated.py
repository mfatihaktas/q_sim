import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import numpy, math
from sim_components import *

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

def test_simplex():
  k, r, t = 2, 2, 1
  arr_rate = 0.9
  num_q = int(1 + t*r)
  qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
  # qmu_list = [0.01, 1, 1]
  # qmu_list = [0.001, 1, 1]
  qmu_list = [0.00001, 1.0, 1.0]
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
  """
  print("\n")
  # print("aq_monitor.polled_state__counter_map= {}".format(pprint.pformat(aq_monitor.polled_state__counter_map) ) ) 
  total_counter = sum([c for rs, c in aq_monitor.polled_state__counter_map.items() ] )
  polled_state__counter_map = {rs:float(c)/total_counter for rs, c in aq_monitor.polled_state__counter_map.items() }
  print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
  
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
  """

def plot_exp_dist(mu_arg, func = None):
  # for mu in [2.0, 1.0, 0.00001]:
  #   exp_density = lambda: random.expovariate(mu)
  #   # exp_density = lambda: -math.log(1.0-random.random() )/mu
  #   data = [exp_density() for i in range(10000) ]
  #   plot.hist(data, bins='auto')
  #   plot.title(r'mu= {}'.format(mu) )
  #   plot.savefig("plot_exp_dist_mu_{}.png".format(mu) )
  
  if func == None:
    data = [random.expovariate(mu_arg) for i in range(10000) ]
  else:
    data = [func() for i in range(10000) ]
  plot.clf()
  plot.hist(data, bins='auto')
  plot.title(r'mu= {}'.format(mu_arg) )
  plot.savefig("plot_exp_dist_mu_{}.png".format(mu_arg) )
  
  # def exp_cum_dist(mu, x):
  #   return 1 - math.exp(-mu*x)
  
  # # mu_list = numpy.linspace(1, 100, 10)
  # mu_list = numpy.logspace(-3, 2, 10)
  # x_list = numpy.arange(0.01, 4, 0.01)
  # color = iter(cm.rainbow(numpy.linspace(0, 1, len(mu_list) ) ) )
  # for mu in mu_list:
  #   F_x_list = [exp_cum_dist(mu, x) for x in x_list]
  #   plot.plot(x_list, F_x_list, '.', color=next(color), label=r'$\mu$={}'.format(mu) )
  #   plot.legend()
  # plot.xlabel("x")
  # plot.ylabel("F(x)")
  # plot.savefig("plot_exp_dist_mu_{}.png".format(mu_arg) )

def prob_comp_start():
  mu = 1.0
  simplex_E_S_p = 0.5/mu
  simplex_E_S_c = 0.665/mu
  simplex_E_S = 0.4*simplex_E_S_p + 0.6*simplex_E_S_c
  simplex_E_S_p_2 = 0.5/(mu**2)
  simplex_E_S_c_2 = 7/(9*(mu**2) )
  simplex_E_S_2 = 0.4*simplex_E_S_p_2 + 0.6*simplex_E_S_c_2
  for arr_rate in numpy.arange(0.1, 1.2, 0.05):
    E_X = 1/arr_rate
    """
    a = simplex_E_S_c - simplex_E_S_p
    b = 2*simplex_E_S_p - simplex_E_S_c - E_X
    c = E_X*(E_X-simplex_E_S_p)/(simplex_E_S_c-simplex_E_S_p) - simplex_E_S_p
    print("a= {}, b= {}, c= {}".format(a, b, c) )
    delta = b**2 - 4*a*c
    if delta < 0:
      r_1_real = -b/(2*a)
      r_1_complex = math.sqrt(abs(delta) )/(2*a)
      print("r_1= {} + i{}".format(r_1_real, r_1_complex) )
    else:
      r_1 = (-b + math.sqrt(delta) )/(2*a)
      print("r_1= {}".format(r_1) )
    """
    """
    a = 1
    b = simplex_E_S_c - E_X
    c = E_X**2
    delta = b**2 - 4*a*c
    print("b**2= {}, 4*a*c= {}, delta= {}".format(b**2, 4*a*c, delta) )
    """
    E_S = E_X
    p_c = (E_S - simplex_E_S_p)/(simplex_E_S_c - simplex_E_S_p)
    print("arr_rate= {}, p_c= {}".format(arr_rate, p_c) )

if __name__ == "__main__":
  random.seed(33)
  test_simplex()
  # plot_exp_dist(0)
  # plot_exp_dist(1)
  # plot_exp_dist(0.1)
  # plot_exp_dist(0.01)
  # plot_exp_dist(0.001)
  # plot_exp_dist(0.0001)
  # plot_exp_dist(0.00001)
  # plot_exp_dist(0.000001)
  