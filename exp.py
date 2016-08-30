import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from cycler import cycler
from random import expovariate
import pprint, numpy, simpy

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

def test_mds_n_k(arr_rate, serv_rate, num_q, k):
  print("test_mds_n_k:: arr_rate= {}, serv_rate= {}, num_q= {}, k= {}"
        .format(arr_rate, serv_rate, num_q, k) )
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
    mdsq = MDSQ("mdsq", env, k, qid_list, qserv_dist_list)
    mdsq_monitor = MDSQMonitor(env, q=mdsq, poll_dist=lambda: 1)
    
    pg.out = mdsq
    
    # env.run(until=5)
    # env.run(until=5000)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    print("num_q= {}, arr_rate= {}, serv_rate= {}".format(num_q, arr_rate, serv_rate) )
    # ro = arr_rate/serv_rate
    # E_T = 1.0/(num_q*serv_rate - arr_rate)
    st_list = mdsq.join_sink.st_list
    if len(st_list) > 0:
      sim_E_T = float(sum(st_list) )/len(st_list)
      # print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
      # print("sim_E[T]= {:.3f}".format(sim_E_T) )
      # diff_list.append(abs(E_T - sim_E_T) )
      return sim_E_T
  # print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )
  print("\n")
  print("mdsq.join_sink.qid__num_win_map= {}".format(pprint.pformat(mdsq.join_sink.qid__num_win_map) ) )
  total_num_wins = sum([n for i, n in mdsq.join_sink.qid__num_win_map.items() ] )
  print("pg.n_sent= {}, total_num_wins= {}".format(pg.n_sent, total_num_wins) )
  qid__win_freq_map = {i:float(n)/total_num_wins for i, n in mdsq.join_sink.qid__num_win_map.items() }
  print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  
  print("\n")
  print("mdsq_monitor.state__counter_map= {}".format(pprint.pformat(mdsq_monitor.state__counter_map) ) ) 
  total_counter = sum([c for rs, c in mdsq_monitor.state__counter_map.items() ] )
  polled_state__counter_map = {rs:float(c)/total_counter for rs, c in mdsq_monitor.state__counter_map.items() }
  print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
  
def test_simplex_q(arr_rate, serv_rate, k, r):
  print("test_simplex_q:: arr_rate= {}, serv_rate= {}, k= {}, r= {}"
        .format(arr_rate, serv_rate, k, r) )
  diff_list = []
  env = simpy.Environment()
  pg = PacketGenerator(env, _id="p_gen",
                       adist=lambda: random.expovariate(arr_rate),
                       sdist=lambda: 1)
  t = 1
  num_q = 1 + t*r
  qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
  qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
  
  a_q = AVQ("a_q", env, k, r, t, qid_list, qserv_dist_list)
  aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 0.1)
  a_q.join_q.out_m = aq_monitor
  
  pg.out = a_q
  
  # env.run(until=5)
  # env.run(until=500)
  env.run(until=50000)
  # env.run(until=200000)
  
  # for num_q= 2
  # print("num_q= {}, arr_rate= {}, serv_rate= {}".format(num_q, arr_rate, serv_rate) )
  # print("a_q= {}".format(a_q) )
  # ro = arr_rate/serv_rate
  # E_T_m_m_1 = 1.0/(serv_rate-arr_rate)
  # E_T_f_j_2 = (12-ro)/8*E_T_m_m_1
  st_list = a_q.join_sink.st_list
  E_T = 0
  if len(st_list) > 0:
    sim_E_T = float(sum(st_list) )/len(st_list)
    # print("E_T_m_m_1= {:.3f}, E_T_f_j_2= {:.3f}".format(E_T_m_m_1, E_T_f_j_2) )
    # print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
    diff_list.append(abs(E_T - sim_E_T) )
    return sim_E_T
  # print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )
  
  print("\n")
  # print("a_q.join_sink.qid__num_win_map= {}".format(pprint.pformat(a_q.join_sink.qid__num_win_map) ) )
  total_num_wins = sum([n for i, n in a_q.join_sink.qid__num_win_map.items() ] )
  print("pg.n_sent= {}, total_num_wins= {}".format(pg.n_sent, total_num_wins) )
  qid__win_freq_map = {i:float(n)/total_num_wins for i, n in a_q.join_sink.qid__num_win_map.items() }
  print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  
  print("\n")
  # print("a_q.join_q.state__num_found_map= {}".format(pprint.pformat(a_q.join_q.state__num_found_map) ) )
  # total_num_founds = sum([n for s, n in a_q.join_q.state__num_found_map.items() ] )
  # state__found_freq_map = {s:float(n)/total_num_founds for s, n in a_q.join_q.state__num_found_map.items() }
  # print("state__found_freq_map= {}".format(pprint.pformat(state__found_freq_map) ) )
  
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
  
def plot_winning_freqs():
  k, r, t = 2, 2, 1
  serv_rate = 1.0
  print("plot_winning_freqs:: k= {}, r= {}, t= {}, serv_rate= {}".format(k, r, t, serv_rate) )
  arr_list = []
  latex_symbol_list = ["\\gamma", "\\alpha", "\\beta"]
  qid__latex_symbol_map = {}
  qid__win_freq_list_map = {}
  # for arr_rate in numpy.arange(0.1, 1.05, 0.1):
  # for arr_rate in numpy.arange(0.1, 0.25, 0.1):
  for arr_rate in numpy.arange(0.05, 1.25, 0.05):
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = 1 + t*r
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
    a_q = AVQ("a_q", env, k, r, t, qid_list, qserv_dist_list)
    
    pg.out = a_q
    # aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 1)
    
    # env.run(until=500)
    env.run(until=50000)
    
    total_num_wins = sum([n for i, n in a_q.join_sink.qid__num_win_map.items() ] )
    qid__win_freq_map = {i:float(n)/total_num_wins for i, n in a_q.join_sink.qid__num_win_map.items() }
    print("arr_rate= {}, qid__win_freq_map= {}".format(arr_rate, pprint.pformat(qid__win_freq_map) ) )
    
    arr_list.append(arr_rate)
    q_counter = 0
    for qid, win_freq in qid__win_freq_map.items():
      if qid not in qid__win_freq_list_map:
        qid__win_freq_list_map[qid] = []
        qid__latex_symbol_map[qid] = latex_symbol_list[q_counter]
        q_counter += 1
      qid__win_freq_list_map[qid].append(win_freq)
    # print("qid__latex_symbol_map= {}".format(pprint.pformat(qid__latex_symbol_map) ) )
  
  ax = plot.gca()
  ax.grid(True)
  ax.set_prop_cycle(cycler('color', ['k', 'r', 'b'] ) ) # + cycler('lw', [1, 1, 1, 1] )
  for qid, win_freq_list in qid__win_freq_list_map.items():
    plot.plot(arr_list, win_freq_list, 'o-', label=r'$w_{}$'.format(qid__latex_symbol_map[qid] ) )
    plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("Winning frequency")
  plot.title("Fraction of the jobs terminated by each server\nfor increasing arrival rate, service rate: {}".format(serv_rate) )
  plot.savefig("plot_winning_freqs.png")

def plot_a_q():
  num_q = 3 # 3
  k, r = 2, 2 # 2, 2
  serv_rate = 1.0 # 5.0 # 1.0
  arr_list, p_lambda_list = [], []
  sim_E_T_mds_n_k_list, sim_E_T_simplex_list, E_T_simplex_list, sim_E_T_mds_n_1_list = [], [], [], []
  # for arr_rate in numpy.arange(0.1, 1.05, 0.1):
  for arr_rate in numpy.arange(0.1, 1.25, 0.1):
    arr_list.append(arr_rate)
    sim_E_T_mds_n_k = test_mds_n_k(arr_rate, serv_rate, num_q, k) # 0 # test_mds_n_k(arr_rate, serv_rate, num_q, k)
    sim_E_T_mds_n_1 = test_mds_n_1(arr_rate, serv_rate, num_q) # 0 # test_mds_n_1(arr_rate, serv_rate, num_q)
    # E_T_simplex = float(1/num_q)*sim_E_T_mds_n_1 + (1-float(1/num_q) )*sim_E_T_mds_n_k
    # E_T_simplex = 2/(3*serv_rate) + 3/5*sim_E_T_mds_n_1 + 2/5*sim_E_T_mds_n_k
    E_S_p = 0.5/serv_rate
    E_S_c = 0.665/serv_rate
    E_S = 0.4*E_S_p + 0.6*E_S_c
    E_S_p_2 = 0.5/(serv_rate**2)
    E_S_c_2 = 7/(9*(serv_rate**2) )
    E_S_2 = 0.4*E_S_p_2 + 0.6*E_S_c_2
    # E_T_simplex = E_S/(1-arr_rate*E_S) # w/ M/M/1 assumption
    E_T_simplex = E_S + (arr_rate/2)*E_S_2/(1-arr_rate*E_S) # w/ M/G/1 assumption
    sim_E_T_simplex = test_simplex_q(arr_rate, serv_rate, k, r)
    
    # p_lambda = (sim_E_T_mds_n_k-sim_E_T_simplex)/(sim_E_T_mds_n_k-sim_E_T_mds_n_1)
    # p_lambda_list.append(p_lambda)
    sim_E_T_mds_n_k_list.append(sim_E_T_mds_n_k)
    E_T_simplex_list.append(E_T_simplex)
    sim_E_T_simplex_list.append(sim_E_T_simplex)
    sim_E_T_mds_n_1_list.append(sim_E_T_mds_n_1)
  
  # plot.plot(qm.t_list, qm.n_list, 'r-')
  # plot.figure(1)
  # plot.subplot(211)
  plot.plot(arr_list, sim_E_T_mds_n_k_list, 'ro', label="mds_{}_{}".format(num_q, k) )
  plot.legend()
  plot.plot(arr_list, E_T_simplex_list, 'go', label="model_simplex")
  plot.legend()
  plot.plot(arr_list, sim_E_T_simplex_list, 'ko', label="simplex")
  plot.legend()
  plot.plot(arr_list, sim_E_T_mds_n_1_list, 'bo', label="mds_{}_1".format(num_q) )
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'Average response time vs. arrival rate, $\mu$= {}'.format(serv_rate) )
  # plot.subplot(212)
  # plot.plot(arr_list, p_lambda_list, 'b-')
  # plot.legend()
  # plot.xlabel(r'$\lambda$')
  # plot.ylabel("p(lambda)")
  # plot.title("")
  # plot.tight_layout()
  # ax = plot.gca()
  # ax.grid(True)
  plot.savefig("plot_a_q.png")
  
def rand_test():
  ar, sr = 0.9, 1.0 # 11.0, 10.0 # 0.9, 1.0
  k, r = 2, 2
  # for arr_rate in numpy.arange(0.1, 1.0, 0.1):
  #   test_simplex_q(arr_rate, sr, k=k, r=r)
  # test_simplex_q(ar, sr, k=k, r=r)
  # 0.75 1st wins
  test_mds_n_k(ar, sr, k+1, k)
  # test_mds_n_1(ar, sr, k+1)
  
if __name__ == "__main__":
  ar, sr = 2, 1.0 # 1.1, 1.0 # 0.9, 1.0
  # eg_pgen_psink()
  # eg_overloaded_m1_q()
  # test_m_m_1()
  # test_fj()
  # test_mds_n_1(ar, sr, 3)
  # test_mds_n_k(ar, sr, 3, 2)
  # test_simplex_q(ar, sr, k=2, r=2)
  # plot_winning_freqs()
  plot_a_q()
  # rand_test()
  
  
