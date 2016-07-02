import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from random import expovariate
import pprint
import simpy

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
  arr_rate = 0.9
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

def test_mds_n_1():
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    arr_rate = 0.5
    serv_rate = 1.0
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = 5
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
    mds_q = MDSQ("mds_q", env, 1, qid_list, qserv_dist_list)
    # qm = QMonitor(env, q=m1_q, dist=lambda: 500)
    
    pg.out = mds_q
    
    # env.run(until=5)
    # env.run(until=100)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    print("num_q= {}, arr_rate= {}, serv_rate= {}".format(num_q, arr_rate, serv_rate) )
    # ro = arr_rate/serv_rate
    E_T = 1.0/(num_q*serv_rate - arr_rate)
    st_list = mds_q.join_sink.st_list
    if len(st_list) > 0:
      sim_E_T = float(sum(st_list) )/len(st_list)
      print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
      diff_list.append(abs(E_T - sim_E_T) )
  print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )

"""
def test_mds_n_k():
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    arr_rate = 0.5
    serv_rate = 1.0
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = 5
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
    mds_q = MDSQ("mds_q", env, 1, qid_list, qserv_dist_list)
    # qm = QMonitor(env, q=m1_q, dist=lambda: 500)
    
    pg.out = mds_q
    
    # env.run(until=5)
    # env.run(until=100)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    print("num_q= {}, arr_rate= {}, serv_rate= {}".format(num_q, arr_rate, serv_rate) )
    # ro = arr_rate/serv_rate
    E_T = 1.0/(num_q*serv_rate - arr_rate)
    st_list = mds_q.join_sink.st_list
    if len(st_list) > 0:
      sim_E_T = float(sum(st_list) )/len(st_list)
      print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
      diff_list.append(abs(E_T - sim_E_T) )
  print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )
"""

def test_simplex_q():
  diff_list = []
  for c in range(1):
    env = simpy.Environment()
    arr_rate = 0.9 # 0.5
    serv_rate = 1.0
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    k, r, t = 2, 2, 1
    num_q = 1 + t*r
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(serv_rate) for i in range(num_q) ]
    a_q = AQ("a_q", env, k, r, t, qid_list, qserv_dist_list)
    
    pg.out = a_q
    
    # env.run(until=5)
    env.run(until=50000)
    # env.run(until=200000)
    
    # for num_q= 2
    print("num_q= {}, arr_rate= {}, serv_rate= {}".format(num_q, arr_rate, serv_rate) )
    print("a_q= {}".format(a_q) )
    ro = arr_rate/serv_rate
    E_T_m_m_1 = 1.0/(serv_rate-arr_rate)
    E_T_f_j_2 = (12-ro)/8*E_T_m_m_1
    st_list = a_q.join_sink.st_list
    E_T = 7/9*E_T_m_m_1 + 2/9*E_T_f_j_2
    if len(st_list) > 0:
      sim_E_T = float(sum(st_list) )/len(st_list)
      print("E_T_m_m_1= {:.3f}, E_T_f_j_2= {:.3f}".format(E_T_m_m_1, E_T_f_j_2) )
      print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
      diff_list.append(abs(E_T - sim_E_T) )
  print("diff_list= [{}]".format("".join("%s, " % d for d in diff_list) ) )
  
  print("a_q.join_sink.qid__num_win_map= {}".format(pprint.pformat(a_q.join_sink.qid__num_win_map) ) )
  total_num_wins = sum([n for i, n in a_q.join_sink.qid__num_win_map.items() ] )
  print("pg.n_sent= {}, total_num_wins= {}".format(pg.n_sent, total_num_wins) )
  qid__win_freq_map = {i:float(n)/total_num_wins for i, n in a_q.join_sink.qid__num_win_map.items() }
  print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  
  print("a_q.group_id__q_map= {}".format(pprint.pformat(a_q.group_id__q_map) ) )

if __name__ == "__main__":
  # eg_pgen_psink()
  # eg_overloaded_m1_q()
  # test_m_m_1()
  # test_fj()
  # test_mds_n_1()
  test_simplex_q()
  