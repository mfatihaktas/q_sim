import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
from cycler import cycler
from random import expovariate
import sys, pprint, math, numpy, simpy, getopt

from sim_components import *

plot_color_list = ["indigo", "darkorange", "yellowgreen", "cyan", "darkmagenta", "darkred", "black", "slateblue", "goldenrod", "darksalmon", "forestgreen", "saddlebrown", "grey"]

def test_mds_n_k(arr_rate, serv_rate, num_q, k):
  print("test_mds_n_k:: arr_rate= {}, serv_rate= {}, num_q= {}, k= {}"
        .format(arr_rate, serv_rate, num_q, k) )
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
  
  # E_T = 1.0/(num_q*serv_rate - arr_rate)
  st_list = mdsq.join_sink.st_list
  if len(st_list) > 0:
    sim_E_T = float(sum(st_list) )/len(st_list)
    # print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, sim_E_T) )
    # diff_list.append(abs(E_T - sim_E_T) )
    return sim_E_T
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
  mu = 1.0
  # Inner bound on the arr_rate for stability
  def beta(x, y):
    return math.gamma(x)*math.gamma(y)/math.gamma(x+y)
  E_S_split_merge = 1/(mu*r)*beta(t+1, 1/r)
  arr_rate_ub = 1/E_S_split_merge
  print("plot_winning_freqs:: k= {}, r= {}, t= {}, mu= {}, arr_rate_ub={}".format(k, r, t, mu, arr_rate_ub) )
  arr_list = []
  qid__win_freq_list_map = {}
  # for arr_rate in numpy.arange(0.1, 1.05, 0.1):
  # for arr_rate in numpy.arange(0.1, 2.0, 0.2):
  # for arr_rate in numpy.arange(0.05, arr_rate_ub*1.25, 0.05):
  for arr_rate in numpy.arange(0.05, arr_rate_ub*1.1, arr_rate_ub/10):
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = 1 + t*r
    qid_list = ["{}".format(i) for i in range(1, num_q + 1) ]
    qserv_dist_list = [lambda: random.expovariate(mu) for i in range(num_q) ]
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
        q_counter += 1
      qid__win_freq_list_map[qid].append(win_freq)
  
  ax = plot.gca()
  ax.grid(True)
  ax.set_prop_cycle(cycler('color', ['k', 'r', 'b'] ) ) # + cycler('lw', [1, 1, 1, 1] )
  qid__latex_symbol_map = {"1":"s", "2":"{1,1}", "3":"{1,2}"} # ["\\gamma", "\\alpha", "\\beta"]
  for qid, win_freq_list in qid__win_freq_list_map.items():
    plot.plot(arr_list, win_freq_list, 'o-', label=r'$w_{}$'.format(qid__latex_symbol_map[qid] ) )
    plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("Winning frequency")
  plot.title(r'$\mu$= {}'.format(mu) )
  plot.savefig("plot_winning_freqs_avq.png")

def mds_n_2_E_T(arr_rate, mu, n): # for only k = 2
  # high-regime assumption
  f_jc = (n - 2)/(n - 1)
  f_jp = 1 - f_jc
  mds_E_S_p = 1/((n - 1)*mu)
  mds_E_S_c = (1/mu)*(harmonic_sum(n) - harmonic_sum(n-2) )
  mds_E_S = f_jp*mds_E_S_p + f_jc*mds_E_S_c
  mds_E_S_p_2 = 2/(((n - 1)**2)*(mu**2) )
  mds_E_S_c_2 = (1/(mu**2) )*(harmonic_2_sum(n) - harmonic_2_sum(n-2) ) + mds_E_S_c**2
  mds_E_S_2 = f_jp*mds_E_S_p_2 + f_jc*mds_E_S_c_2
  # mds_E_T = mds_E_S/(1 - arr_rate*mds_E_S) # w/ M/M/1 assumption
  mds_E_T = mds_E_S + (arr_rate/2)*mds_E_S_2/(1 - arr_rate*mds_E_S) # w/ M/G/1 assumption
  return mds_E_T
  
def adjustable_mds_n_2_E_T(arr_rate, mu, n): # for only k = 2
  # gamma = min(arr_rate, mu)
  # gamma = min(arr_rate*mu/(arr_rate+mu), \
  #             mu/(harmonic_sum(n) - harmonic_sum(n-2) ) )
  gamma = mu
  ro = gamma/mu
  p_0 = 1/(1 + n*ro/(n - 1 - ro) )
  def p_i(i):
    return (ro/(n - 1) )**i * n*p_0
  
  nu = (n-1)*mu + gamma
  sum_for_pi = p_0*n*gamma + nu*(1-p_0)
  pi_0 = p_0*n*gamma/sum_for_pi
  pi_1 = p_i(1)*nu/sum_for_pi
  
  f_jd = (n-1)*mu/nu*(1-pi_0)
  f_c = pi_1*(n-1)*mu/nu
  f_jc = f_c/f_jd
  f_jp = 1 - f_jc
  
  mds_E_S_p = 1/((n - 1)*mu)
  mds_E_S_c = (1/mu)*(harmonic_sum(n) - harmonic_sum(n-2) )
  mds_E_S = f_jp*mds_E_S_p + f_jc*mds_E_S_c
  mds_E_S_p_2 = 2/(((n - 1)**2)*(mu**2) )
  mds_E_S_c_2 = (1/(mu**2) )*(harmonic_2_sum(n) - harmonic_2_sum(n-2) ) + mds_E_S_c**2
  mds_E_S_2 = f_jp*mds_E_S_p_2 + f_jc*mds_E_S_c_2
  mds_E_T = mds_E_S + (arr_rate/2)*mds_E_S_2/(1 - arr_rate*mds_E_S)
  return mds_E_T

def split_merge_mds_n_k_E_T(arr_rate, mu, n, k): # works for k >= 1
  harmonic_diff = harmonic_sum(n) - harmonic_sum(n-k)
  if 1/arr_rate < harmonic_diff/mu:
    return None
  ro = arr_rate/mu
  return harmonic_diff/mu + \
    arr_rate*(harmonic_2_sum(n)-harmonic_2_sum(n-k) + harmonic_diff**2)/(2*mu**2 * (1 - ro*harmonic_diff) )

def adj_split_merge_mds_n_k_E_T(arr_rate, mu, n, k):
  # gamma = arr_rate
  # p_k_1 = 1/((n-k+1)*mu/gamma * (harmonic_sum(n) - harmonic_sum(n-k+1) ) + 1)
  # avg_active_n = (n-k+1)*(k-1)*mu/gamma * p_k_1
  # mu_ = n/avg_active_n * mu
  # return split_merge_mds_n_k_E_T(arr_rate, mu_, n, k)
  # return split_merge_mds_n_k_E_T(arr_rate, mu, n, k-1) + 1/((n-k+1)*mu)
  return split_merge_mds_n_k_E_T(arr_rate, mu, n, k-1) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)

def adj_2_split_merge_mds_n_k_E_T(arr_rate, mu, n, k):
  # return split_merge_mds_n_k_E_T(arr_rate, mu, n, k-2) + 1/((n-k+1)*mu) + 1/((n-k+1)*mu)
  return split_merge_mds_n_k_E_T(arr_rate, mu, n, k-2) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+2, 1) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)

def recur_split_merge_mds_n_k_E_T(arr_rate, mu, n, k):
  if n > k:
    if k > 1:
      # return recur_split_merge_mds_n_k_E_T(arr_rate, mu, n, k-1) + 1/((n-k+1)*mu)
      return recur_split_merge_mds_n_k_E_T(arr_rate, mu, n, k-1) + split_merge_mds_n_k_E_T(arr_rate, mu, n-k+1, 1)
    elif k == 1:
      return split_merge_mds_n_k_E_T(arr_rate, mu, n, k)
    else:
      log(ERROR, "Unexpected k= {}".format(k) )
      sys.exit(1)
  else:
    log(ERROR, "Unexpected n= {} <= k= {}".format(n, k) )
    return split_merge_mds_n_k_E_T(arr_rate, mu, n, k)

# def w_exact_start_setup_prob_simplex_E_T(arr_rate, mu):
#   simplex_E_S_p = 0.5/mu
#   simplex_E_S_c = 0.665/mu
#   # computing p_c:Pr{a job makes a complete-start}
#   E_X = 1/arr_rate
#   delta = (E_X-simplex_E_S_p)**2 - 4*()

def inner_bound_on_arr_rate(n, k, r, t, mu):
  def beta(x, y):
    return math.gamma(x)*math.gamma(y)/math.gamma(x+y)
  simplex_E_S_split_merge = 1/(mu*r)*beta(t+1, 1/r)
  
  mds_n_k_E_S_split_merge = 1/mu * (harmonic_sum(n) - harmonic_sum(n-k) )
  
  # return min(1/simplex_E_S_split_merge, 1/mds_n_k_E_S_split_merge)
  return 1/mds_n_k_E_S_split_merge

def plot_a_q(num_q):
  k, r, t = 2, 2, 1 # 2, 2
  mu = 1.0 # 5.0 # 1.0
  arr_rate_ub = inner_bound_on_arr_rate(num_q, k, r, t, mu)
  log(WARNING, "k= {}, r= {}, t= {}, mu= {}, arr_rate_ub={}".format(k, r, t, mu, arr_rate_ub) )
  
  arr_list, p_lambda_list = [], []
  fj_k_k_E_T_list = []
  varki_mds_n_k_E_T_list, mds_n_k_E_T_list, adj_mds_n_k_E_T_list, sim_mds_n_k_E_T_list = [], [], [], []
  split_merge_mds_n_k_E_T_list, adj_split_merge_mds_n_k_E_T_list, adj_2_split_merge_mds_n_k_E_T_list = [], [], []
  recur_split_merge_mds_n_k_E_T_list = []
  simplex_E_T_list, sim_simplex_E_T_list, sim_mds_n_1_E_T_list = [], [], []
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, 0.1):
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
    arr_list.append(arr_rate)
    # varki
    ro = float(arr_rate/mu)
    varki_mds_n_k_E_T = 1/mu * (harmonic_sum(num_q) - harmonic_sum(num_q-k) ) + \
      1/mu * ro*(gen_harmonic_sum(num_q, ro) - gen_harmonic_sum(num_q-k, ro) )
    varki_mds_n_k_E_T_list.append(varki_mds_n_k_E_T)
    
    mds_n_k_E_T_list.append(mds_n_2_E_T(arr_rate, mu, num_q) )
    adj_mds_n_k_E_T_list.append(adjustable_mds_n_2_E_T(arr_rate, mu, num_q) )
    
    split_merge_mds_n_k_E_T_list.append(split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    adj_split_merge_mds_n_k_E_T_list.append(adj_split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    adj_2_split_merge_mds_n_k_E_T_list.append(adj_2_split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    recur_split_merge_mds_n_k_E_T_list.append(recur_split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    # sim
    sim_mds_n_k_E_T_f_sum = 0
    for f in range(3):
      sim_mds_n_k_E_T_f_sum += test_mds_n_k(arr_rate, mu, num_q, k)
    sim_mds_n_k_E_T = sim_mds_n_k_E_T_f_sum/3
    # sim_mds_n_k_E_T = test_mds_n_k(arr_rate, mu, num_q, k)
    sim_mds_n_k_E_T_list.append(sim_mds_n_k_E_T)
    fj_k_k_E_T = 0 # test_mds_n_k(arr_rate, mu, k, k) # 0
    fj_k_k_E_T_list.append(fj_k_k_E_T)
    
    sim_mds_n_1_E_T = test_mds_n_k(arr_rate, mu, num_q, 1)
    sim_mds_n_1_E_T_list.append(sim_mds_n_1_E_T)
    
    # simplex_E_T = float(1/num_q)*sim_mds_n_1_E_T + (1-float(1/num_q) )*sim_mds_n_k_E_T
    # simplex_E_T = 2/(3*mu) + 3/5*sim_mds_n_1_E_T + 2/5*sim_mds_n_k_E_T
    simplex_E_S_p = 0.5/mu
    simplex_E_S_c = 0.665/mu
    simplex_E_S = 0.4*simplex_E_S_p + 0.6*simplex_E_S_c
    simplex_E_S_p_2 = 0.5/(mu**2)
    simplex_E_S_c_2 = 7/(9*(mu**2) )
    simplex_E_S_2 = 0.4*simplex_E_S_p_2 + 0.6*simplex_E_S_c_2
    # simplex_E_T = simplex_E_S/(1 - arr_rate*simplex_E_S) # w/ M/M/1 assumption
    simplex_E_T = simplex_E_S + (arr_rate/2)*simplex_E_S_2/(1 - arr_rate*simplex_E_S) # w/ M/G/1 assumption
    simplex_E_T_list.append(simplex_E_T)
    # sim_simplex_E_T = test_simplex_q(arr_rate, mu, k, r)
    # sim_simplex_E_T_f_sum = 0
    # for f in range(3):
    #   sim_simplex_E_T_f_sum += test_simplex_q(arr_rate, mu, k, r)
    sim_simplex_E_T = 0 # sim_simplex_E_T_f_sum/3
    sim_simplex_E_T_list.append(sim_simplex_E_T)
    # p_lambda = (sim_mds_n_k_E_T-sim_simplex_E_T)/(sim_mds_n_k_E_T-simmds_n_1_E_T)
    # p_lambda_list.append(p_lambda)
    
  # plot.plot(qm.t_list, qm.n_list, 'r-')
  # plot.figure(1)
  # plot.subplot(211)
  plot.plot(arr_list, varki_mds_n_k_E_T_list, 'mo', label="varki_mds_{}_{}".format(num_q, k) )
  plot.legend()
  # plot.plot(arr_list, mds_n_k_E_T_list, 'go', label="model_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  plot.plot(arr_list, adj_mds_n_k_E_T_list, 'go', label="model_mds_{}_{}".format(num_q, k) )
  plot.legend()
  
  plot.plot(arr_list, split_merge_mds_n_k_E_T_list, 'ro', label="sm_mds_{}_{}".format(num_q, k) )
  plot.legend()
  # plot.plot(arr_list, adj_split_merge_mds_n_k_E_T_list, 'bo', label="adj_sm_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  # plot.plot(arr_list, adj_2_split_merge_mds_n_k_E_T_list, 'co', label="adj_2_sm_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  # plot.plot(arr_list, recur_split_merge_mds_n_k_E_T_list, 'yo', label="recur_sm_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  
  plot.plot(arr_list, sim_mds_n_k_E_T_list, 'ko', label="mds_{}_{}".format(num_q, k) )
  plot.legend()
  # plot.plot(arr_list, sim_simplex_E_T_list, 'ko', label="simplex")
  plot.legend()
  # plot.plot(arr_list, simplex_E_T_list, 'go', label="model_simplex")
  # plot.legend()
  # plot.plot(arr_list, fj_k_k_E_T_list, 'bo', label="fj_{}_{}".format(k, k) )
  # plot.legend()
  # plot.plot(arr_list, sim_mds_n_1_E_T_list, 'bo', label="mds_{}_1".format(num_q) )
  # plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, k= {}, homogeneous servers $\mu$= {}'.format(num_q, k, mu) )
  # plot.subplot(212)
  # plot.plot(arr_list, p_lambda_list, 'b-')
  # plot.legend()
  # plot.xlabel(r'$\lambda$')
  # plot.ylabel("p(lambda)")
  # plot.title("")
  # plot.tight_layout()
  # ax = plot.gca()
  # ax.grid(True)
  plot.savefig("plot_a_q__n_{}_k_{}.png".format(num_q, k) )
  log(WARNING, "done.")

def plot_complex_simplex(num_q):
  k, r, t = 2, 2, 2 # 2, 2
  mu = 1.0
  arr_rate_ub = inner_bound_on_arr_rate(num_q, k, r, t, mu)
  print("plot_a_q:: k= {}, r= {}, t= {}, mu= {}, arr_rate_ub={}".format(k, r, t, mu, arr_rate_ub) )
  
  arr_list, p_lambda_list = [], []
  fj_k_k_E_T_list = []
  varki_mds_n_k_E_T_list, mds_n_k_E_T_list, adj_mds_n_k_E_T_list, sim_mds_n_k_E_T_list = [], [], [], []
  split_merge_mds_n_k_E_T_list, adj_split_merge_mds_n_k_E_T_list, adj_2_split_merge_mds_n_k_E_T_list = [], [], []
  recur_split_merge_mds_n_k_E_T_list = []
  simplex_E_T_list, sim_simplex_E_T_list, sim_mds_n_1_E_T_list = [], [], []
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, 0.1):
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
    arr_list.append(arr_rate)
    # varki
    ro = float(arr_rate/mu)
    varki_mds_n_k_E_T = 1/mu * (harmonic_sum(num_q) - harmonic_sum(num_q-k) ) + \
      1/mu * ro*(gen_harmonic_sum(num_q, ro) - gen_harmonic_sum(num_q-k, ro) )
    varki_mds_n_k_E_T_list.append(varki_mds_n_k_E_T)
    
    mds_n_k_E_T_list.append(mds_n_2_E_T(arr_rate, mu, num_q) )
    adj_mds_n_k_E_T_list.append(adjustable_mds_n_2_E_T(arr_rate, mu, num_q) )
    
    split_merge_mds_n_k_E_T_list.append(split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    adj_split_merge_mds_n_k_E_T_list.append(adj_split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    adj_2_split_merge_mds_n_k_E_T_list.append(adj_2_split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    recur_split_merge_mds_n_k_E_T_list.append(recur_split_merge_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    # sim
    sim_mds_n_k_E_T_f_sum = 0
    for f in range(3):
      sim_mds_n_k_E_T_f_sum += test_mds_n_k(arr_rate, mu, num_q, k)
    sim_mds_n_k_E_T = sim_mds_n_k_E_T_f_sum/3
    # sim_mds_n_k_E_T = test_mds_n_k(arr_rate, mu, num_q, k)
    sim_mds_n_k_E_T_list.append(sim_mds_n_k_E_T)
    fj_k_k_E_T = 0 # test_mds_n_k(arr_rate, mu, k, k) # 0
    fj_k_k_E_T_list.append(fj_k_k_E_T)
    
    sim_mds_n_1_E_T = test_mds_n_k(arr_rate, mu, num_q, 1)
    sim_mds_n_1_E_T_list.append(sim_mds_n_1_E_T)
    
    # simplex_E_T = float(1/num_q)*sim_mds_n_1_E_T + (1-float(1/num_q) )*sim_mds_n_k_E_T
    # simplex_E_T = 2/(3*mu) + 3/5*sim_mds_n_1_E_T + 2/5*sim_mds_n_k_E_T
    simplex_E_S_p = 0.5/mu
    simplex_E_S_c = 0.665/mu
    simplex_E_S = 0.4*simplex_E_S_p + 0.6*simplex_E_S_c
    simplex_E_S_p_2 = 0.5/(mu**2)
    simplex_E_S_c_2 = 7/(9*(mu**2) )
    simplex_E_S_2 = 0.4*simplex_E_S_p_2 + 0.6*simplex_E_S_c_2
    # simplex_E_T = simplex_E_S/(1 - arr_rate*simplex_E_S) # w/ M/M/1 assumption
    simplex_E_T = simplex_E_S + (arr_rate/2)*simplex_E_S_2/(1 - arr_rate*simplex_E_S) # w/ M/G/1 assumption
    simplex_E_T_list.append(simplex_E_T)
    # sim_simplex_E_T = test_simplex_q(arr_rate, mu, k, r)
    # sim_simplex_E_T_f_sum = 0
    # for f in range(3):
    #   sim_simplex_E_T_f_sum += test_simplex_q(arr_rate, mu, k, r)
    sim_simplex_E_T = 0 # sim_simplex_E_T_f_sum/3
    sim_simplex_E_T_list.append(sim_simplex_E_T)
    # p_lambda = (sim_mds_n_k_E_T-sim_simplex_E_T)/(sim_mds_n_k_E_T-simmds_n_1_E_T)
    # p_lambda_list.append(p_lambda)
    
  # plot.plot(qm.t_list, qm.n_list, 'r-')
  # plot.figure(1)
  # plot.subplot(211)
  plot.plot(arr_list, varki_mds_n_k_E_T_list, 'mo', label="varki_mds_{}_{}".format(num_q, k) )
  plot.legend()
  # plot.plot(arr_list, mds_n_k_E_T_list, 'go', label="model_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  plot.plot(arr_list, adj_mds_n_k_E_T_list, 'go', label="model_mds_{}_{}".format(num_q, k) )
  plot.legend()
  
  plot.plot(arr_list, split_merge_mds_n_k_E_T_list, 'ro', label="sm_mds_{}_{}".format(num_q, k) )
  plot.legend()
  # plot.plot(arr_list, adj_split_merge_mds_n_k_E_T_list, 'bo', label="adj_sm_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  # plot.plot(arr_list, adj_2_split_merge_mds_n_k_E_T_list, 'co', label="adj_2_sm_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  # plot.plot(arr_list, recur_split_merge_mds_n_k_E_T_list, 'yo', label="recur_sm_mds_{}_{}".format(num_q, k) )
  # plot.legend()
  
  plot.plot(arr_list, sim_mds_n_k_E_T_list, 'ko', label="mds_{}_{}".format(num_q, k) )
  plot.legend()
  # plot.plot(arr_list, sim_simplex_E_T_list, 'ko', label="simplex")
  plot.legend()
  # plot.plot(arr_list, simplex_E_T_list, 'go', label="model_simplex")
  # plot.legend()
  # plot.plot(arr_list, fj_k_k_E_T_list, 'bo', label="fj_{}_{}".format(k, k) )
  # plot.legend()
  # plot.plot(arr_list, sim_mds_n_1_E_T_list, 'bo', label="mds_{}_1".format(num_q) )
  # plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, k= {}, homogeneous servers $\mu$= {}'.format(num_q, k, mu) )
  # plot.subplot(212)
  # plot.plot(arr_list, p_lambda_list, 'b-')
  # plot.legend()
  # plot.xlabel(r'$\lambda$')
  # plot.ylabel("p(lambda)")
  # plot.title("")
  # plot.tight_layout()
  # ax = plot.gca()
  # ax.grid(True)
  plot.savefig("plot_a_q__n_{}_k_{}.png".format(num_q, k) )
  log(WARNING, "done.")

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

def rand_test():
  ar, sr = 0.9, 1.0 # 11.0, 10.0 # 0.9, 1.0
  k, r = 2, 2
  # for arr_rate in numpy.arange(0.1, 1.0, 0.1):
  #   test_simplex_q(arr_rate, sr, k=k, r=r)
  # test_simplex_q(ar, sr, k=k, r=r)
  # 0.75 1st wins
  test_mds_n_k(ar, sr, k+1, k)
  # test_mds_n_k(ar, sr, k+1)

def get_opts(argv):
  opt_map = {}
  try:
    opts, args = getopt.getopt(argv, '', ['num_q='] )
  except getopt.GetoptError:
    log(ERROR, "Unexpected command line arg, expecting: exp.py --num_q=<>")
    sys.exit(1)
  
  for opt, arg in opts:
    opt_map[opt] = arg
  return opt_map

if __name__ == "__main__":
  opt_map = get_opts(sys.argv[1:] )
  log(WARNING, "opt_map= {}".format(pprint.pformat(opt_map) ) )
  num_q = int(opt_map["--num_q"] )
  ar, sr = 1.0, 1.0 # 1.1, 1.0 # 0.9, 1.0
  # eg_pgen_psink()
  # eg_overloaded_m1_q()
  # test_m_m_1()
  # test_fj()
  # test_mds_n_k(ar, sr, 3)
  # test_mds_n_k(ar, sr, 3, 2)
  # test_simplex_q(ar, sr, k=2, r=2)
  # plot_winning_freqs()
  plot_a_q(num_q)
  # prob_comp_start()
  # plot_split_merge_approx_error(ar, sr, num_q)
  # rand_test()
  
  # sm_E_T = split_merge_mds_n_k_E_T(arr_rate=1, mu=1, n=num_q, k=4)
  # print("split_merge_mds_{}_{}_E_T= {}".format(num_q, 4, sm_E_T) )
  
