import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
from cycler import cycler
from random import expovariate
import sys, pprint, math, numpy, simpy, getopt, itertools

from sim_components import *
from simplex_models import *
# plot_color_l = ["indigo", "darkorange", "yellowgreen", "cyan", "darkmagenta", "darkred", "black", "slateblue", "goldenrod", "darksalmon", "forestgreen", "saddlebrown", "grey"]

def test_mds_n_k(num_f_run, arr_rate, mu, n, k):
  sim_E_T_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "arr_rate= {}, mu= {}, n= {}, k= {}".format(arr_rate, mu, n, k) )
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    qid_l = ["{}".format(i) for i in range(1, n + 1) ]
    # qserv_dist_l = [lambda: random.expovariate(mu) for i in range(n) ]
    mdsq = MDSQ("mdsq", env, k, qid_l, qserv_rate_l=[mu for i in range(n) ] )
    mdsq_monitor = MDSQMonitor(env, q=mdsq, poll_dist=lambda: 1)
    pg.out = mdsq
    env.run(until=50000) # env.run(until=5000)
    
    st_l = mdsq.join_sink.st_l
    if len(st_l) > 0:
      sim_E_T_f_sum += float(sum(st_l) )/len(st_l)
      continue
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
  return sim_E_T_f_sum/num_f_run
  
def test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l=[], w_sys=True):
  sim_E_T_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "w_sys= {}, arr_rate= {}, mu= {}, k= {}, r= {}, t= {}, qmu_l= {}". \
      format(w_sys, arr_rate, mu, k, r, t, pprint.pformat(qmu_l) ) )
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = int(1 + t*r) if w_sys else int(t*r)
    qid_l = ["{}".format(i) for i in range(1, num_q + 1) ]
    if len(qmu_l) == 0:
      qmu_l = [mu for i in range(num_q) ]
    a_q = AVQ("a_q", env, k, r, t, qid_l, qmu_l, w_sys=w_sys)
    aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 0.1)
    a_q.join_q.out_m = aq_monitor
    pg.out = a_q
    env.run(until=50000) # env.run(until=500)
    
    st_l = a_q.join_sink.st_l
    if len(st_l) > 0:
      sim_E_T_f_sum += float(sum(st_l) )/len(st_l)
      # continue
    # print("a_q.join_sink.qid__num_win_map= {}".format(pprint.pformat(a_q.join_sink.qid__num_win_map) ) )
    total_num_wins = sum([n for i, n in a_q.join_sink.qid__num_win_map.items() ] )
    print("pg.n_sent= {}, total_num_wins= {}".format(pg.n_sent, total_num_wins) )
    qid__win_freq_map = {i:float(n)/total_num_wins for i, n in a_q.join_sink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
    """
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
    """
  return sim_E_T_f_sum/num_f_run
  
def plot_winning_freqs():
  k, r, t = 2, 2, 1
  mu = 1.0
  # Inner bound on the arr_rate for stability
  def beta(x, y):
    return math.gamma(x)*math.gamma(y)/math.gamma(x+y)
  E_S_sm = 1/(mu*r)*beta(t+1, 1/r)
  arr_rate_ub = 1/E_S_sm
  print("plot_winning_freqs:: k= {}, r= {}, t= {}, mu= {}, arr_rate_ub={}".format(k, r, t, mu, arr_rate_ub) )
  arr_rate_l = []
  qid__win_freq_l_map = {}
  # for arr_rate in numpy.arange(0.1, 1.05, 0.1):
  # for arr_rate in numpy.arange(0.1, 2.0, 0.2):
  # for arr_rate in numpy.arange(0.05, arr_rate_ub*1.25, 0.05):
  for arr_rate in numpy.arange(0.05, arr_rate_ub*1.1, arr_rate_ub/10):
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="p_gen",
                         adist=lambda: random.expovariate(arr_rate),
                         sdist=lambda: 1)
    num_q = 1 + t*r
    qid_l = ["{}".format(i) for i in range(1, num_q + 1) ]
    a_q = AVQ("a_q", env, k, r, t, qid_l, qserv_dist_l=[mu for i in range(num_q) ] )
    
    pg.out = a_q
    # aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 1)
    
    # env.run(until=500)
    env.run(until=50000)
    
    total_num_wins = sum([n for i, n in a_q.join_sink.qid__num_win_map.items() ] )
    qid__win_freq_map = {i:float(n)/total_num_wins for i, n in a_q.join_sink.qid__num_win_map.items() }
    print("arr_rate= {}, qid__win_freq_map= {}".format(arr_rate, pprint.pformat(qid__win_freq_map) ) )
    
    arr_rate_l.append(arr_rate)
    q_counter = 0
    for qid, win_freq in qid__win_freq_map.items():
      if qid not in qid__win_freq_l_map:
        qid__win_freq_l_map[qid] = []
        q_counter += 1
      qid__win_freq_l_map[qid].append(win_freq)
  
  ax = plot.gca()
  ax.grid(True)
  ax.set_prop_cycle(cycler('color', ['k', 'r', 'b'] ) ) # + cycler('lw', [1, 1, 1, 1] )
  qid__latex_symbol_map = {"1":"s", "2":"{1,1}", "3":"{1,2}"} # ["\\gamma", "\\alpha", "\\beta"]
  for qid, win_freq_l in qid__win_freq_l_map.items():
    plot.plot(arr_rate_l, win_freq_l, 'o-', label=r'$w_{}$'.format(qid__latex_symbol_map[qid] ) )
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("Winning frequency")
  plot.title(r'$\mu$= {}'.format(mu) )
  plot.savefig("plot_winning_freqs_avq.png")

def plot_mds(num_q):
  k = 2
  mu = 1.0
  arr_rate_ub = mds_inner_bound_on_arr_rate(num_q, k, mu)
  log(WARNING, "k= {}, mu= {}, arr_rate_ub={}".format(k, mu, arr_rate_ub) )
  
  arr_rate_l, p_lambda_l = [], []
  sim_fj_k_k_E_T_l = []
  varki_mds_n_k_E_T_l, mds_n_k_E_T_l, adj_mds_n_k_E_T_l, sim_mds_n_k_E_T_l = [], [], [], []
  sm_mds_n_k_E_T_l, adj_sm_mds_n_k_E_T_l, adj_2_sm_mds_n_k_E_T_l = [], [], []
  recur_sm_mds_n_k_E_T_l = []
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, 0.1):
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
    arr_rate_l.append(arr_rate)
    
    sm_mds_n_k_E_T_l.append(sm_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    adj_sm_mds_n_k_E_T_l.append(adj_sm_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    adj_2_sm_mds_n_k_E_T_l.append(adj_2_sm_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    recur_sm_mds_n_k_E_T_l.append(recur_sm_mds_n_k_E_T(arr_rate, mu, num_q, k) )
    
    mds_n_k_E_T_l.append(mds_n_2_E_T(arr_rate, mu, num_q) )
    adj_mds_n_k_E_T_l.append(adjustable_mds_n_2_E_T(arr_rate, mu, num_q) )
    
    ro = float(arr_rate/mu)
    varki_mds_n_k_E_T = 1/mu * (harmonic_sum(num_q) - harmonic_sum(num_q-k) ) + \
      1/mu * ro*(gen_harmonic_sum(num_q, ro) - gen_harmonic_sum(num_q-k, ro) )
    varki_mds_n_k_E_T_l.append(varki_mds_n_k_E_T)
    # sim
    num_f_run = 1
    sim_mds_n_k_E_T = test_mds_n_k(num_f_run, arr_rate, mu, num_q, k)
    sim_mds_n_k_E_T_l.append(sim_mds_n_k_E_T)
    
    sim_fj_k_k_E_T = 0 #  test_mds_n_k(num_f_run, arr_rate, mu, k, k)
    sim_fj_k_k_E_T_l.append(sim_fj_k_k_E_T)
    
    sim_mds_n_1_E_T = 0 # test_mds_n_k(num_f_run, arr_rate, mu, num_q, 1)
    sim_mds_n_1_E_T_l.append(sim_mds_n_1_E_T)
  plot.plot(arr_rate_l, sm_mds_n_k_E_T_l, 'ro', label="SM-MDS({},{})".format(num_q, k) )
  # plot.plot(arr_rate_l, adj_sm_mds_n_k_E_T_l, 'bo', label="adj_sm_MDS({},{})".format(num_q, k) )
  # plot.plot(arr_rate_l, adj_2_sm_mds_n_k_E_T_l, 'co', label="adj_2_sm_MDS({},{})".format(num_q, k) )
  # plot.plot(arr_rate_l, recur_sm_mds_n_k_E_T_l, 'yo', label="recur_sm_MDS({},{})".format(num_q, k) )
  plot.plot(arr_rate_l, sim_mds_n_k_E_T_l, 'ko', label="MDS({},{})".format(num_q, k) )
  plot.plot(arr_rate_l, adj_mds_n_k_E_T_l, 'go', label="LB-MDS({},{})".format(num_q, k) )
  plot.plot(arr_rate_l, varki_mds_n_k_E_T_l, 'mo', label="Varki-MDS({},{})".format(num_q, k) )
  
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, k= {}, identical servers $\mu$= {}'.format(num_q, k, mu) )
  plot.savefig("plot_mds__n_{}_k_{}.png".format(num_q, k) )
  log(WARNING, "done; n= {}, k= {}".format(n, k) )

def plot_simplex(num_q):
  w_sys = True # False
  t, r, k = 1, 2, 2
  mu = 1
  gamma = mu
  arr_rate_ub = 1.65 # simplex_inner_bound_on_arr_rate(r, t, mu, w_sys)
  log(WARNING, "w_sys= {}, t= {}, r= {}, k= {}, mu= {}, arr_rate_ub={}".format(w_sys, t, r, k, mu, arr_rate_ub) )
      
  arr_rate_l = []
  sim_mds_n_k_E_T_l, sim_mds_n_1_E_T_l = [], []
  sim_hetero_simplex_E_T_l, sim_hetero2_simplex_E_T_l, sim_hetero3_simplex_E_T_l = [], [], []
  fj_2_2_E_T_l, sim_fj_2_E_T_l = [], []
  
  simplex_sm_E_T_l, sim_simplex_E_T_l, simplex_E_T_l, simplex_alt_E_T_l, E_T_simplex_matrix_analytic_l = [], [], [], [], []
  E_T_simplex_ub_l, E_T_simplex_lb_l, E_T_simplex_naive_lb_l, E_T_simplex_varki_gauri_lb_l = [], [], [], []
  simplex_trial_E_T_l, simplex_trial2_E_T_l, simplex_trial3_E_T_l, simplex_trial4_E_T_l = [], [], [], []
  
  simplex_wo_sys_sm_E_T_l = []
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7):
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
  # for arr_rate in numpy.arange(0.05, 1.26, 0.2):
    arr_rate_l.append(arr_rate)
    # sim
    num_f_run = 2
    # gamma=mu= 1, for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7)
    # if w_sys and t == 1:
    #   pass deneme
    if w_sys and t == 1:
      # sim_simplex_E_T_l= [
      #   0.6906089199666947,
      #   0.7708763120409886,
      #   0.9118019138920423,
      #   1.0672253191648233,
      #   1.3378310917617424,
      #   1.7793492983131396,
      #   2.8576296831667237]
      #   # 2.7576296831667237]
      sim_simplex_E_T_l= [
        0.690182518823904,
        0.7157334332190215,
        0.7564590917460736,
        0.7967851721092581,
        0.8349148470704952,
        0.9039070102346713,
        0.9633550609535413,
        1.008759189066751,
        1.093635647677166,
        1.1765798382640023,
        1.2715006428459072,
        1.4319781498055488,
        1.5838561388495074,
        1.769614490873982,
        2.026993960656419,
        2.4789818006986994,
        3.187262080396768,
        4.005550778820975,
        7.333009029591799,
        14.343792390544774]
    elif w_sys and t == 2:
      sim_simplex_E_T_l= [
        0.5538279298569435,
        0.6165688809324565,
        0.6982546908320676,
        0.8219708418135749,
        0.9997765522284852,
        1.3121746276083508,
        1.9956790864269085]
        # 1.9856790864269085]
    elif w_sys and t == 3:
      sim_simplex_E_T_l= [
        0.4692185744911441,
        0.5214332293554798,
       	0.5905958611717899,
      	0.6859006122530505,
      	0.8446347683684426,
      	1.0891266908697737,
      	1.6007724635530007]
    elif w_sys and t == 4:
      sim_simplex_E_T_l= [
        0.41308720307512653,
        0.46198842372979687,
        0.5172919002827522,
        0.5956587704334568,
        0.7313172804953793,
        0.946993937225123,
        1.4346171015038325]
    elif w_sys and t == 5:
      sim_simplex_E_T_l= [
        0.3773439036595755,
        0.4157220154779594,
        0.46831693074677627,
        0.5443438205065009,
        0.6560641337080825,
        0.8544174638869552,
        1.312970149903467]
    else:
      sim_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, w_sys=w_sys)
      sim_simplex_E_T_l.append(sim_simplex_E_T)
    
    marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
    def plot_split_to_all():
      # C = num_q*mu
      # def qmu_l(c):
      #   mu_ = C/2/(c+1)
      #   gamma = c*2*mu_
      #   return [gamma, mu_, mu_]
      # hetero_simplex_c = 0.001
      # sim_hetero_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(hetero_simplex_c) )
      # sim_hetero_simplex_E_T_l.append(sim_hetero_simplex_E_T)
      # hetero2_simplex_c = 0.3
      # sim_hetero2_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(hetero2_simplex_c) )
      # sim_hetero2_simplex_E_T_l.append(sim_hetero2_simplex_E_T)
      # hetero3_simplex_c = 0.6
      # sim_hetero3_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(hetero3_simplex_c) )
      # sim_hetero3_simplex_E_T_l.append(sim_hetero3_simplex_E_T)
      
      # fj_2_2_E_T_l.append(fj_2_2_E_T(arr_rate, C/2) )
      # fj_2_2_E_T_l.append(fj_2_2_E_T(arr_rate, mu) )
      # sim_fj_2_E_T_l.append(test_mds_n_k(num_f_run, arr_rate, C/2, n=2, k=2))
      
      # sim_mds_n_1_E_T_l.append(test_mds_n_k(num_f_run, arr_rate, mu, num_q, 1) )
      
      simplex_sm_E_T_l.append(simplex_sm_E_T(t, arr_rate, mu) )
      # simplex_wo_sys_sm_E_T_l.append(simplex_wo_sys_sm_E_T(t, arr_rate, mu) )
      
      if t == 1:
        simplex_E_T_l.append(simplex_w_one_repair__E_T(arr_rate, mu) )
        E_T_simplex_matrix_analytic_l.append(simplex_w_one_repair__E_T_matrix_analytic(t, arr_rate, mu) )
      elif t == 2:
        if w_sys:
          simplex_alt_E_T_l.append(simplex_w_two_repair__E_T(arr_rate, mu, M=2) )
          simplex_E_T_l.append(simplex_w_two_repair__E_T(arr_rate, mu, M=5) )
        else:
          simplex_E_T_l.append(simplex_wo_sys_w_two_repair__E_T(arr_rate, mu) )
      E_T_simplex_ub_l.append(E_T_simplex_lb(t, arr_rate, gamma, mu, ub=True) )
      E_T_simplex_lb_l.append(E_T_simplex_lb(t, arr_rate, gamma, mu) )
      E_T_simplex_naive_lb_l.append(E_T_simplex_lb(t, arr_rate, gamma, mu, naive=True) )
      E_T_simplex_varki_gauri_lb_l.append(E_T_simplex_varki_gauri_lb(t, arr_rate, gamma, mu) )
      # simplex_trial_E_T_l.append(simplex_w_one_repair__E_T_trial(1, t, arr_rate, mu) )
      # simplex_trial2_E_T_l.append(simplex_w_one_repair__E_T_trial(1.2, t, arr_rate, mu) )
      # simplex_trial3_E_T_l.append(simplex_w_one_repair__E_T_trial(1.5, t, arr_rate, mu) )
      # simplex_trial4_E_T_l.append(simplex_w_one_repair__E_T_trial(1.8, t, arr_rate, mu) )
      # plot.plot(arr_rate_l, sim_mds_n_k_E_T_l, 'ro', label="MDS({},{})".format(num_q, k) )
      plot.plot(arr_rate_l, simplex_sm_E_T_l, 'r', label=r'$E[\hat{T}_{SM}]$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, E_T_simplex_ub_l, 'm', label=r'$E[\hat{T}(\mathbf{\rho})]$', marker=next(marker), linestyle='', mew=2)
      # plot.plot(arr_rate_l, simplex_wo_sys_sm_E_T_l, 'ro', label="simplex_wo_sys_sm_t_{}".format(t) )
      # plot.plot(arr_rate_l, E_T_simplex_matrix_analytic_l, 'y', label=r'$E[\hat{T}_{MA}]$', marker=next(marker), linestyle='', mew=2)
      log(WARNING, "sim_simplex_E_T_l= {}".format(pprint.pformat(sim_simplex_E_T_l) ) )
      # plot.plot(arr_rate_l, simplex_alt_E_T_l, 'b', label=r'$E[\hat{T}_{LB}], M=2$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, sim_simplex_E_T_l, 'k', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
      # plot.plot(arr_rate_l, simplex_E_T_l, 'g', label=r'$E[\hat{T}_{LB}], M=5$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, E_T_simplex_lb_l, 'g', label=r'$E[\hat{T}(\mathbf{\rho_{max}})]$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, E_T_simplex_naive_lb_l, 'b', label=r'$E[\hat{T}(\mathbf{1})]$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, E_T_simplex_varki_gauri_lb_l, 'c', label=r'$E[\hat{T}_{fast-serial}]$', marker=next(marker), linestyle='', mew=2)
      color = iter(cm.rainbow(numpy.linspace(0, 2, 4) ) )
      # plot.plot(arr_rate_l, sim_hetero_simplex_E_T_l, 'o', color=next(color), label="hetero_simplex_c_{}".format(hetero_simplex_c) )
      # plot.plot(arr_rate_l, sim_hetero2_simplex_E_T_l, 'o', color=next(color), label="hetero_simplex_c_{}".format(hetero2_simplex_c) )
      # plot.plot(arr_rate_l, sim_hetero3_simplex_E_T_l, 'o', color=next(color), label="hetero_simplex_c_{}".format(hetero3_simplex_c) )
      # plot.plot(arr_rate_l, sim_fj_2_E_T_l, 'o', color=next(color), label="fj_2")
      # plot.plot(arr_rate_l, fj_2_2_E_T_l, 'o', color=next(color), label="model_fj_2")
      # plot.plot(arr_rate_l, simplex_trial_E_T_l, 'bo', label="trial_Simplex(t:{})".format(t) )
      # color = iter(cm.rainbow(numpy.linspace(0, 1, 4) ) )
      # plot.plot(arr_rate_l, simplex_trial2_E_T_l, 'o', color=next(color), label="trial_simplex2_t_{}".format(t) )
      # plot.plot(arr_rate_l, simplex_trial3_E_T_l, 'o', color=next(color), label="trial_simplex3_t_{}".format(t) )
      # plot.plot(arr_rate_l, simplex_trial4_E_T_l, 'bo', color=next(color), label="trial_simplex4_t_{}".format(t) )
      # plot.plot(arr_rate_l, sim_mds_n_1_E_T_l, 'bo', label="MDS({},1)".format(num_q) )
  def plot_split_to_one():
    arr_rate_ub = 0.9*arr_rate_ub_simplex_split_to_one(t, mu)
    log(WARNING, "w_sys= {}, t= {}, mu= {}, arr_rate_ub={}".format(w_sys, t, r, k, mu, arr_rate_ub) )
    
    arr_rate_sto_l = []
    E_T_simplex_sto_l = []
    
    for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
      arr_rate_sto_l.append(arr_rate)
      E_T_simplex_sto_l.append(E_T_simplex_split_to_one(t, arr_rate, mu) )
    log(WARNING, "sim_simplex_E_T_l= {}".format(pprint.pformat(sim_simplex_E_T_l) ) )
    plot.plot(arr_rate_l, sim_simplex_E_T_l, 'k', label=r'$E[T]$, split-to-all', marker=next(marker), linestyle='', mew=2)
    plot.plot(arr_rate_sto_l, E_T_simplex_sto_l, 'b', label=r'$E[T]$, split-to-one', marker=next(marker), linestyle='', mew=2)
  # plot_split_to_all()
  plot_split_to_one()
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  # plot.title(r't= {}, r= {}, k= {}, $\gamma$= {}, $\mu$= {}'.format(t, r, k, gamma, mu) )
  plot.title(r't= {}, $\gamma$= {}, $\mu$= {}'.format(t, gamma, mu) )
  plot.savefig("plot_simplex__t_{}.png".format(t) )
  log(WARNING, "done; t= {}, r= {}, k= {}".format(t, r, k) )

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
    simplex_sm_E_T_l, sim_hetero_simplex_E_T_l = [], []
    hetero_simplex_E_T_l = []
    for c in numpy.linspace(0.25, 5, 10):
      c_l.append(c)
      # sim
      simplex_sm_E_T_l.append(simplex_sm_E_T(t, arr_rate, mu, c) )
      
      num_f_run = 3
      sim_hetero_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(c) )
      sim_hetero_simplex_E_T_l.append(sim_hetero_simplex_E_T)
      
      hetero_simplex_E_T_l.append(simplex_w_one_repair__E_T(arr_rate, mu, c) )
    plot.plot(c_l, simplex_sm_E_T_l, 'o', color=next(color), label="SM-Simplex $\lambda={0:0.2f}$".format(arr_rate) )
    plot.plot(c_l, sim_hetero_simplex_E_T_l, 'o', color=next(color), label=r'Simplex $\lambda={0:0.2f}$'.format(arr_rate) )
    plot.plot(c_l, hetero_simplex_E_T_l, 'o', color=next(color), label=r'LB-Simplex $\lambda={0:0.2f}$'.format(arr_rate) )
    
  plot.xlabel("c")
  plot.ylabel("E[T] (s)")
  plot.title(r'Simplex(t:{}), heterogeneous servers; $c=\gamma/\mu$'.format(t) )
  plot.savefig("plot_simplex_w_varying_serv_rate_alloc__t_{}.png".format(t) )
  log(WARNING, "done; n= {} k= {} t= {}".format(num_q, k, t) )

def plot_avq():
  w_sys = True
  simplex_t, simplex_r, simplex_k = 3, 2, 2
  t, r, k = 1, simplex_t*simplex_r, simplex_k
  mu = 1.0
  gamma = mu
  arr_rate_ub = simplex_inner_bound_on_arr_rate(r, t, mu)
  log(WARNING, "w_sys= {}, t= {}, r= {}, k= {}, mu= {}, arr_rate_ub={}".format(w_sys, t, r, k, mu, arr_rate_ub) )
  
  arr_rate_l = []
  sim_E_T_avq_l, E_T_avq_l = [], []
  sim_E_T_simplex_l, E_T_simplex_sm_l= [], []
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/10):
  for arr_rate in numpy.arange(0.05, 1.26, 0.2):
  # for arr_rate in numpy.arange(0.05, 1.0, 0.2):
    arr_rate_l.append(arr_rate)
    # sim
    num_f_run = 1
    # sim_E_T_avq_l.append(test_simplex_q(num_f_run, arr_rate, mu, k, r, t, w_sys=w_sys) )
    sim_E_T_simplex_l.append(test_simplex_q(num_f_run, arr_rate, mu, simplex_k, simplex_r, simplex_t, w_sys=w_sys) )
    
    E_T_simplex_sm_l.append(simplex_sm_E_T(simplex_t, arr_rate, mu) )
    E_T_avq_l.append(E_T_avq_sys__mds_r_2(arr_rate, gamma, mu, r) )
  
  plot.plot(arr_rate_l, E_T_simplex_sm_l, 'ro', label="SM-Simplex(t:{},r:{},k:{})".format(simplex_t, simplex_r, simplex_k) )
  plot.plot(arr_rate_l, sim_E_T_simplex_l, 'ko', label="Simplex(t:{},r:{},k:{})".format(simplex_t, simplex_r, simplex_k) )
  plot.plot(arr_rate_l, E_T_avq_l, 'go', label="UB-AVQ(t:{},r:{},k:{})".format(t, r, k) )
  # plot.plot(arr_rate_l, sim_E_T_avq_l, 'ko', label="avq [t:{},r:{},k:{}]".format(t, r, k) )
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r't= {}, r= {}, k= {}, $\gamma$= {}, $\mu$= {}'.format(t, r, k, gamma, mu) )
  plot.savefig("plot_avq__t_{}_r_{}.png".format(t, r) )
  log(WARNING, "done; t= {}, r= {}, k= {}".format(t, r, k) )

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
  # plot_winning_freqs()
  plot_simplex(num_q)
  # plot_avq()
  # plot_simplex_w_varying_serv_rate_alloc(num_q)
