import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
from cycler import cycler
from random import expovariate
import sys, pprint, math, numpy, simpy, getopt, itertools

from simplex_sim_components import *
from simplex_models import *
# from mds_exp import *
# plot_color_l = ["indigo", "darkorange", "yellowgreen", "cyan", "darkmagenta", "darkred", "black", "slateblue", "goldenrod", "darksalmon", "forestgreen", "saddlebrown", "grey"]

def test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l=[], w_sys=True, mixed_traff=False):
  sim_E_T_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "arr_rate= {}, mu= {}, k= {}, r= {}, t= {}, qmu_l= {}, w_sys= {}, mixed_traff= {}". \
      format(arr_rate, mu, k, r, t, pprint.pformat(qmu_l), w_sys, mixed_traff) )
    
    num_q = int(1 + t*r) if w_sys else int(t*r)
    qid_l = [i for i in range(num_q) ]
    if len(qmu_l) == 0:
      qmu_l = [mu for i in range(num_q) ]
    
    env = simpy.Environment()
    if mixed_traff:
      sym__rgroup_l_map = {}
      sym_l = None
      if t < 3 or t == 4:
        sym_l = ['a', 'b']
      if t == 3:
        sym_l = ['a', 'b', 'c']
      for sym in sym_l:
        rgroup_l = []
        if w_sys and t == 1:
          if sym == 'a':
            rgroup_l.append([0] )
            rgroup_l.append([1, 2] )
            # for g in range(1, t+1):
            #   rgroup_l.append([2*g-1, 2*g] )
          elif sym == 'b':
            rgroup_l.append([1] )
            rgroup_l.append([0, 2] )
        elif w_sys and t == 2:
          if sym == 'a':
            rgroup_l.append([0] )
            rgroup_l.append([1, 2] )
            rgroup_l.append([3, 4] )
          elif sym == 'b':
            rgroup_l.append([1] )
            rgroup_l.append([0, 2] )
            rgroup_l.append([3, 4] )
        elif w_sys and t == 3:
          if sym == 'a':
            rgroup_l.append([0] )
            rgroup_l.append([1, 2] )
            rgroup_l.append([3, 4] )
            rgroup_l.append([5, 6] )
          elif sym == 'b':
            rgroup_l.append([1] )
            rgroup_l.append([0, 2] )
            rgroup_l.append([3, 5] )
            rgroup_l.append([4, 6] )
          elif sym == 'c':
            rgroup_l.append([3] )
            rgroup_l.append([0, 4] )
            rgroup_l.append([1, 5] )
            rgroup_l.append([2, 6] )
        elif w_sys and t == 4:
          if sym == 'a':
            rgroup_l.append([0] )
            rgroup_l.append([1, 2] )
            rgroup_l.append([3, 4] )
            rgroup_l.append([5, 6] )
            rgroup_l.append([7, 8] )
          elif sym == 'b':
            rgroup_l.append([1] )
            rgroup_l.append([0, 2] )
            rgroup_l.append([3, 4] )
            rgroup_l.append([5, 6] )
            rgroup_l.append([7, 8] )
        sym__rgroup_l_map[sym] = rgroup_l
      pg = MT_PacketGenerator(env, _id="p_gen",
                             adist=lambda: random.expovariate(arr_rate),
                             sdist=lambda: 1,
                             sym_l=sym_l)
        
      log(WARNING, "sym__rgroup_l_map=\n {}".format(pprint.pformat(sym__rgroup_l_map) ) )
      a_q = MT_AvQ("cds_q", env, qid_l, qmu_l, sym__rgroup_l_map)
      # aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 0.1)
      # a_q.join_q.out_m = aq_monitor
      pg.out = a_q
    else:
      pg = PacketGenerator(env, _id="p_gen",
                           adist=lambda: random.expovariate(arr_rate),
                           sdist=lambda: 1)
      a_q = AVQ("a_q", env, k, r, t, qid_l, qmu_l, w_sys=w_sys)
      # aq_monitor = AVQMonitor(env, aq=a_q, poll_dist=lambda: 0.1)
      # a_q.join_q.out_m = aq_monitor
      pg.out = a_q
      pg.init()
    env.run(until=50000) # env.run(until=50000) # env.run(until=500)
    if mixed_traff:
      print("pg.sym__n_sent= {}".format(pprint.pformat(pg.sym__n_sent) ) )
    
    
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
  print("sim_E_T_f_sum/num_f_run= {}".format(sim_E_T_f_sum/num_f_run) )
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
    pg.init()
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

def plot_simplex(num_q):
  w_sys = True # False
  t, r, k = 1, 2, 2
  mu = 1
  gamma = mu
  arr_rate_ub = simplex_inner_bound_on_arr_rate(r, t, mu, w_sys)
  log(WARNING, "w_sys= {}, t= {}, r= {}, k= {}, mu= {}, arr_rate_ub={}".format(w_sys, t, r, k, mu, arr_rate_ub) )
  
  arr_rate_l = []
  E_T_mds_n_k_sim_l, E_T_sim_mds_n_1_l = [], []
  sim_hetero_E_T_simplex_l, sim_hetero2_E_T_simplex_l, sim_hetero3_E_T_simplex_l = [], [], []
  E_T_fj_2_l, E_T_sim_fj_2_l = [], []
  
  E_T_simplex_sm_l, E_T_sim_simplex_l, E_T_simplex_l, E_T_simplex_alt_l, E_T_simplex_matrix_analytic_l = [], [], [], [], []
  E_T_simplex_approx_l, E_T_simplex_lb_l, E_T_simplex_naive_lb_l, E_T_simplex_varki_gauri_lb_l = [], [], [], []
  simplex_trial_E_T_l, simplex_trial2_E_T_l, simplex_trial3_E_T_l, simplex_trial4_E_T_l = [], [], [], []
  
  E_T_simplex_wo_sys_sm_l = []
  E_T_sim_simplex_mixed_traff_l = []
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7):
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
  # for arr_rate in numpy.arange(0.05, 1.26, 0.2):
    arr_rate_l.append(arr_rate)
    # sim
    num_f_run = 1
    # gamma=mu= 1, for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7)
    # if w_sys and t == 1:
    #   pass deneme
    if w_sys and t == 1:
      # E_T_sim_simplex_l= [
      #   0.6906089199666947,
      #   0.7708763120409886,
      #   0.9118019138920423,
      #   1.0672253191648233,
      #   1.3378310917617424,
      #   1.7793492983131396,
      #   2.8576296831667237]
      #   # 2.7576296831667237]
      E_T_sim_simplex_l= [
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
      E_T_sim_simplex_l= [
        0.5538279298569435,
        0.6165688809324565,
        0.6982546908320676,
        0.8219708418135749,
        0.9997765522284852,
        1.3121746276083508,
        1.9956790864269085]
        # 1.9856790864269085]
    elif w_sys and t == 3:
      E_T_sim_simplex_l= [
        0.4692185744911441,
        0.5214332293554798,
       	0.5905958611717899,
      	0.6859006122530505,
      	0.8446347683684426,
      	1.0891266908697737,
      	1.6007724635530007]
    elif w_sys and t == 4:
      E_T_sim_simplex_l= [
        0.41308720307512653,
        0.46198842372979687,
        0.5172919002827522,
        0.5956587704334568,
        0.7313172804953793,
        0.946993937225123,
        1.4346171015038325]
    elif w_sys and t == 5:
      E_T_sim_simplex_l= [
        0.3773439036595755,
        0.4157220154779594,
        0.46831693074677627,
        0.5443438205065009,
        0.6560641337080825,
        0.8544174638869552,
        1.312970149903467]
    else:
      E_T_sim_simplex_l.append(test_simplex_q(num_f_run, arr_rate, mu, k, r, t, w_sys=w_sys) )
    # mixed_traff
    if w_sys and t == 1:
      E_T_sim_simplex_mixed_traff_l= [
        0.6907408144088172,
        0.7534786674301083,
        0.8670549437839371,
        1.0157725875780679,
        1.2238618265783803,
        1.5909891765709532,
        2.286087812745277]
    elif w_sys and t == 2:
      E_T_sim_simplex_mixed_traff_l= [
        0.532665362364591,
        0.6086401739792854,
        0.6731769751198722,
        0.7854682740033085,
        0.9373880168977933,
        1.2287967062609655,
        1.6920865129907166]
    elif w_sys and t == 3:
      E_T_sim_simplex_mixed_traff_l= [
        0.45739370981666344,
        0.5113558732199666,
        0.5709895671208652,
        0.6539680930186696,
        0.782540234323826,
        0.976091978810901,
        1.3603386072899841]
    else:
      E_T_sim_simplex_mixed_traff_l.append(test_simplex_q(num_f_run, arr_rate, mu, k, r, t, w_sys=w_sys, mixed_traff=True) )
    
    # C = num_q*mu
    # def qmu_l(c):
    #   mu_ = C/2/(c+1)
    #   gamma = c*2*mu_
    #   return [gamma, mu_, mu_]
    # hetero_simplex_c = 0.001
    # sim_hetero_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(hetero_simplex_c) )
    # sim_hetero_E_T_simplex_l.append(sim_hetero_simplex_E_T)
    # hetero2_simplex_c = 0.3
    # sim_hetero2_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(hetero2_simplex_c) )
    # sim_hetero2_E_T_simplex_l.append(sim_hetero2_simplex_E_T)
    # hetero3_simplex_c = 0.6
    # sim_hetero3_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(hetero3_simplex_c) )
    # sim_hetero3_E_T_simplex_l.append(sim_hetero3_simplex_E_T)
    
    # E_T_fj_2_l.append(fj_2_2_E_T(arr_rate, C/2) )
    # E_T_fj_2_l.append(fj_2_2_E_T(arr_rate, mu) )
    # E_T_sim_fj_2_l.append(test_mds_n_k(num_f_run, arr_rate, C/2, n=2, k=2))
    
    # E_T_sim_mds_n_1_l.append(test_mds_n_k(num_f_run, arr_rate, mu, num_q, 1) )
    
    E_T_simplex_sm_l.append(simplex_sm_E_T(t, arr_rate, mu) )
    # E_T_simplex_wo_sys_sm_l.append(simplex_wo_sys_sm_E_T(t, arr_rate, mu) )
    
    if t == 1:
      E_T_simplex_l.append(simplex_w_one_repair__E_T(arr_rate, mu) )
      E_T_simplex_matrix_analytic_l.append(simplex_w_one_repair__E_T_matrix_analytic(t, arr_rate, mu) )
    elif t == 2:
      if w_sys:
        E_T_simplex_alt_l.append(simplex_w_two_repair__E_T(arr_rate, mu, M=2) )
        E_T_simplex_l.append(simplex_w_two_repair__E_T(arr_rate, mu, M=5) )
      else:
        E_T_simplex_l.append(simplex_wo_sys_w_two_repair__E_T(arr_rate, mu) )
    E_T_simplex_approx_l.append(E_T_simplex_lb(t, arr_rate, gamma, mu, ub=True) )
    E_T_simplex_lb_l.append(E_T_simplex_lb(t, arr_rate, gamma, mu) )
    E_T_simplex_naive_lb_l.append(E_T_simplex_lb(t, arr_rate, gamma, mu, naive=True) )
    E_T_simplex_varki_gauri_lb_l.append(E_T_simplex_varki_gauri_lb(t, arr_rate, gamma, mu) )
    # simplex_trial_E_T_l.append(simplex_w_one_repair__E_T_trial(1, t, arr_rate, mu) )
    # simplex_trial2_E_T_l.append(simplex_w_one_repair__E_T_trial(1.2, t, arr_rate, mu) )
    # simplex_trial3_E_T_l.append(simplex_w_one_repair__E_T_trial(1.5, t, arr_rate, mu) )
    # simplex_trial4_E_T_l.append(simplex_w_one_repair__E_T_trial(1.8, t, arr_rate, mu) )
    
    marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
    def plot_split_to_all():
      # plot.plot(arr_rate_l, E_T_mds_n_k_sim_l, 'ro', label="MDS({},{})".format(num_q, k) )
      # plot.plot(arr_rate_l, E_T_simplex_sm_l, 'r', label=r'$E[\hat{T}_{SM}]$', marker=next(marker), linestyle='', mew=2)
      # plot.plot(arr_rate_l, E_T_simplex_approx_l, 'm', label=r'$E[\hat{T}(\mathbf{\hat{\rho}})]$', marker=next(marker), linestyle='', mew=2)
      # plot.plot(arr_rate_l, E_T_simplex_wo_sys_sm_l, 'ro', label="simplex_wo_sys_sm_t_{}".format(t) )
      # plot.plot(arr_rate_l, E_T_simplex_matrix_analytic_l, 'y', label=r'$E[\hat{T}_{MA}]$', marker=next(marker), linestyle='', mew=2)
      # log(WARNING, "E_T_sim_simplex_l= {}".format(pprint.pformat(E_T_sim_simplex_l) ) )
      # plot.plot(arr_rate_l, E_T_sim_simplex_l, 'k', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, E_T_sim_simplex_l, 'k', label=r'$E[T], fixed-arrival$', marker=next(marker), linestyle='', mew=2)
      # log(WARNING, "E_T_sim_simplex_mixed_traff_l= {}".format(pprint.pformat(E_T_sim_simplex_mixed_traff_l) ) )
      plot.plot(arr_rate_l, E_T_sim_simplex_mixed_traff_l, color='brown', label=r'$E[T], mixed-arrival$', marker=next(marker), linestyle='', mew=2)
      # plot.plot(arr_rate_l, E_T_simplex_alt_l, 'b', label=r'$E[\hat{T}], M=2$', marker=next(marker), linestyle='', mew=2)
      # plot.plot(arr_rate_l, E_T_simplex_l, 'g', label=r'$E[\hat{T}], M=5$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, E_T_simplex_lb_l, 'g', label=r'$E[\hat{T}(\hat{\rho})]$', marker=next(marker), linestyle='', mew=2)
      plot.plot(arr_rate_l, E_T_simplex_naive_lb_l, 'b', label=r'$E[\hat{T}(1)]$', marker=next(marker), linestyle='', mew=2)
      # plot.plot(arr_rate_l, E_T_simplex_varki_gauri_lb_l, 'c', label=r'$E[\hat{T}_{fast-serial}]$', marker=next(marker), linestyle='', mew=2)
      color = iter(cm.rainbow(numpy.linspace(0, 2, 4) ) )
      # plot.plot(arr_rate_l, sim_hetero_E_T_simplex_l, 'o', color=next(color), label="hetero_simplex_c_{}".format(hetero_simplex_c) )
      # plot.plot(arr_rate_l, sim_hetero2_E_T_simplex_l, 'o', color=next(color), label="hetero_simplex_c_{}".format(hetero2_simplex_c) )
      # plot.plot(arr_rate_l, sim_hetero3_E_T_simplex_l, 'o', color=next(color), label="hetero_simplex_c_{}".format(hetero3_simplex_c) )
      # plot.plot(arr_rate_l, E_T_sim_fj_2_l, 'o', color=next(color), label="fj_2")
      # plot.plot(arr_rate_l, E_T_fj_2_l, 'o', color=next(color), label="model_fj_2")
      # plot.plot(arr_rate_l, simplex_trial_E_T_l, 'bo', label="trial_Simplex(t:{})".format(t) )
      # color = iter(cm.rainbow(numpy.linspace(0, 1, 4) ) )
      # plot.plot(arr_rate_l, simplex_trial2_E_T_l, 'o', color=next(color), label="trial_simplex2_t_{}".format(t) )
      # plot.plot(arr_rate_l, simplex_trial3_E_T_l, 'o', color=next(color), label="trial_simplex3_t_{}".format(t) )
      # plot.plot(arr_rate_l, simplex_trial4_E_T_l, 'bo', color=next(color), label="trial_simplex4_t_{}".format(t) )
      # plot.plot(arr_rate_l, E_T_sim_mds_n_1_l, 'bo', label="MDS({},1)".format(num_q) )
  def plot_split_to_one():
    arr_rate_ub = 0.9*arr_rate_ub_simplex_split_to_one(t, mu)
    log(WARNING, "w_sys= {}, t= {}, mu= {}, arr_rate_ub={}".format(w_sys, t, r, k, mu, arr_rate_ub) )
    
    arr_rate_sto_l = []
    E_T_simplex_sto_l = []
    
    for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
      arr_rate_sto_l.append(arr_rate)
      E_T_simplex_sto_l.append(E_T_simplex_split_to_one(t, arr_rate, mu) )
    log(WARNING, "E_T_sim_simplex_l= {}".format(pprint.pformat(E_T_sim_simplex_l) ) )
    plot.plot(arr_rate_l, E_T_sim_simplex_l, 'k', label=r'$E[T]$, replicate-to-all', marker=next(marker), linestyle='', mew=2)
    plot.plot(arr_rate_sto_l, E_T_simplex_sto_l, 'b', label=r'$E[T]$, select-one', marker=next(marker), linestyle='', mew=2)
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
    E_T_simplex_sm_l, sim_hetero_E_T_simplex_l = [], []
    hetero_E_T_simplex_l = []
    for c in numpy.linspace(0.25, 5, 10):
      c_l.append(c)
      # sim
      E_T_simplex_sm_l.append(simplex_sm_E_T(t, arr_rate, mu, c) )
      
      num_f_run = 3
      sim_hetero_simplex_E_T = test_simplex_q(num_f_run, arr_rate, mu, k, r, t, qmu_l(c) )
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
  
  # test_simplex_q(num_f_run=1, arr_rate=0.1, mu=1, k=2, r=2, t=1, mixed_traff=True)
  
  # plot_winning_freqs()
  # plot_simplex(num_q)
  plot_mds(num_q)
  # plot_mds_n_2(num_q)
  # plot_avq()
  # plot_simplex_w_varying_serv_rate_alloc(num_q)
