import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
from cycler import cycler
from random import expovariate
import sys, pprint, math, numpy, simpy, getopt, itertools

from mds_sim_components import *
from mds_models import *

def test_mds_n_k(num_f_run, mds_arr_rate, mu, n, k, r=None, sys_arr_rate=None, preempt=False):
  sim_E_T_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "mds_arr_rate= {}, mu= {}, n= {}, k= {}, r= {}, sys_arr_rate= {}".format(mds_arr_rate, mu, n, k, r, sys_arr_rate) )
    env = simpy.Environment()
    
    sys_arr_dist = None if sys_arr_rate is None else lambda: random.expovariate(sys_arr_rate)
    pg = MDS_PacketGenerator(env, _id="p_gen",
                             mds_arr_dist=lambda: random.expovariate(mds_arr_rate),
                             size_dist=lambda: 1, sys_arr_dist=sys_arr_dist)
    qid_l = ["{}".format(i) for i in range(n) ]
    # qserv_dist_l = [lambda: random.expovariate(mu) for i in range(n) ]
    mdsq = MDSQ("mdsq", env, k, qid_l, qserv_rate_l=[mu for i in range(n) ], r=r, preempt=preempt)
    mdsq_monitor = MDSQMonitor(env, q=mdsq, poll_dist=lambda: 1)
    pg.out = mdsq
    env.run(until=10) # env.run(until=50000) # env.run(until=5000)
    
    st_l = mdsq.join_sink.st_l
    if len(st_l) > 0:
      sim_E_T_f_sum += float(sum(st_l) )/len(st_l)
      # continue
    print("\n")
    print("mdsq.join_sink.qid__num_win_map= {}".format(pprint.pformat(mdsq.join_sink.qid__num_win_map) ) )
    total_num_wins = sum([n for i, n in mdsq.join_sink.qid__num_win_map.items() ] )
    print("total_num_wins= {}, \n pg= {}".format(total_num_wins, pg) )
    qid__win_freq_map = {i:float(n)/total_num_wins for i, n in mdsq.join_sink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
    
    # print("\n")
    # print("mdsq_monitor.state__counter_map= {}".format(pprint.pformat(mdsq_monitor.state__counter_map) ) )
    # total_counter = sum([c for rs, c in mdsq_monitor.state__counter_map.items() ] )
    # polled_state__counter_map = {rs:float(c)/total_counter for rs, c in mdsq_monitor.state__counter_map.items() }
    # print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
  return sim_E_T_f_sum/num_f_run

def plot_mds(num_q):
  n = 5 # num_q
  k = 4
  mu = 1.0
  arr_rate_ub = mds_exact_bound_on_arr_rate(mu, n, k) # mds_inner_bound_on_arr_rate
  log(WARNING, "n= {}, k= {}, mu= {}, arr_rate_ub={}".format(n, k, mu, arr_rate_ub) )
  
  arr_rate_l, p_lambda_l = [], []
  sim_fj_k_k_E_T_l = []
  varki_E_T_mds_n_k_l, E_T_mds_n_k_l, adj_E_T_mds_n_k_l, E_T_mds_n_k_sim_l = [], [], [], []
  E_T_mds_n_k_sm_l, adj_E_T_mds_n_k_sm_l, adj_2_E_T_mds_n_k_sm_l = [], [], []
  recur_E_T_mds_n_k_sm_l = []
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, 0.1):
  # for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/20):
  for arr_rate in numpy.arange(arr_rate_ub-0.05, arr_rate_ub+0.05, 0.05):
    arr_rate_l.append(arr_rate)
    
    E_T_mds_n_k_sm_l.append(E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    # adj_E_T_mds_n_k_sm_l.append(adj_E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    # adj_2_E_T_mds_n_k_sm_l.append(adj_2_E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    # recur_E_T_mds_n_k_sm_l.append(recur_E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    
    E_T_mds_n_k_l.append(E_T_mds_n_2(arr_rate, mu, n) )
    adj_E_T_mds_n_k_l.append(E_T_mds_n_2_adjustable(arr_rate, mu, n) )
    
    ro = float(arr_rate/mu)
    varki_mds_n_k_E_T = 1/mu * (harmonic_sum(n) - harmonic_sum(n-k) ) + \
      1/mu * ro*(gen_harmonic_sum(n, ro) - gen_harmonic_sum(n-k, ro) )
    varki_E_T_mds_n_k_l.append(varki_mds_n_k_E_T)
    # sim
    num_f_run = 1
    sim_mds_n_k_E_T = test_mds_n_k(num_f_run, arr_rate, mu, n, k)
    E_T_mds_n_k_sim_l.append(sim_mds_n_k_E_T)
    
    sim_fj_k_k_E_T = 0 #  test_mds_n_k(num_f_run, arr_rate, mu, k, k)
    sim_fj_k_k_E_T_l.append(sim_fj_k_k_E_T)
    
    # sim_mds_n_1_E_T = 0 # test_mds_n_k(num_f_run, arr_rate, mu, n, 1)
    # E_T_sim_mds_n_1_l.append(sim_mds_n_1_E_T)
  plot.plot(arr_rate_l, E_T_mds_n_k_sm_l, 'ro', label="SM-MDS({},{})".format(n, k) )
  # plot.plot(arr_rate_l, adj_E_T_mds_n_k_sm_l, 'bo', label="adj_sm_MDS({},{})".format(n, k) )
  # plot.plot(arr_rate_l, adj_2_E_T_mds_n_k_sm_l, 'co', label="adj_2_sm_MDS({},{})".format(n, k) )
  # plot.plot(arr_rate_l, recur_E_T_mds_n_k_sm_l, 'yo', label="recur_sm_MDS({},{})".format(n, k) )
  print("E_T_mds_n_k_sim_l= {}".format(E_T_mds_n_k_sim_l) )
  plot.plot(arr_rate_l, E_T_mds_n_k_sim_l, 'ko', label="MDS({},{})".format(n, k) )
  # plot.plot(arr_rate_l, adj_E_T_mds_n_k_l, 'go', label="LB-MDS({},{})".format(n, k) )
  plot.plot(arr_rate_l, varki_E_T_mds_n_k_l, 'mo', label="Varki-MDS({},{})".format(n, k) )
  
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, k= {}, identical servers $\mu$= {}'.format(n, k, mu) )
  plot.savefig("plot_mds__n_{}_k_{}.png".format(n, k) )
  log(WARNING, "done; n= {}, k= {}".format(n, k) )

def plot_mds_n_2(n):
  n = 3
  k = 2
  mu = 1
  sys_arr_rate = 0.5 # None
  # gamma = mu
  arr_rate_ub = mds_inner_bound_on_arr_rate(mu, n, k) # mds_exact_bound_on_arr_rate(mu, n, k)
  if sys_arr_rate is not None:
    arr_rate_ub = arr_rate_ub - sys_arr_rate
  log(WARNING, "n= {}, k= {}, mu= {}, sys_arr_rate= {}, arr_rate_ub={}".format(n, k, mu, sys_arr_rate, arr_rate_ub) )
  
  arr_rate_l = []
  E_T_mds_n_k_sm_l, E_T_mds_n_k_sim_l, E_T_mds_n_k_lb_l, E_T_mds_n_k_varki_gauri_lb_l = [], [], [], []
  # for arr_rate in [*numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7), arr_rate_ub-0.1]:
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7):
  # for arr_rate in numpy.arange(arr_rate_ub-0.15, arr_rate_ub+0.05, 0.05):
    arr_rate_l.append(arr_rate)
    # sim
    num_f_run = 1
    if sys_arr_rate is None and n == 3:
      E_T_mds_n_k_sim_l= [
        0.8556037128377268,
        0.9607942308333743,
        1.0536232157615888,
        1.2124623964086967,
        1.4178943036232037,
        1.7094303177886272,
        2.2350543957767584]
    elif sys_arr_rate is None and n == 4:
      E_T_mds_n_k_sim_l= [
        0.5812755345723262,
        0.6569461385017236,
        0.7425261002623366,
        0.8560773825930664,
        1.032387913703661,
        1.294593549223201,
        1.8578248960696804]
    elif sys_arr_rate is None and n == 5:
      E_T_mds_n_k_sim_l= [
        0.46425896730877336,
        0.5146692509223326,
        0.5789923822694713,
        0.6703636918860807,
        0.8008080993505126,
        1.0426485991760615,
        1.5271327221577]
    else:
      E_T_mds_n_k_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu, n, k, sys_arr_rate=sys_arr_rate) )
    
    E_T_mds_n_k_sm_l.append(E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    if k == 2:
      E_T_mds_n_k_lb_l.append(E_T_mds_n_2(arr_rate, mu, n) )
    E_T_mds_n_k_varki_gauri_lb_l.append(E_T_mds_n_k_varki_gauri_lb(arr_rate, mu, n, k) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(arr_rate_l, E_T_mds_n_k_sm_l, color='red', label=r'$E[\hat{T}_{SM}]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_mds_n_k_sim_l= {}".format(pprint.pformat(E_T_mds_n_k_sim_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_k_sim_l, color='black', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_k_lb_l, color='green', label=r'$E[\hat{T}_{LB}]$', marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_k_varki_gauri_lb_l, color='blue', label=r'$E[\hat{T}_{Gauri}]$', marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  # plot.title(r't= {}, r= {}, k= {}, $\gamma$= {}, $\mu$= {}'.format(t, r, k, gamma, mu) )
  plot.title(r'n= {}, k= {}, $\mu$= {}'.format(n, k, mu) )
  plot.savefig("plot_mds_{}_2.png".format(n) )
  log(WARNING, "done; n= {}, sys_arr_rate= {}".format(n, sys_arr_rate) )

def plot_mds_n_r_2(n):
  n = 4
  r = 3
  k = 2
  mu = 1
  arr_rate_ub = n/r * mds_exact_bound_on_arr_rate(mu, r, k)
  log(WARNING, "n= {}, r= {}, k= {}, mu= {}, arr_rate_ub={}".format(n, r, k, mu, arr_rate_ub) )
  
  arr_rate_l = []
  E_T_mds_n_r_k_sim_l, E_T_mds_n_r_k_lb_l = [], []
  E_T_mds_n_k_sm_l, E_T_mds_n_k_sim_l, E_T_mds_n_k_lb_l, E_T_mds_n_k_varki_gauri_lb_l = [], [], [], []
  
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7):
    arr_rate_l.append(arr_rate)
    # sim
    num_f_run = 1
    if n == 5:
      E_T_mds_n_k_sim_l= [
        0.46277525699566846,
        0.5129005496906122,
        0.5988096492061943,
        0.7127299561917352,
        0.8968951989714525,
        1.3006600764338618,
        2.5795503865912557]
    elif n == 10:
      E_T_mds_n_k_sim_l= [
        0.2077709453034781,
        0.23974071436935418,
        0.2735596331276213,
        0.33250455607082485,
        0.42319272770694766,
        0.5998126329293456,
        1.2145654562282062]
    else:
      E_T_mds_n_k_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu, n, k) )
    #
    sim = False
    if r == 3:
      if n == 5:
        E_T_mds_n_r_k_sim_l= [
          0.8481985778068147,
          0.9504370765921353,
          1.0894126844007486,
          1.3057841854966015,
          1.5967802446041406,
          2.2072872739194938,
          3.7146063520660775]
      elif n == 10:
        E_T_mds_n_r_k_sim_l= [
          0.8428394336727548,
          0.9436568150475411,
          1.0617479941774777,
          1.2454907412899359,
          1.5196422689122135,
          1.9965888962169007,
          3.0864873347481456]
      else:
        sim = True
    else:
      sim = True
    if sim:
      E_T_mds_n_r_k_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu, n, k, r=r) )
    
    E_T_mds_n_k_sm_l.append(E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    if k == 2:
      E_T_mds_n_k_lb_l.append(E_T_mds_n_2(arr_rate, mu, n) )
    
    E_T_mds_n_r_k_lb_l.append(E_T_mds_n_r_2(arr_rate, mu, n, r) )
    E_T_mds_n_k_varki_gauri_lb_l.append(E_T_mds_n_k_varki_gauri_lb(arr_rate, mu, n, k) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  # plot.plot(arr_rate_l, E_T_mds_n_k_sm_l, color='red', label=r'$E[\hat{T}_{SM}]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_mds_n_k_sim_l= {}".format(pprint.pformat(E_T_mds_n_k_sim_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_k_sim_l, color='black', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_mds_n_r_k_sim_l= {}".format(pprint.pformat(E_T_mds_n_r_k_sim_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_r_k_sim_l, color='gray', label=r'$E[T], r: {}$'.format(r), marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_k_lb_l, color='green', label=r'$E[\hat{T}_{LB}]$', marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_r_k_lb_l, color='slateblue', label=r'$E[\hat{T}_{LB}], r$', marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_k_varki_gauri_lb_l, color='blue', label=r'$E[\hat{T}_{Gauri}]$', marker=next(marker), linestyle='', mew=2)
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, r= {}, k= {}, $\mu$= {}'.format(n, r, k, mu) )
  plot.savefig("plot_mds_{}_{}_{}.png".format(n, r, k) )
  log(WARNING, "done; n= {}, r= {}, k= {}".format(n, r, k) )

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
  
  test_mds_n_k(num_f_run=1, mds_arr_rate=0.5, mu=1, n=3, k=2, r=None, sys_arr_rate=0.5, preempt=True)
  
  # plot_mds(num_q)
  # plot_mds_n_2(num_q)
  # plot_mds_n_r_2(num_q)
  