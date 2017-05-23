import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
from cycler import cycler
from random import expovariate
import sys, pprint, math, numpy, simpy, getopt, itertools

from mds_sim_components import *
from mds_models import *

def test_m_m_1(num_f_run, arr_rate, mu, long_sim=False):
  E_T_sim_f_sum = 0
  if long_sim:
    num_f_run = 5
  for f in range(num_f_run):
    log(WARNING, "arr_rate= {}, mu= {}".format(arr_rate, mu) )
    env = simpy.Environment()
    pg = PacketGenerator(env, _id="pg",
                         arr_dist=lambda: random.expovariate(arr_rate) )
    m1_q = S1_Q(_id=0, env=env, serv_rate=mu)
    # qm = QMonitor(env, q=m1_q, dist=lambda: 0.5)
    pg.out = m1_q
    pg.init()
    env.run(until=10**(long_sim) * 50000)
    
    E_T_sim_f_sum += sum(m1_q.qt_l)/len(m1_q.qt_l)
  return E_T_sim_f_sum/num_f_run
  
  # print("received: {}, dropped {}, sent {}".format(m1_q.n_recved, m1_q.n_dropped, pg.n_sent) )
  # print("loss rate: {}".format(float(m1_q.n_dropped)/m1_q.n_recved) )
  
  # print("arr_rate= {}, serv_rate= {}".format(arr_rate, serv_rate) )
  # E_N = arr_rate/(serv_rate - arr_rate)
  # print("E[N]= {}, sim_E[N]= {:.3f}".format(E_N, float(sum(qm.n_l) )/len(qm.n_l) ) )
  # E_T = E_N/arr_rate
  # print("E[T]= {:.3f}, sim_E[T]= {:.3f}".format(E_T, float(sum(m1_q.qt_l) )/len(m1_q.qt_l) ) )
  # E_W = E_T - 1/serv_rate
  # print("E[W]= {:.3f}, sim_E[W]= {:.3f}".format(E_W, float(sum(m1_q.wt_l) )/len(m1_q.wt_l) ) )

def plot_m_m_1():
  mu = 1
  log(WARNING, "mu= {}".format(mu) )
  
  arr_rate_l = []
  E_T_m_m_1_sim_l, E_T_m_m_1_l = [], []
  
  num_f_run = 1
  for arr_rate in numpy.arange(0.05, mu, mu/20):
    arr_rate_l.append(arr_rate)
  
    E_T_m_m_1_sim_l.append(test_m_m_1(num_f_run, arr_rate, mu, long_sim=(arr_rate >= 0.8) ) )
    E_T_m_m_1_l.append(1/(mu-arr_rate) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  print("E_T_m_m_1_sim_l= {}".format(pprint.pformat(E_T_m_m_1_sim_l) ) )
  plot.plot(arr_rate_l, E_T_m_m_1_sim_l, color='black', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_m_m_1_l= {}".format(pprint.pformat(E_T_m_m_1_l) ) )
  plot.plot(arr_rate_l, E_T_m_m_1_l, color='green', label=r'$E[\hat{T}]$', marker=next(marker), linestyle='', mew=2)
  
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'$\mu$= {}'.format(mu) )
  plot.savefig("plot_m_m_1.png")
  log(WARNING, "done.")

def test_mds_n_k(num_f_run, arr_rate, mu, n, k, r=None, l_s=None, preempt=False, long_sim=False,
                 p_i_l=[] ):
  E_T_sim_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "arr_rate= {}, mu= {}, n= {}, k= {}, r= {}".format(arr_rate, mu, n, k, r) )
    env = simpy.Environment()
    # sys_arr_dist = None if l_s is None else lambda: random.expovariate(n*l_s)
    pg = MDS_PG(env, _id="pg",
                inter_arr_dist=lambda: random.expovariate(arr_rate) )
    mdsq = MDSQ("mdsq", env, k, qid_l=range(n), qmu_l=n*[mu], r=r, preempt=preempt)
    monitor = MDSQMonitor(env, q=mdsq, poll_dist=lambda: 1)
    pg.out = mdsq
    env.run(until=60**(long_sim) * 50000) # env.run(until=5000)
    
    st_l = mdsq.join_sink.st_l
    if len(st_l) > 0:
      E_T_sim_f_sum += float(sum(st_l) )/len(st_l)
    
    total_num_wins = sum([n for i, n in mdsq.join_sink.qid__num_win_map.items() ] )
    qid__win_freq_map = {i:float(n)/total_num_wins for i, n in mdsq.join_sink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
    total_n_types = sum(mdsq.type__num_m)
    p_i_l[:] = [n/total_n_types for t,n in enumerate(mdsq.type__num_m) ]
    print("p_i_l= {}".format(p_i_l) )
    
    # print("\n")
    # # print("monitor.state__counter_map= {}".format(pprint.pformat(monitor.state__counter_map) ) )
    # total_counter = sum([c for rs, c in monitor.state__counter_map.items() ] )
    # polled_state__freq_map = {rs:float(c)/total_counter for rs, c in monitor.state__counter_map.items() }
    # print("polled_state__freq_map= {}".format(pprint.pformat(polled_state__freq_map) ) )
    # print("----------------------------------------")
  E_T = E_T_sim_f_sum/num_f_run
  print(">> E_T= {}".format(E_T) )
  if E_T > 30: return None
  return E_T

def plot_mds():
  n, k = 4, 2
  mu = 1
  arr_rate_ub = mds_exact_bound_on_arr_rate(mu, n, k) - 0.05 # mds_inner_bound_on_arr_rate
  log(WARNING, "n= {}, k= {}, mu= {}, arr_rate_ub={}".format(n, k, mu, arr_rate_ub) )
  
  arr_rate_l = []
  E_T_mds_n_k_sim_l, E_T_mds_n_k_approx_l = [], []
  E_T_mds_n_k_sim_based_approx_l = []
  E_T_mds_n_k_sm_l, E_T_mds_n_k_varki_gauri_lb_l = [], []
  
  sim_mds_n_k = False
  if k == 3 and n == 44:
    E_T_mds_n_k_sim_l= [
      1.1289106189418618,
      1.2378560824181688,
      1.460233615664515,
      1.7249196226715797,
      2.204430925193168,
      3.279193514129549,
      3.6529696076728753,
      4.136404228840242,
      4.946227654840564,
      6.081240521624994,
      7.540744921118187,
      16.511235591950822]
  else:
    sim_mds_n_k = True
  
  num_f_run = 1
  for arr_rate in [*numpy.linspace(0.05, 0.8*arr_rate_ub, 5, endpoint=False), *numpy.linspace(0.8*arr_rate_ub, arr_rate_ub, 7) ]:
    arr_rate_l.append(arr_rate)
    
    E_T_mds_n_k_sm_l.append(E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    p_i_l = []
    if sim_mds_n_k:
      E_T_mds_n_k_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu, n, k, p_i_l=p_i_l) )
    E_T_mds_n_k_sim_based_approx_l.append(E_T_mds_n_k(arr_rate, mu, n, k, p_i_l=p_i_l) )
    E_T_mds_n_k_approx_l.append(E_T_mds_n_k(arr_rate, mu, n, k, incremental=True) )
    
    E_T_mds_n_k_varki_gauri_lb_l.append(E_T_mds_n_k_varki_gauri_lb(arr_rate, mu, n, k) )
  plot.plot(arr_rate_l, E_T_mds_n_k_sm_l, label=r'Split-merge upper-bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  print("E_T_mds_n_k_sim_l= {}".format(pprint.pformat(E_T_mds_n_k_sim_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_k_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_k_sim_based_approx_l, label=r'Sim-based approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_k_approx_l, label=r'Best approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_k_varki_gauri_lb_l, label=r'Fast-serial lower-bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  
  plot.legend()
  plot.xlabel(r'Arrival rate $\lambda$')
  plot.ylabel("Average download time E[T] (s)")
  plot.title(r'$n= {}, k= {}, \mu= {}$'.format(n, k, mu) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mds_n_{}_k_{}.pdf".format(n, k) )
  log(WARNING, "done; n= {}, k= {}".format(n, k) )

def plot_mds_n_2(n):
  n = 3
  k = 2
  mu = 1
  l_s = None # 0.1 # None
  arr_rate_ub = mds_exact_bound_on_arr_rate(mu, n, k) # 0.9* # mds_inner_bound_on_arr_rate mds_exact_bound_on_arr_rate(mu, n, k)
  long_sim_threshold = mds_inner_bound_on_arr_rate(mu, n, k)
  if l_s is not None:
    arr_rate_ub = mds_exact_bound_on_arr_rate(mu-l_s, n, k) # arr_rate_ub - l_s
    long_sim_threshold = 0.25
  log(WARNING, "n= {}, k= {}, mu= {}, l_s= {}, arr_rate_ub={}".format(n, k, mu, l_s, arr_rate_ub) )
  
  arr_rate_l = []
  E_T_mds_n_2_sm_l, E_T_mds_n_2_sim_l, E_T_mds_n_2_lb_l, E_T_mds_n_2_varki_gauri_lb_l = [], [], [], []
  E_T_mds_n_2_ref_sim_l = []
  E_T_fj_n_2_l, E_T_mds_n_2_no_task_cancel_lb_l = [], []
  # for arr_rate in [*numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7), arr_rate_ub-0.1]:
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/15):
  # for arr_rate in numpy.arange(arr_rate_ub-0.4, arr_rate_ub+0.1, 0.05):
  # for arr_rate in numpy.arange(1.4, arr_rate_ub, arr_rate_ub/10):
    arr_rate_l.append(arr_rate)
    # sim
    num_f_run = 1
    # arr_rate_ub = mds_inner_bound_on_arr_rate(mu, n, k)
    if l_s is None and n == 22:
      E_T_mds_n_2_sim_l= [
        1.5513603858625389,
        1.8384528867954009,
        2.1617336375623855,
        2.7540010811170306,
        3.7524068561252295,
        5.851931858663693,
        15.213311485849891]
    elif l_s is None and n == 33:
      E_T_mds_n_2_sim_l= [
        0.8556037128377268,
        0.9607942308333743,
        1.0536232157615888,
        1.2124623964086967,
        1.4178943036232037,
        1.7094303177886272,
        2.2350543957767584]
    elif l_s is None and n == 44:
      E_T_mds_n_2_sim_l= [
        0.5812755345723262,
        0.6569461385017236,
        0.7425261002623366,
        0.8560773825930664,
        1.032387913703661,
        1.294593549223201,
        1.8578248960696804]
    elif l_s is None and n == 55:
      E_T_mds_n_2_sim_l= [
        0.46425896730877336,
        0.5146692509223326,
        0.5789923822694713,
        0.6703636918860807,
        0.8008080993505126,
        1.0426485991760615,
        1.5271327221577]
    # arr_rate_ub = mds_exact_bound_on_arr_rate(mu-l_s, n, k)
    elif l_s == 0.5 and n == 2:
      E_T_mds_n_2_sim_l= [
        3.311723509041599,
        3.9064807930294005,
        4.7163872409392775,
        6.025729982696663,
        9.094495734021644,
        15.142935526573915,
        63.52972948214127]
    elif l_s == 0.5 and n == 3:
      E_T_mds_n_2_sim_l= [
        1.746180491658557,
        1.9809069915526545,
        2.3429080103826627,
        2.882124001524471,
        3.7103667759692947,
        5.595304302324676,
        15.40908652933348]
    elif l_s == 0.1 and n == 22:
      E_T_mds_n_2_sim_l= [
        1.7623854828635468,
        2.0332476336801144,
        2.438552758699048,
        3.14359003472221,
        4.231802664045872,
        6.691411227717583,
        16.71678713361064]
    elif l_s == 0.1 and n == 3:
      E_T_mds_n_2_sim_l= [
        0.9465711129456883,
        1.0855081905157227,
        1.2517454120269855,
        1.5100066194076383,
        1.9085846851322814,
        2.853524260299928,
        5.65779334183252]
    elif l_s == 0.8 and n == 2:
      E_T_mds_n_2_sim_l= [
        9.823807771880826,
        12.503365770563832,
        16.291943558996547,
        20.725176924660058,
        42.25979927352205,
        224.22674504841316]
    elif l_s == 0.8 and n == 3:
      E_T_mds_n_2_sim_l= [
        4.796265556856205,
        5.81356900419174,
        6.902232910641226,
        8.477814813825821,
        12.754118189436356,
        24.813176858705937]
    else:
      E_T_mds_n_2_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu, n, k, l_s=l_s, long_sim=(arr_rate >= long_sim_threshold) ) )

    #
    # if l_s == 0.11 and n == 3:
    #   pass
    # else:
    #   E_T_mds_n_2_ref_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu-l_s, n, k) )
    E_T_mds_n_2_sm_l.append(E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    E_T_mds_n_2_lb_l.append(E_T_mds_n_2(arr_rate, mu, n, l_s) )
    E_T_mds_n_2_varki_gauri_lb_l.append(E_T_mds_n_k_varki_gauri_lb(arr_rate, mu, n, k) )
    # E_T_fj_n_2_l.append(E_T_fj_n_2(arr_rate, mu, n, l_s) )
    # E_T_mds_n_2_no_task_cancel_lb_l.append(E_T_mds_n_2_no_task_cancel_lb(arr_rate, mu, n, l_s) )
  print("arr_rate_l= {}".format(pprint.pformat(arr_rate_l) ) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  plot.plot(arr_rate_l, E_T_mds_n_2_sm_l, color='red', label=r'$E[\hat{T}_{SM}]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_mds_n_2_sim_l= {}".format(pprint.pformat(E_T_mds_n_2_sim_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_2_sim_l, color='black', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_mds_n_2_lb_l= {}".format(pprint.pformat(E_T_mds_n_2_lb_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_2_lb_l, color='green', label=r'$E[\hat{T}_{LB}]$', marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_2_varki_gauri_lb_l, color='blue', label=r'$E[\hat{T}_{serial-model}]$', marker=next(marker), linestyle='', mew=2)
  # plot.plot(arr_rate_l, E_T_fj_n_2_l, color='slateblue', label=r'$E[\hat{T}_{FJ(n,2)}]$', marker=next(marker), linestyle='', mew=2)
  # plot.plot(arr_rate_l, E_T_mds_n_2_no_task_cancel_lb_l, color='palegreen', label=r'$E[\hat{T}_{LB}]$, no-task-cancel', marker=next(marker), linestyle='', mew=2)
  # plot.plot(arr_rate_l, E_T_mds_n_2_ref_sim_l, color='darkviolet', label=r'$E[T]$, mu={}'.format(mu-l_s), marker=next(marker), linestyle='', mew=2)
  
  plot.legend()
  plot.xlabel(r'$\lambda$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, k= {}, $\mu$= {}'.format(n, k, mu) )
  # plot.title(r'n= {}, k= {}, $\mu$= {}, $\lambda_s$= {}'.format(n, k, mu, l_s) )
  plot.savefig("plot_mds_{}_2.png".format(n) )
  log(WARNING, "done; n= {}, l_s= {}".format(n, l_s) )

def plot_mds_n_r_2(n):
  n = 4
  r = 3
  k = 2
  mu = 1
  arr_rate_ub = n/r * mds_exact_bound_on_arr_rate(mu, r, k)
  log(WARNING, "n= {}, r= {}, k= {}, mu= {}, arr_rate_ub={}".format(n, r, k, mu, arr_rate_ub) )
  
  arr_rate_l = []
  E_T_mds_n_r_2_sim_l, E_T_mds_n_r_2_lb_l = [], []
  E_T_mds_n_2_sm_l, E_T_mds_n_2_sim_l, E_T_mds_n_2_lb_l, E_T_mds_n_2_varki_gauri_lb_l = [], [], [], []
  
  for arr_rate in numpy.arange(0.05, arr_rate_ub, arr_rate_ub/7):
    arr_rate_l.append(arr_rate)
    # sim
    num_f_run = 1
    if n == 5:
      E_T_mds_n_2_sim_l= [
        0.46277525699566846,
        0.5129005496906122,
        0.5988096492061943,
        0.7127299561917352,
        0.8968951989714525,
        1.3006600764338618,
        2.5795503865912557]
    elif n == 10:
      E_T_mds_n_2_sim_l= [
        0.2077709453034781,
        0.23974071436935418,
        0.2735596331276213,
        0.33250455607082485,
        0.42319272770694766,
        0.5998126329293456,
        1.2145654562282062]
    else:
      E_T_mds_n_2_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu, n, k) )
    #
    sim = False
    if r == 3:
      if n == 5:
        E_T_mds_n_r_2_sim_l= [
          0.8481985778068147,
          0.9504370765921353,
          1.0894126844007486,
          1.3057841854966015,
          1.5967802446041406,
          2.2072872739194938,
          3.7146063520660775]
      elif n == 10:
        E_T_mds_n_r_2_sim_l= [
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
      E_T_mds_n_r_2_sim_l.append(test_mds_n_k(num_f_run, arr_rate, mu, n, k, r=r) )
    
    E_T_mds_n_2_sm_l.append(E_T_mds_n_k_sm(arr_rate, mu, n, k) )
    if k == 2:
      E_T_mds_n_2_lb_l.append(E_T_mds_n_2(arr_rate, mu, n) )
    
    E_T_mds_n_r_2_lb_l.append(E_T_mds_n_r_2(arr_rate, mu, n, r) )
    E_T_mds_n_2_varki_gauri_lb_l.append(E_T_mds_n_k_varki_gauri_lb(arr_rate, mu, n, k) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  # plot.plot(arr_rate_l, E_T_mds_n_2_sm_l, color='red', label=r'$E[\hat{T}_{SM}]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_mds_n_k_sim_l= {}".format(pprint.pformat(E_T_mds_n_k_sim_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_k_sim_l, color='black', label=r'$E[T]$', marker=next(marker), linestyle='', mew=2)
  print("E_T_mds_n_r_2_sim_l= {}".format(pprint.pformat(E_T_mds_n_r_2_sim_l) ) )
  plot.plot(arr_rate_l, E_T_mds_n_r_2_sim_l, color='gray', label=r'$E[T], r: {}$'.format(r), marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_2_lb_l, color='green', label=r'$E[\hat{T}_{LB}]$', marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_r_2_lb_l, color='slateblue', label=r'$E[\hat{T}_{LB}], r$', marker=next(marker), linestyle='', mew=2)
  plot.plot(arr_rate_l, E_T_mds_n_2_varki_gauri_lb_l, color='blue', label=r'$E[\hat{T}_{serial-model}]$', marker=next(marker), linestyle='', mew=2)
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
  # log(WARNING, "opt_map= {}".format(pprint.pformat(opt_map) ) )
  # num_q = int(opt_map["--num_q"] )
  
  # test_mds_n_k(num_f_run=1, arr_rate=0.5, mu=1, n=3, k=2, r=None, l_s=0.5, preempt=True)
  
  plot_mds()
  # plot_mds_n_2()
  # plot_m_m_1()
  # plot_mds_n_r_2()
  
