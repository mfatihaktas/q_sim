import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# matplotlib.rcParams['ps.useafm'] = True
# matplotlib.rcParams['pdf.use14corefonts'] = True
# matplotlib.rcParams['text.usetex'] = True
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
from random import expovariate
import sys, pprint, math, numpy, simpy, getopt, itertools

from simplex_sim_components import *
from simplex_models import *
from mds_models import mds_exact_bound_on_arr_rate
from mds_exp import test_mds_n_k

def plot_simplex_steady_state_prob_hist():
  k, r, t = 2, 2, 1
  num_q = int(1 + t*r)
  qid_l = ["{}".format(i) for i in range(1, num_q + 1) ]
  qmu_l = [1, 1, 1]
  def get_state_prob_map(arr_rate):
    log(WARNING, "arr_rate= {}, k= {}, r= {}, t= {}, qmu_l= {}".format(arr_rate, k, r, t, pprint.pformat(qmu_l) ) )
    env = simpy.Environment()
    pg = PG(env, _id="pg",
            inter_arr_dist=lambda: random.expovariate(arr_rate) )
    avq = AVQ("avq", env, k, r, t, qid_l, qmu_l=qmu_l)
    monitor = AVQMonitor(env, avq, poll_dist=lambda: 0.1)
    avq.join_q.out_m = monitor
    pg.out = avq
    env.run(until=50000)
    
    # print("monitor.polled_state__counter_map= {}".format(pprint.pformat(monitor.polled_state__counter_map) ) )
    total_counter = sum([c for rs, c in monitor.polled_state__counter_map.items() ] )
    state_prob_map = {rs:float(c)/total_counter for rs, c in monitor.polled_state__counter_map.items() }
    # print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
    return state_prob_map # ['0,(0,0)']
  # for arr_rate in numpy.arange(0.05, 1.2, 0.1):
  color = iter(cm.rainbow(numpy.linspace(0, 1, 20) ) )
  plot.figure(figsize=(20,10) )
  for arr_rate in numpy.arange(0.05, 1.3, 0.1):
  # for arr_rate in numpy.arange(0.05, 0.1, 0.1):
    state_prob_map = get_state_prob_map(arr_rate)
    
    def state(kp, i, j):
      return "{},({},{})".format(kp, i, j)
    i__tau_l_map = {}
    for i in range(10):
      if i not in i__tau_l_map:
        i__tau_l_map[i] = []
      for kp in range(i, 10):
        s_u, s_l = state(kp, i, 0), state(kp+1, i, 0)
        if s_u in state_prob_map and s_l in state_prob_map:
          i__tau_l_map[i].append(state_prob_map[s_l]/state_prob_map[s_u] )
        # if state(k+1, 0, i) in state_prob_map:
        #   i__tau_l_map[i].append(state_prob_map[state(k+1, 0, i) ] /state_prob_map[state(k, 0, i) ] )
    log(WARNING, "i__tau_l_map=\n {}".format(pprint.pformat(i__tau_l_map) ) )
    #
    wing_cutoff_i = 2
    wing_cutoff_sum = 0
    for s, p in state_prob_map.items():
      split_l = s.split(",")
      if int(split_l[1].split("(")[1] ) > wing_cutoff_i or int(split_l[2].split(")")[0] ) > wing_cutoff_i:
        wing_cutoff_sum += p
      
    s_l, p_l = [], []
    for s, p in state_prob_map.items():
      if p > 0.01:
        s_l.append(s)
        p_l.append(p)
    plot.bar(range(len(p_l) ), p_l, color=next(color) )
    plot.xticks([i+0.5 for i in range(len(s_l) ) ], s_l, size='small')
    plot.xlabel("State")
    plot.ylabel("Steady-state probability")
    plot.title(r't= {}, $\lambda$= {}, [$\alpha$, $\beta$, $\gamma$]= {}, sum_on_plot= {}, wing_cutoff_sum= {}'. \
      format(t, "{0:.2f}".format(arr_rate), pprint.pformat(qmu_l), "{0:.2f}".format(sum(p_l)), "{0:.2f}".format(wing_cutoff_sum) ) )
    plot.savefig("plot_simplex_steady_state_prob_hist_ar_{0:.2f}.png".format(arr_rate) )
    plot.clf()

def test_avq(num_f_run, arr_rate, mu, k, r, t, qmu_l=[], w_sys=True, mixed_traff=False,
             p_i_l= [] ):
  E_T_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "arr_rate= {}, mu= {}, k= {}, r= {}, t= {}, qmu_l= {}, w_sys= {}, mixed_traff= {}". \
        format(arr_rate, mu, k, r, t, pprint.pformat(qmu_l), w_sys, mixed_traff) )
    
    num_q = int(1 + t*r) if w_sys else int(t*r)
    qid_l = range(num_q)
    if len(qmu_l) == 0: qmu_l = num_q*[mu]
    
    env = simpy.Environment()
    if mixed_traff:
      sym__rgroup_l_map = {}
      if t == 1: sym_l = ['a', 'b']
      elif t == 3: sym_l = ['a', 'b', 'c']
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
        sym__rgroup_l_map[sym] = rgroup_l
      pg = MT_PG(env, _id="pg", inter_arr_dist=lambda: random.expovariate(arr_rate), sym_l=sym_l)
      log(WARNING, "sym__rgroup_l_map=\n {}".format(pprint.pformat(sym__rgroup_l_map) ) )
      avq = MT_AVQ("mt_avq", env, qid_l, qmu_l, sym__rgroup_l_map)
      # monitor = AVQMonitor(env, aq=avq, poll_dist=lambda: 0.1)
      # avq.join_q.out_m = monitor
      pg.out = avq
    else:
      pg = PG(env, _id="pg", inter_arr_dist=lambda: random.expovariate(arr_rate) )
      avq = AVQ("avq", env, k, r, t, qid_l, qmu_l, w_sys=w_sys)
      # monitor = AVQMonitor(env, aq=avq, poll_dist=lambda: 0.1)
      # avq.join_q.out_m = monitor
      pg.out = avq
      pg.init()
    env.run(until=50000)
    if mixed_traff:
      print("pg.sym__n_sent= {}".format(pprint.pformat(pg.sym__n_sent) ) )
    
    st_l = avq.join_sink.st_l
    if len(st_l) > 0:
      E_T_f_sum += float(sum(st_l) )/len(st_l)
      # continue
    # print("avq.join_sink.qid__num_win_map= {}".format(pprint.pformat(avq.join_sink.qid__num_win_map) ) )
    total_n_wins = sum([n for i, n in avq.join_sink.qid__num_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.join_sink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
    if not mixed_traff:
      total_n_types = sum(avq.type__num_m)
      p_i_l[:] = [n/total_n_types for t,n in enumerate(avq.type__num_m) ]
      print("p_i_l= {}".format(p_i_l) )
    """
    print("\n")
    # print("avq.join_q.state__num_found_map= {}".format(pprint.pformat(avq.join_q.state__num_found_map) ) )
    # total_num_founds = sum([n for s, n in avq.join_q.state__num_found_map.items() ] )
    # state__found_freq_map = {s:float(n)/total_num_founds for s, n in avq.join_q.state__num_found_map.items() }
    # print("state__found_freq_map= {}".format(pprint.pformat(state__found_freq_map) ) )
    
    print("\n")
    # print("monitor.polled_state__counter_map= {}".format(pprint.pformat(monitor.polled_state__counter_map) ) )
    total_counter = sum([c for rs, c in monitor.polled_state__counter_map.items() ] )
    polled_state__counter_map = {rs:float(c)/total_counter for rs, c in monitor.polled_state__counter_map.items() }
    print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
    
    print("\n")
    # print("monitor.state__num_found_by_job_departed_map= {}".format(pprint.pformat(monitor.state__num_found_by_job_departed_map) ) )
    total_counter = sum([c for rs, c in monitor.state__num_found_by_job_departed_map.items() ] )
    state__freq_found_by_job_departed_map = {rs:float(c)/total_counter for rs, c in monitor.state__num_found_by_job_departed_map.items() }
    print("state__freq_found_by_job_departed_map= {}".format(pprint.pformat(state__freq_found_by_job_departed_map) ) )
    
    print("\n")
    # print("monitor.start_setup__num_found_by_job_departed_map= {}".format(pprint.pformat(monitor.start_setup__num_found_by_job_departed_map) ) )
    total_counter = sum([c for rs, c in monitor.start_setup__num_found_by_job_departed_map.items() ] )
    start_setup__freq_found_by_job_departed_map = {rs:float(c)/total_counter for rs, c in monitor.start_setup__num_found_by_job_departed_map.items() }
    print("start_setup__freq_found_by_job_departed_map= {}".format(pprint.pformat(start_setup__freq_found_by_job_departed_map) ) )
    """
  E_T = E_T_f_sum/num_f_run
  print(">> E_T= {}".format(E_T) )
  if E_T > 30: return None
  return E_T
  
def plot_winning_freqs():
  k, r, t = 2, 2, 1
  mu = 1
  arr_rate_ub = simplex_inner_bound_on_arr_rate(r, t, mu, w_sys=True)
  log(WARNING, "k= {}, r= {}, t= {}, mu= {}, arr_rate_ub={}".format(k, r, t, mu, arr_rate_ub) )
  arr_rate_l = []
  qid__win_freq_l_map = {}
  for arr_rate in numpy.linspace(0.05, arr_rate_ub*1.1, 20):
    env = simpy.Environment()
    pg = PG(env, _id="pg",
            inter_arr_dist=lambda: random.expovariate(arr_rate) )
    num_q = 1 + t*r
    qid_l = ["{}".format(i) for i in range(num_q) ]
    avq = AVQ("avq", env, k, r, t, qid_l, qmu_l=[mu for i in range(num_q) ] )
    
    pg.out = avq
    pg.init()
    # monitor = AVQMonitor(env, aq=avq, poll_dist=lambda: 1)
    
    env.run(until=50000)
    
    total_n_wins = sum([n for i, n in avq.join_sink.qid__num_win_map.items() ] )
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in avq.join_sink.qid__num_win_map.items() }
    print("arr_rate= {}, qid__win_freq_map= {}".format(arr_rate, pprint.pformat(qid__win_freq_map) ) )
    
    arr_rate_l.append(arr_rate)
    for qid, win_freq in qid__win_freq_map.items():
      if qid not in qid__win_freq_l_map:
        qid__win_freq_l_map[qid] = []
      qid__win_freq_l_map[qid].append(win_freq)
  
  plot.axhline(y=0.6, label=r'lower-bound, $w_0$', c=next(dark_color), lw=2, ls='--')
  plot.axhline(y=0.4, label=r'upper-bound, $w_1$ or $w_2$', c=next(dark_color), lw=2, ls='--')
  counter = 0
  for qid, win_freq_l in qid__win_freq_l_map.items():
    if counter == 0:
      plot.plot(arr_rate_l, win_freq_l, label=r'sim,$w_0$', color=next(dark_color), marker=next(marker), ms=8, mew=2, ls=':')
    else:
      plot.plot(arr_rate_l, win_freq_l, label=r'sim,$w_1$ or $w_2$', color=next(dark_color), marker=next(marker), ms=8, mew=2, ls=':')
    counter += 1
  plot.legend()
  plot.xlabel(r'Arrival rate $\lambda$', fontsize=12)
  plot.ylabel("Winning frequency", fontsize=12)
  plot.title(r'Simplex($t=1$), $\gamma=\alpha=\beta$= {}'.format(mu) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  # plot.savefig("plot_winning_freqs.png", bbox_inches='tight')
  plot.savefig("plot_winning_freqs.pdf", dpi=fig.dpi)
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_simplex_vs_rep():
  t = 3
  mu = 1 # serv rate at each server for simplex
  if t == 1: arr_rate_ub = 1.6
  elif t == 3: arr_rate_ub = 2.4
  elif t == 7:
    arr_rate_ub = simplex_inner_bound_on_arr_rate(r=2, t=t, mu=mu, w_sys=True)
    arr_rate_ub = float(1.1*arr_rate_ub)
  mixed_traff = False # True
  if mixed_traff: arr_rate_ub = 1.1*arr_rate_ub
  log(WARNING, "t= {}, mu= {}, arr_rate_ub= {}, mixed_traff= {}".format(t, mu, arr_rate_ub, mixed_traff) )
  
  n = 2*t + 1
  n_sym = int(numpy.log2(n+1) )
  # # Same distance
  # n_rep = t + 1
  # n_total_rep = n_sym*n_rep
  # mu_rep = n*mu/n_total_rep
  
  # n_mds = n_sym + t
  # k_mds = n_sym
  # mu_mds = (2*t+1)*mu/n_mds
  # arr_rate_ub_mds = None
  # if t == 3 and not mixed_traff: arr_rate_ub_mds = arr_rate_ub + 0.15 # mds_exact_bound_on_arr_rate(mu_mds, n_mds, k_mds)
  
  # Preserving hot-cold data mix
  # n_rep = t + 1
  # n_total_rep = n_rep
  # arr_rate_ub_mds = None
  
  # Same repair bandwidth
  n_rep = t + 1
  n_total_rep = int(n_sym*(t+1)/2)
  mu_rep = n*mu/n_total_rep if not mixed_traff else n*mu/n_total_rep/n_sym
  arr_rate_ub_mds = None
  
  arr_rate_ub_rep = n_rep*mu_rep
  
  sim_simplex_reqed = False
  E_T_sim_simplex_l = []
  if not mixed_traff and t == 1:
    E_T_sim_simplex_l= [
      0.6775872854372559,
      0.7909557937247363,
      0.9486987202221493,
      1.166209238915134,
      1.5685720588787688,
      2.478342315521276,
      2.6376081306859107,
      2.906788473547391,
      3.263700392764921,
      3.5974807041868426,
      4.289127887822366,
      4.794525358984301,
      5.896928018871929,
      8.099664758903687,
      12.74155958739236]
  elif mixed_traff and t == 1:
    E_T_sim_simplex_mixed_traff_l= [
      0.6795142458623882,
      0.7748927520953908,
      0.9120551663968248,
      1.1017354073281063,
      1.4008309793905753,
      2.0319166972531395,
      2.3461415096416802,
      2.617752845887241,
      2.931842457820586,
      3.3957906721917803,
      4.275140545352988,
      5.384652265631004,
      8.289396804081276,
      None, # 21.85423973012918,
      None]
  elif not mixed_traff and t == 3:
    E_T_sim_simplex_l= [
      0.4676519075931255,
      0.5247256264186801,
      0.6230081386991332,
      0.775814486873029,
      1.0207917160021767,
      1.6244613243247372,
      1.7481208563178903,
      1.9667165686859327,
      2.163968348080258,
      2.5923594863306776,
      3.0700378671376627,
      3.796384731111067,
      4.841880170965622,
      6.610367379250164,
      13.559429107437742]
  elif mixed_traff and t == 3:
    E_T_sim_simplex_mixed_traff_l= [
      0.46628732795742817,
      0.5184094604634668,
      0.5975473670434864,
      0.7272615729604553,
      0.9228862984361961,
      1.3432430706439402,
      1.5297012938889547,
      1.7382202900329649,
      2.006828591863818,
      2.409746021676913,
      2.9987862815607667,
      4.1494167022302415,
      6.7589082110731376,
      None,
      None]
  elif not mixed_traff and t == 7:
    E_T_sim_simplex_l= [
      0.31868938934489865,
      0.3650196292881234,
      0.4281058344507201,
      0.5206469367259021,
      0.6957249200007437,
      1.1325417176453465,
      1.2307386079673424,
      1.3867025010207843,
      1.5768489395874896,
      1.865829597118924,
      2.1844400783734677,
      2.89287730113055,
      4.276904798075734,
      6.184072327220002,
      None]
  else:
    sim_simplex_reqed = True
  
  sim_mds_reqed = False
  E_T_sim_mds_l = []
  if t == 3:
    E_T_sim_mds_l= [
      0.4291382378049635,
      0.4859752967032978,
      0.5573834220518918,
      0.6504572423217563,
      0.7912534680581111,
      1.0617796194912665,
      1.1173955998468372,
      1.1864819039768486,
      1.3132561853089193,
      1.4183354786680833,
      1.5441924947724337,
      1.6800188501504796,
      1.97388257061194,
      2.365205967704707,
      2.552714259149294]
  else:
    sim_mds_reqed = True
  
  sim_mds_split_to_one_reqed = False
  E_T_sim_split_to_one_mds_l = []
  if t == 3:
    E_T_sim_split_to_one_mds_l= [
      0.77365082603341717,
      0.82440222647912942,
      0.88499585518811741,
      0.95059809100622572,
      1.026735997953014,
      1.1276811830357545,
      1.2540326440649683,
      1.4212608769595043,
      1.6517287453133336,
      1.9954850953566452,
      2.5853499093220909,
      3.8254183518878659,
      8.5337611351281506,
      None,
      None]
  else:
    sim_mds_split_to_one_reqed = True
  
  mew, ms = 3, 8
  num_f_run = 2
  def plot_rep_to_all():
    # Simplex
    arr_rate_simplex_l = []
    for arr_rate in [*numpy.linspace(0.05, 0.8*arr_rate_ub, 5, endpoint=False), *numpy.linspace(0.8*arr_rate_ub, arr_rate_ub, 10) ]:
      arr_rate_simplex_l.append(arr_rate)
      if sim_simplex_reqed:
        E_T_sim_simplex_l.append(test_avq(num_f_run, arr_rate, mu, k=2, r=2, t=t, w_sys=True, mixed_traff=mixed_traff) )
    print("E_T_sim_simplex_l= {}".format(pprint.pformat(E_T_sim_simplex_l) ) )
    # E_T_sim_simplex_l = [e*n for e in E_T_sim_simplex_l]
    c = next(dark_color)
    label = 'Simplex' # if t != 1 else 'Simplex or MDS'
    plot.plot(arr_rate_simplex_l, E_T_sim_simplex_l, label=label, color=c, marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # stab_lim = E_T_simplex_approx(t, arr_rate, gamma=mu, mu=mu, incremental=True, arr_rate_ub=True)
    # plot.axvline(stab_lim, label="Simplex stability", color=c, linestyle='--')
    # Rep
    arr_rate_rep_l, E_T_rep_n_1_l = [], []
    for arr_rate in numpy.linspace(0.05, arr_rate_ub_rep-0.05, 20):
      arr_rate_rep_l.append(arr_rate)
      E_T_rep_n_1_l.append(E_T_rep_n_1(arr_rate, mu_rep, n_rep) )
    # E_T_rep_n_1_l = [e*n_rep for e in E_T_rep_n_1_l]
    c = next(dark_color)
    plot.plot(arr_rate_rep_l, E_T_rep_n_1_l, label=r'Replication', color=c, marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.axvline(arr_rate_ub_rep, label="Rep stability", color=c, linestyle='--')
    # # MDS
    # if arr_rate_ub_mds is not None:
    #   arr_rate_mds_l = []
    #   for arr_rate in [*numpy.linspace(0.05, 0.7*arr_rate_ub_mds, 5, endpoint=False), *numpy.linspace(0.7*arr_rate_ub_mds, arr_rate_ub, 10, endpoint=False) ]:
    #   # for arr_rate in numpy.linspace(arr_rate_ub_mds, arr_rate_ub_mds, 1):
    #     arr_rate_mds_l.append(arr_rate)
    #     if sim_mds_reqed:
    #       E_T_sim_mds_l.append(test_avq(num_f_run, arr_rate, mu_mds, k=k_mds, r=n_mds-1, t=1, w_sys=True) )
    #   print("E_T_sim_mds_l= {}".format(pprint.pformat(E_T_sim_mds_l) ) )
    #   plot.plot(arr_rate_mds_l, E_T_sim_mds_l, label=r'MDS', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  def plot_split_to_one():
    # Simplex
    arr_rate_ub = arr_rate_ub_simplex_split_to_one(t, mu) + 0.1
    log(WARNING, "arr_rate_ub= {}".format(arr_rate_ub) )
    arr_rate_l, E_T_simplex_l = [], []
    for arr_rate in numpy.linspace(0.05, arr_rate_ub, 20):
      arr_rate_l.append(arr_rate)
      E_T_simplex_l.append(E_T_simplex_split_to_one_min(t, arr_rate, mu) )
    label = 'Simplex' # if t != 1 else 'Simplex or MDS'
    plot.plot(arr_rate_l, E_T_simplex_l, label=label, color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # Rep
    arr_rate_ub_rep = n_rep*mu_rep
    arr_rate_l, E_T_rep_l = [], []
    for arr_rate in numpy.linspace(0.05, arr_rate_ub_rep-0.2, 20):
      arr_rate_l.append(arr_rate)
      E_T_rep_l.append(E_T_rep_n_1_split_to_one(arr_rate, mu_rep, n_rep) )
    plot.plot(arr_rate_l, E_T_rep_l, label=r'Replication', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # # MDS
    # if t != 1:
    #   arr_rate_ub_mds = 1.1*arr_rate_ub_rep
    #   arr_rate_l = []
    #   for arr_rate in numpy.linspace(0.05, arr_rate_ub_mds, 15):
    #   # for arr_rate in numpy.linspace(arr_rate_ub, arr_rate_ub_mds, 2):
    #     arr_rate_l.append(arr_rate)
    #     if sim_mds_split_to_one_reqed:
    #       sys = 1/(mu_mds-0.5*arr_rate)
    #       mds = test_mds_n_k(num_f_run, 0.5*arr_rate, mu_mds, n_mds-1, k_mds)
    #       if sys < 0 or sys > 30 or mds is None: E_T = None
    #       else: E_T = 0.5*sys + 0.5*mds
    #       E_T_sim_split_to_one_mds_l.append(E_T)
    #   print("E_T_sim_split_to_one_mds_l= {}".format(pprint.pformat(E_T_sim_split_to_one_mds_l) ) )
    #   plot.plot(arr_rate_l, E_T_sim_split_to_one_mds_l, label=r'MDS', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  plot_rep_to_all()
  scheduling = "Replicate-to-all"
  # plot_split_to_one()
  # scheduling = "Split-to-one"
  plot.legend(prop={'size':12})
  plot.xlabel(r'Arrival rate $\lambda$ (Request/sec)', fontsize=12)
  plot.ylabel(r'Average download time (sec)', fontsize=12)
  # plot.title(r'$t={}, \mu={}$'.format(t, mu) )
  plot.title(r'{} scheduling, $t= {}$'.format(scheduling, t) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_simplex_vs_rep_t_{}_{}.pdf".format(t, scheduling) )
  fig.clear()
  # Energy
  # arr_rate_simplex_l, Energy_simplex_l = [], []
  # for arr_rate in numpy.linspace(0.1, arr_rate_ub, 20):
  #   arr_rate_simplex_l.append(arr_rate)
  #   Energy_simplex_l.append(n/arr_rate)
  # arr_rate_rep_l, Energy_rep_l = [], []
  # for arr_rate in numpy.linspace(0.1, arr_rate_ub_rep, 20):
  #   arr_rate_rep_l.append(arr_rate)
  #   Energy_rep_l.append(n_total_rep/arr_rate)
  # plot.plot(arr_rate_simplex_l, Energy_simplex_l, label='Simplex', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # plot.plot(arr_rate_rep_l, Energy_rep_l, label='Rep', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # plot.legend()
  # plot.xlabel(r'Arrival rate $\lambda$', fontsize=12)
  # plot.ylabel(r'Unit of energy per request', fontsize=12)
  # plot.title(r'$t={}, \mu={}$'.format(t, mu) )
  # fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1., def_size[1]/1.)
  # fig.tight_layout()
  # plot.savefig("plot_simplex_vs_rep_t_{}_energy.pdf".format(t) )
  # fig.clear()
  log(WARNING, "done; scheduling= {}, t= {}".format(scheduling, t) )

def plot_simplex():
  w_sys = True
  t, r, k = 1, 2, 2
  mu = 1
  gamma = mu
  if t == 1:
    arr_rate_ub = 1.6
  elif t == 3:
    arr_rate_ub = 2.4
  elif t == 7:
    arr_rate_ub = simplex_inner_bound_on_arr_rate(r, t, mu, w_sys=True)
    arr_rate_ub = float(1.1*arr_rate_ub)
  else:
    arr_rate_ub = simplex_inner_bound_on_arr_rate(r, t, mu, w_sys)
  mixed_traff = False
  log(WARNING, "w_sys= {}, t= {}, r= {}, k= {}, mu= {}, arr_rate_ub= {}, mixed_traff= {}".format(w_sys, t, r, k, mu, arr_rate_ub, mixed_traff) )
  
  E_T_simplex_sm_l, E_T_sim_simplex_l, E_T_simplex_l, E_T_simplex_alt_l, E_T_simplex_matrix_analytic_l = [], [], [], [], []
  E_T_simplex_best_approx_l, E_T_simplex_better_approx_l, E_T_simplex_naive_approx_l, E_T_simplex_varki_gauri_lb_l = [], [], [], []
  E_T_simplex_sim_based_approx_l = []
  E_T_sim_simplex_mixed_traff_l = []
  
  sim_simplex = False
  if w_sys and t == 1:
    E_T_sim_simplex_l= [
      0.6775872854372559,
      0.7909557937247363,
      0.9486987202221493,
      1.166209238915134,
      1.5685720588787688,
      2.478342315521276,
      2.6376081306859107,
      2.906788473547391,
      3.263700392764921,
      3.5974807041868426,
      4.289127887822366,
      4.794525358984301,
      5.896928018871929,
      8.099664758903687,
      12.74155958739236]
  elif w_sys and t == 3:
    E_T_sim_simplex_l= [
      0.4676519075931255,
      0.5247256264186801,
      0.6230081386991332,
      0.775814486873029,
      1.0207917160021767,
      1.6244613243247372,
      1.7481208563178903,
      1.9667165686859327,
      2.163968348080258,
      2.5923594863306776,
      3.0700378671376627,
      3.796384731111067,
      4.841880170965622,
      6.610367379250164,
      13.559429107437742]
  elif w_sys and t == 7:
    E_T_sim_simplex_l= [
      0.31868938934489865,
      0.3650196292881234,
      0.4281058344507201,
      0.5206469367259021,
      0.6957249200007437,
      1.1325417176453465,
      1.2307386079673424,
      1.3867025010207843,
      1.5768489395874896,
      1.865829597118924,
      2.1844400783734677,
      2.89287730113055,
      4.276904798075734,
      6.184072327220002,
      None]
  else:
    sim_simplex = True
  # mixed_traff
  sim_simplex_mixed_traff = False
  if w_sys and t == 1:
    E_T_sim_simplex_mixed_traff_l= [
      0.6795142458623882,
      0.7748927520953908,
      0.9120551663968248,
      1.1017354073281063,
      1.4008309793905753,
      2.0319166972531395,
      2.3461415096416802,
      2.617752845887241,
      2.931842457820586,
      3.3957906721917803,
      4.275140545352988,
      5.384652265631004,
      8.289396804081276,
      None, # 21.85423973012918,
      None]
  elif w_sys and t == 3:
    E_T_sim_simplex_mixed_traff_l= [
      0.46628732795742817,
      0.5184094604634668,
      0.5975473670434864,
      0.7272615729604553,
      0.9228862984361961,
      1.3432430706439402,
      1.5297012938889547,
      1.7382202900329649,
      2.006828591863818,
      2.409746021676913,
      2.9987862815607667,
      4.1494167022302415,
      6.7589082110731376,
      None,
      None]
  else:
    sim_simplex_mixed_traff = True
  
  num_f_run = 2
  arr_rate_l, arr_rate_approx_l = [], []
  for arr_rate in [*numpy.linspace(0.05, 0.8*arr_rate_ub, 5, endpoint=False), *numpy.linspace(0.8*arr_rate_ub, arr_rate_ub, 10) ]:
    arr_rate_l.append(arr_rate)
    
    p_i_l= []
    if sim_simplex:
      E_T_sim_simplex_l.append(test_avq(num_f_run, arr_rate, mu, k, r, t, w_sys=w_sys, p_i_l=p_i_l) )
    # E_T_simplex_sim_based_approx_l.append(E_T_simplex_approx(t, arr_rate, gamma, mu, p_i_l=p_i_l) )
    
    E_T_simplex_sm_l.append(E_T_simplex_splitmerge(t, arr_rate, mu) )
    if t == 1:
      E_T_simplex_l.append(simplex_w_one_repair__E_T(arr_rate, mu) )
      E_T_simplex_matrix_analytic_l.append(E_T_simplex_w_one_repair__matrix_analytic(t, arr_rate, mu) )
    elif t == 2:
      if w_sys:
        E_T_simplex_alt_l.append(simplex_w_two_repair__E_T(arr_rate, mu, M=2) )
        E_T_simplex_l.append(simplex_w_two_repair__E_T(arr_rate, mu, M=5) )
      else:
        E_T_simplex_l.append(simplex_wo_sys_w_two_repair__E_T(arr_rate, mu) )
    E_T_simplex_naive_approx_l.append(E_T_simplex_approx(t, arr_rate, gamma, mu, naive=True) )
    E_T_simplex_better_approx_l.append(E_T_simplex_approx(t, arr_rate, gamma, mu) )
    E_T_simplex_best_approx_l.append(E_T_simplex_approx(t, arr_rate, gamma, mu, incremental=True) )
    # E_T_simplex_varki_gauri_lb_l.append(E_T_simplex_varki_gauri_lb(t, arr_rate, gamma, mu) )
  arr_rate_mixed_traff_l = []
  for arr_rate in [*numpy.linspace(0.05, 0.8*arr_rate_ub, 5, endpoint=False), *numpy.linspace(0.8*arr_rate_ub, 1.1*arr_rate_ub, 10) ]:
  # for arr_rate in numpy.linspace(arr_rate_ub, 1.1*arr_rate_ub, 3):
    arr_rate_mixed_traff_l.append(arr_rate)
    if sim_simplex_mixed_traff:
      E_T_sim_simplex_mixed_traff_l.append(test_avq(num_f_run, arr_rate, mu, k, r, t, w_sys=w_sys, mixed_traff=True) )
  mew, ms = 3, 8
  def plot_poster():
    # for better looking plot
    arr_rate_approx_l = list(arr_rate_l)
    
    ar = arr_rate_ub+0.03
    arr_rate_approx_l.append(ar)
    E_T_simplex_best_approx_l.append(E_T_simplex_approx(t, ar, gamma, mu, incremental=True) )
    
    plot.plot(arr_rate_l, E_T_sim_simplex_l, label="Rep-to-all, simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(arr_rate_approx_l, E_T_simplex_best_approx_l, label="Rep-to-all, approx", zorder=2, marker=next(marker), color='black', linestyle=':', mew=mew, ms=ms)
  def plot_rep_to_all():
    log(WARNING, "E_T_sim_simplex_l= {}".format(pprint.pformat(E_T_sim_simplex_l) ) )
    label = 'Simulation, fixed-arrivals' if mixed_traff else 'Simulation'
    plot.plot(arr_rate_l, E_T_sim_simplex_l, label=label, marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    # plot.plot(arr_rate_l, E_T_simplex_sim_based_approx_l, label=r'Sim-based approximation', marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    if mixed_traff:
      log(WARNING, "E_T_sim_simplex_mixed_traff_l= {}".format(pprint.pformat(E_T_sim_simplex_mixed_traff_l) ) )
      plot.plot(arr_rate_mixed_traff_l, E_T_sim_simplex_mixed_traff_l, label=r'Simulation, mixed-arrivals', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    else:
      plot.plot(arr_rate_l, E_T_simplex_sm_l, label=r'Split-merge upper-bound', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
      # if t == 1:
      #   plot.plot(arr_rate_l, E_T_simplex_matrix_analytic_l, label=r'Matrix-analytic upper-bound', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
      #   plot.plot(arr_rate_l, E_T_simplex_l, label=r'High-traffic approximation', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
      plot.plot(arr_rate_l, E_T_simplex_naive_approx_l, label=r'Naive approximation', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
      plot.plot(arr_rate_l, E_T_simplex_better_approx_l, label=r'Better approximation', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
      # E_T_simplex_best_approx_l[-1] = None
      plot.plot(arr_rate_l, E_T_simplex_best_approx_l, label=r'Best approximation', zorder=2, marker=next(marker), color='black', linestyle=':', mew=mew, ms=ms)
      # plot.plot(arr_rate_l, E_T_simplex_varki_gauri_lb_l, label=r'$E[\hat{T}_{fast-serial}]$', marker=next(marker), color=next(dark_color), linestyle=':', mew=mew)
    # stab_lim = E_T_simplex_approx(t, arr_rate, gamma, mu, incremental=True, arr_rate_ub=True)
    # plot.axvline(stab_lim, label="Stability limit", color='black', linestyle='--')
    # plot.gca().set_xlim([0, stab_lim+0.1] )
  def plot_split_to_one():
    arr_rate_ub = 0.9*arr_rate_ub_simplex_split_to_one(t, mu)
    log(WARNING, "w_sys= {}, t= {}, mu= {}, arr_rate_ub={}".format(w_sys, t, r, k, mu, arr_rate_ub) )
    
    arr_rate_sto_l, E_T_simplex_sto_l = [], []
    for arr_rate in numpy.linspace(0.05, arr_rate_ub, 20):
      arr_rate_sto_l.append(arr_rate)
      E_T_simplex_sto_l.append(E_T_simplex_split_to_one_min(t, arr_rate, mu) )
    # log(WARNING, "E_T_sim_simplex_l= {}".format(pprint.pformat(E_T_sim_simplex_l) ) )
    # plot.plot(arr_rate_l, E_T_sim_simplex_l, 'k', label=r'Replicate-to-all', marker=next(marker), linestyle=':', mew=mew, ms=ms)
    plot.plot(arr_rate_sto_l, E_T_simplex_sto_l, 'b', label=r'Select-one', marker=next(marker), linestyle=':', mew=mew, ms=ms)
  # plot_rep_to_all()
  plot_poster()
  plot_split_to_one()
  plot.legend(prop={'size':11})
  plot.xlabel(r'Arrival rate $\lambda$ (Request/sec)', fontsize=12)
  plot.ylabel(r'Average download time (sec)', fontsize=13)
  plot.title(r'Servers $\sim Exp(\mu={})$, availability $t={}$'.format(mu, t) )
  # plot.title(r'$\mu= {}, t= {}$'.format(mu, t) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_simplex_t_{}.pdf".format(t) )
  # plot.savefig("plot_simplex__t_{}.png".format(t) )
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
    # sim_E_T_avq_l.append(test_avq(num_f_run, arr_rate, mu, k, r, t, w_sys=w_sys) )
    sim_E_T_simplex_l.append(test_avq(num_f_run, arr_rate, mu, simplex_k, simplex_r, simplex_t, w_sys=w_sys) )
    
    E_T_simplex_sm_l.append(E_T_simplex_sn(simplex_t, arr_rate, mu) )
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
  # opt_map = get_opts(sys.argv[1:] )
  # log(WARNING, "opt_map= {}".format(pprint.pformat(opt_map) ) )
  # num_q = int(opt_map["--num_q"] )
  
  # test_avq(num_f_run=1, arr_rate=0.1, mu=1, k=2, r=2, t=1, mixed_traff=True)
  
  # plot_winning_freqs()
  # plot_simplex()
  plot_simplex_vs_rep()
  
  # plot_mds(num_q)
  # plot_n_mds_2(num_q)
  # plot_avq()
  # plot_simplex_w_varying_serv_rate_alloc(num_q)
