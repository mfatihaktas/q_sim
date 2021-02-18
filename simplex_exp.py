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

from simplex_sim import *
from simplex_models import *
from mds_models import mds_exactbound_on_ar
# from mds_exp import sim_mds_nk

from scipy.interpolate import UnivariateSpline

def plot_reptoall_steadystate_probhist():
  t, r, k = 1, 2, 2
  def get_state_prob_m(ar):
    log(WARNING, "ar= {}, t= {}, r= {}, k= {}".format(ar, t, r, k) )
    env = simpy.Environment()
    pg = PG(env, "pg", ar)
    avq = AVQ("avq", env, t, r, k, serv="Exp", servdist_m={'mu': 1} )
    # monitor = AVQMonitor(env, avq, poll_dist=lambda: 0.1)
    # avq.join_q.out_m = monitor
    pg.out = avq
    env.run(until=50000)
    
    # print("monitor.polled_state__counter_map= {}".format(pprint.pformat(monitor.polled_state__counter_map) ) )
    total_counter = sum([c for rs, c in monitor.polled_state__counter_map.items() ] )
    state_prob_m = {rs:float(c)/total_counter for rs, c in monitor.polled_state__counter_map.items() }
    # print("polled_state__counter_map= {}".format(pprint.pformat(polled_state__counter_map) ) )
    return state_prob_m # ['0,(0,0)']
  # for ar in numpy.arange(0.05, 1.2, 0.1):
  color = iter(cm.rainbow(numpy.linspace(0, 1, 20) ) )
  plot.figure(figsize=(20,10) )
  for ar in numpy.arange(0.05, 1.3, 0.1):
  # for ar in numpy.arange(0.05, 0.1, 0.1):
    state_prob_m = get_state_prob_m(ar)
    
    def state(kp, i, j):
      return "{},({},{})".format(kp, i, j)
    i__tau_l_map = {}
    for i in range(10):
      if i not in i__tau_l_map:
        i__tau_l_map[i] = []
      for kp in range(i, 10):
        s_u, s_l = state(kp, i, 0), state(kp+1, i, 0)
        if s_u in state_prob_m and s_l in state_prob_m:
          i__tau_l_map[i].append(state_prob_m[s_l]/state_prob_m[s_u] )
        # if state(k+1, 0, i) in state_prob_m:
        #   i__tau_l_map[i].append(state_prob_m[state(k+1, 0, i) ] /state_prob_m[state(k, 0, i) ] )
    log(WARNING, "i__tau_l_map=\n {}".format(pprint.pformat(i__tau_l_map) ) )
    #
    wing_cutoff_i = 2
    wing_cutoff_sum = 0
    for s, p in state_prob_m.items():
      split_l = s.split(",")
      if int(split_l[1].split("(")[1] ) > wing_cutoff_i or int(split_l[2].split(")")[0] ) > wing_cutoff_i:
        wing_cutoff_sum += p
      
    s_l, p_l = [], []
    for s, p in state_prob_m.items():
      if p > 0.01:
        s_l.append(s)
        p_l.append(p)
    plot.bar(range(len(p_l) ), p_l, color=next(color) )
    plot.xticks([i+0.5 for i in range(len(s_l) ) ], s_l, size='small')
    plot.xlabel("State")
    plot.ylabel("Steady-state probability")
    plot.title(r't= {}, $\lambda$= {}, sum_on_plot= {}, wing_cutoff_sum= {}'. \
      format(t, "{0:.2f}".format(ar), "{0:.2f}".format(sum(p_l)), "{0:.2f}".format(wing_cutoff_sum) ) )
    plot.savefig("plot_reptoall_steadystate_probhist_ar_{0:.2f}.png".format(ar) )
    plot.clf()

def test_avq(nf, ar, t, r, k, serv="Exp", servdist_m=None,
             w_sys=True, mixed_traff=False, sching="rep-to-all", p_i_l= [] ):
  E_T_f_sum = 0
  for f in range(nf):
    log(WARNING, "ar= {}, t= {}, r= {}, k= {}, servdist_m= {}, w_sys= {}, mixed_traff= {}, sching= {}". \
        format(ar, t, r, k, servdist_m, w_sys, mixed_traff, sching) )
    
    env = simpy.Environment()
    if mixed_traff:
      sym_l, sym__rgroup_l_m = simplex_sym_l__sym__rgroup_l_m(t)
      log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
      pg = MT_PG(env, "pg", ar, sym_l)
      avq = MT_AVQ("mt_avq", env, t, sym__rgroup_l_m, serv, servdist_m)
      # monitor = AVQMonitor(env, aq=avq, poll_dist=lambda: 0.1)
      # avq.join_q.out_m = monitor
    else:
      psize = None
      if serv == "Bern*Pareto":
        psize = "Pareto"
        serv = "Bern"
      pg = PG(env, "pg", ar, psize=psize, psize_dist_m=servdist_m)
      avq = AVQ("avq", env, t, r, k, servdist_m, sching, w_sys=w_sys)
      # monitor = AVQMonitor(env, aq=avq, poll_dist=lambda: 0.1)
      # avq.join_q.out_m = monitor
    pg.out = avq
    pg.init()
    c = 3 if serv == "Pareto" or serv == "Bern" else 1
    env.run(until=c*50000) # 20
    
    if mixed_traff:
      print("pg.sym__n_sent= {}".format(pprint.pformat(pg.sym__n_sent) ) )
    st_l = avq.jsink.st_l
    if len(st_l) > 0:
      E_T_f_sum += float(sum(st_l) )/len(st_l)
      # continue
    # print("avq.jsink.qid__num_win_map= {}".format(pprint.pformat(avq.jsink.qid__num_win_map) ) )
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    qid_winfreq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_win_map.items() }
    print("qid_winfreq_map= {}".format(pprint.pformat(qid_winfreq_map) ) )
    # if not mixed_traff:
    #   total_n_types = sum(avq.servtype__num_m)
    #   p_i_l[:] = [n/total_n_types for t, n in enumerate(avq.servtype__num_m) ]
    #   print("p_i_l= {}".format(p_i_l) )
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
  E_T = E_T_f_sum/nf
  print(">> E_T= {}".format(E_T) )
  if E_T > 100: return None
  return E_T
  
def plot_winning_freqs():
  t, r, k = 1, 2, 2
  mu = 1
  servdist_m = {'dist': 'Exp', 'mu': mu}
  ar_ub = reptoall_innerbound_on_ar(t, servdist_m)
  log(WARNING, "t= {}, servdist_m= {}, ar_ub={}".format(t, servdist_m, ar_ub) )
  ar_l = []
  qid__winfreq_l_map = {}
  for ar in numpy.linspace(0.05, ar_ub*1.1, 20):
    env = simpy.Environment()
    pg = PG(env, "pg", ar)
    avq = AVQ("avq", env, t, r, k, servdist_m, "rep-to-all")
    pg.out = avq
    pg.init()
    # monitor = AVQMonitor(env, aq=avq, poll_dist=lambda: 1)
    env.run(until=50000)
    
    total_n_wins = sum([n for i, n in avq.jsink.qid__num_win_map.items() ] )
    qid_winfreq_map = {i:float(n)/total_n_wins for i, n in avq.jsink.qid__num_win_map.items() }
    print("ar= {}, qid_winfreq_map= {}".format(ar, pprint.pformat(qid_winfreq_map) ) )
    
    ar_l.append(ar)
    for qid, win_freq in qid_winfreq_map.items():
      if qid not in qid__winfreq_l_map:
        qid__winfreq_l_map[qid] = []
      qid__winfreq_l_map[qid].append(win_freq)
  
  plot.axhline(y=0.6, label=r'Lower-bound, $w_s$', c=next(dark_color), lw=2, ls='--')
  plot.axhline(y=0.4, label=r'Upper-bound, $w_r$', c=next(dark_color), lw=2, ls='--')
  counter = 0
  for qid, win_freq_l in qid__winfreq_l_map.items():
    if counter == 0:
      plot.plot(ar_l, win_freq_l, label=r'Simulation, $w_s$', color=next(dark_color), marker=next(marker), ms=8, mew=2, ls=':')
    else:
      plot.plot(ar_l, win_freq_l, label=r'Simulation, $w_r$', color=next(dark_color), marker=next(marker), ms=8, mew=2, ls=':')
    counter += 1

  fontsize = 16
  plot.legend(fontsize=13)
  plot.xlabel(r'Arrival rate $\lambda$', fontsize=fontsize)
  plot.ylabel("Fraction of request completions", fontsize=fontsize)
  plot.title(r'Replicate-to-all $t=1$, $\gamma=\alpha=\beta= {}$'.format(mu), fontsize=fontsize)
  fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.set_size_inches(6, 4)
  fig.tight_layout()
  # plot.savefig("plot_winning_freqs.png", bbox_inches='tight')
  plot.savefig("plot_winning_freqs.pdf", dpi=fig.dpi)
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_simplex_vs_rep():
  t, r, k = 3, 2, 2
  serv = "Exp"
  mu = 1
  servdist_m['mu'] = mu
  if t == 1: ar_ub = 1.6
  elif t == 3: ar_ub = 2.4
  elif t == 7:
    ar_ub = float(1.1*reptoall_innerbound_on_ar(mu, t, r, w_sys=True) )
  mixed_traff = False
  if mixed_traff: ar_ub = 1.1*ar_ub
  log(WARNING, "t= {}, ar_ub= {}, serv= {}, servdist_m= {}, mixed_traff= {}".format(t, ar_ub, serv, servdist_m, mixed_traff) )
  
  n = 2*t + 1
  n_sym = int(numpy.log2(n+1) )
  # # Same distance
  # n_rep = t + 1
  # n_total_rep = n_sym*n_rep
  # mu_rep = n*mu/n_total_rep
  
  # n_mds = n_sym + t
  # k_mds = n_sym
  # mu_mds = (2*t+1)*mu/n_mds
  # ar_ub_mds = None
  # if t == 3 and not mixed_traff: ar_ub_mds = ar_ub + 0.15 # mds_exactbound_on_ar(mu_mds, n_mds, k_mds)
  
  # Preserving hot-cold data mix
  # n_rep = t + 1
  # n_total_rep = n_rep
  # ar_ub_mds = None
  
  # Same repair bandwidth
  n_rep = t + 1
  n_total_rep = int(n_sym*(t+1)/2)
  mu_rep = n*mu/n_total_rep if not mixed_traff else n*mu/n_total_rep/n_sym
  ar_ub_mds = None
  
  ar_ub_rep = n_rep*mu_rep
  
  sim_simplex_reqed = False
  ET_sim_l = []
  if not mixed_traff and t == 1:
    ET_sim_l= [
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
    ET_sim_mixedtraff_l= [
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
    ET_sim_l= [
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
    ET_sim_mixedtraff_l= [
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
    ET_sim_l= [
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
  nf = 2
  def plot_reptoall():
    # Simplex
    ar_simplex_l = []
    for ar in [*numpy.linspace(0.05, 0.8*ar_ub, 5, endpoint=False), *numpy.linspace(0.8*ar_ub, ar_ub, 10) ]:
      ar_simplex_l.append(ar)
      if sim_simplex_reqed:
        ET_sim_l.append(test_avq(nf, ar, t, r, k, serv, servdist_m, w_sys=True, mixed_traff=mixed_traff) )
    c = next(dark_color)
    label = 'Simplex' # if t != 1 else 'Simplex or MDS'
    print("ET_sim_l= {}".format(pprint.pformat(ET_sim_l) ) )
    plot.plot(ar_simplex_l, ET_sim_l, label=label, color=c, marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # stab_lim = ET_simplex_approx(t, ar, servdist_m, incremental=True, ar_ub=True)
    # plot.axvline(stab_lim, label="Simplex stability", color=c, linestyle='--')
    # Rep
    ar_rep_l, E_T_rep_n_1_l = [], []
    for ar in numpy.linspace(0.05, ar_ub_rep-0.05, 20):
      ar_rep_l.append(ar)
      E_T_rep_n_1_l.append(E_T_rep_n_1(ar, mu_rep, n_rep) )
    # E_T_rep_n_1_l = [e*n_rep for e in E_T_rep_n_1_l]
    c = next(dark_color)
    plot.plot(ar_rep_l, E_T_rep_n_1_l, label=r'Replication', color=c, marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.axvline(ar_ub_rep, label="Rep stability", color=c, linestyle='--')
    # # MDS
    # if ar_ub_mds is not None:
    #   ar_mds_l = []
    #   for ar in [*numpy.linspace(0.05, 0.7*ar_ub_mds, 5, endpoint=False), *numpy.linspace(0.7*ar_ub_mds, ar_ub, 10, endpoint=False) ]:
    #   # for ar in numpy.linspace(ar_ub_mds, ar_ub_mds, 1):
    #     ar_mds_l.append(ar)
    #     if sim_mds_reqed:
    #       E_T_sim_mds_l.append(test_avq(nf, ar, t=1, r, k, serv, {'mu': mu_mds}, w_sys=True) )
    #   print("E_T_sim_mds_l= {}".format(pprint.pformat(E_T_sim_mds_l) ) )
    #   plot.plot(ar_mds_l, E_T_sim_mds_l, label=r'MDS', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  def plot_selectone():
    # Simplex
    ar_ub = arub_simplex_selectone(t, mu) + 0.1
    log(WARNING, "ar_ub= {}".format(ar_ub) )
    ar_l, ET_l = [], []
    for ar in numpy.linspace(0.05, ar_ub, 20):
      ar_l.append(ar)
      ET_l.append(ET_selectone(t, ar, mu) )
    label = 'Simplex' # if t != 1 else 'Simplex or MDS'
    plot.plot(ar_l, ET_l, label=label, color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # Rep
    ar_ub_rep = n_rep*mu_rep
    ar_l, E_T_rep_l = [], []
    for ar in numpy.linspace(0.05, ar_ub_rep-0.2, 20):
      ar_l.append(ar)
      E_T_rep_l.append(E_T_rep_n_1_split_to_one(ar, mu_rep, n_rep) )
    plot.plot(ar_l, E_T_rep_l, label=r'Replication', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  plot_reptoall()
  scheduling = "Replicate-to-all"
  # plot_selectone()
  # scheduling = "Split-to-one"
  plot.legend(prop={'size':12})
  plot.xlabel(r'Arrival rate $\lambda$ (Request/s)', fontsize=12)
  plot.ylabel(r'Average download time (s)', fontsize=12)
  # plot.title(r'$t={}, \mu={}$'.format(t, mu) )
  plot.title(r'{} scheduling, $t= {}$'.format(scheduling, t) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_simplex_vs_rep_t_{}_{}.pdf".format(t, scheduling) )
  fig.clear()
  # Energy
  # ar_simplex_l, Energy_simplex_l = [], []
  # for ar in numpy.linspace(0.1, ar_ub, 20):
  #   ar_simplex_l.append(ar)
  #   Energy_simplex_l.append(n/ar)
  # ar_rep_l, Energy_rep_l = [], []
  # for ar in numpy.linspace(0.1, ar_ub_rep, 20):
  #   ar_rep_l.append(ar)
  #   Energy_rep_l.append(n_total_rep/ar)
  # plot.plot(ar_simplex_l, Energy_simplex_l, label='Simplex', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # plot.plot(ar_rep_l, Energy_rep_l, label='Rep', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
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

def plot_reptoall():
  mixed_traff, w_sys = False, True
  t, r, k = 1, 2, 2
  serv = "Exp" # "Bern" # "Bern*Pareto" # "Pareto" # "Dolly"
  mu = 1
  # loc, a = 1, 2
  # U, L, p, loc, a = 1, 8, 0.2, 0.1, 1.5 # 1, 8, 0.2, 1, 3
  U, L, p, loc, a = 1, 10, 0.3, 0.1, 1.5 # 1, 8, 0.2, 1, 3
  # For rep-to-all
  if serv == "Exp":
    servdist_m = {'dist': serv, 'mu': mu}
    if t == 1: ar_ub = 1.6
    elif t == 3: ar_ub = 2.4
    elif t == 7: ar_ub = float(1.1*reptoall_innerbound_on_ar(t, servdist_m) )
    else: ar_ub = reptoall_innerbound_on_ar(t, servdist_m)
  elif serv == "Pareto":
    servdist_m = {'dist': serv, 'loc': loc, 'a': a}
    ar_ub = reptoall_innerbound_on_ar(t, servdist_m)
  elif serv == "TPareto":
    servdist_m = {'dist': serv, 'l': l, 'u': u, 'a': a}
    ar_ub = reptoall_innerbound_on_ar(t, servdist_m)
  elif serv == "Bern" or serv == "Bern*Pareto":
    servdist_m = {'dist': serv, 'U': U, 'L': L, 'p': p, 'loc': loc, 'a': a}
    ar_ub = reptoall_innerbound_on_ar(t, servdist_m)
  elif serv == "Dolly":
    servdist_m = None
    if t == 1: ar_ub = 0.28
    elif t == 3: ar_ub = 0.4
  log(WARNING, "w_sys= {}, t= {}, r= {}, k= {}, servdist_m= {}, ar_ub= {}, mixed_traff= {}".format(w_sys, t, r, k, servdist_m, ar_ub, mixed_traff) )
  
  ET_sm_l, ET_sim_l, ET_l, ET_lb_l = [], [], [], []
  ET_alt_l, ET_matrixanalytic_l = [], []
  ET_bestapprox_l, ET_betterapprox_l, ET_naiveapprox_l, ET_varkigauri_lb_l = [], [], [], []
  ET_simbasedapprox_l = []
  ET_sim_mixedtraff_l = []
  
  # All below w_sys=True
  nf = 3
  sim_simplex = False
  if serv == "Exp":
    if t == 1:
      ET_sim_l= [
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
    elif t == 3:
      ET_sim_l= [
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
    else: sim_simplex = True
  elif serv == "Pareto":
    if loc == 1 and a == 2:
      if t == 1:
        ET_sim_l= [
          1.5299993522735693,
          1.7233577876041122,
          1.8952577131712123,
          2.2418712080584897,
          2.853623528849504,
          4.2208097489868,
          4.586420599121132,
          5.191481636572133,
          5.6340499086639815,
          5.9712033727746,
          7.94309766204549,
          9.599736059102067,
          13.280357368839619,
          17.20104661693977,
          25.449711725024084]
      elif t == 3:
        ET_sim_l= [
          1.3221090353539466,
          1.4459274633541828,
          1.6229349092564267,
          1.9043964678064051,
          2.4154300633936936,
          3.6666730405584844,
          3.9217550909479577,
          4.256167164955279,
          4.717366068731679,
          5.891743883842969,
          6.04468767433355,
          8.073514650754076,
          9.880581947509592,
          15.816118977624845,
          28.433468299774272]
      else: sim_simplex = True
    elif loc == 1 and a == 5:
      if t == 3:
        ET_sim_l= [
          1.1276007604818075,
          1.240550592912947,
          1.3862061325608057,
          1.645653757532261,
          2.0688083303883276,
          3.2115831386711813,
          3.2986018954384835,
          3.8148027478966227,
          4.033705086448495,
          5.448028336643181,
          5.697392211154507,
          9.053323168666376,
          10.17868048265699,
          23.644561610837382,
          None] # 93.02644300031747
      else: sim_simplex = True
    else: sim_simplex = True
  elif serv == "Bern":
    if U == 1 and L == 8 and p == 0.2:
      if t == 1:
        # nf = 3
        ET_sim_l= [
          1.6376474738985423,
          1.9851446427827089,
          2.4840795375267626,
          3.1829054073054217,
          4.39332366216294,
          7.063110373762194,
          7.4445330550351665,
          8.208129233744382,
          9.309321611480481,
          10.747520637423975,
          12.460023568734707,
          15.038255521201348,
          18.778687793661728,
          23.582209372296532,
          36.21619587757658]
      elif t == 3:
        # nf = 1
        ET_sim_l= [
          1.1072895175117927,
          1.2582695204803385,
          1.4572200912301614,
          1.8340775367273732,
          2.4430722742069184,
          4.053853819806121,
          4.4494192069988605,
          5.061922101782603,
          5.883304533639656,
          6.705043861319703,
          8.307668993372534,
          11.041651319984396,
          17.564101468045756,
          33.184482866801716,
          None]
      else: sim_simplex = True
    else: sim_simplex = True
  elif serv == "Bern*Pareto":
    if U == 1 and L == 8 and p == 0.2 and loc == 1 and a == 3:
      if t == 11:
        # nf = 3
        ET_sim_l= [
          2.142631836594827,
          2.5302711620514966,
          2.941315337537391,
          3.8773353598252345,
          4.550420407107853,
          6.649089020276313,
          7.000687768519389,
          7.681497353358071,
          8.058275694322152,
          9.541434770613856,
          10.136837383356713,
          11.027889242435874,
          14.072462480848941,
          18.721889173565945,
          29.85022801496356]
      elif t == 33:
        pass
      else: sim_simplex = True
    else: sim_simplex = True
  else: sim_simplex = True
  
  # Mixed traff
  sim_simplex_mixed_traff = False
  if mixed_traff:
    if serv == "Exp":
      if t == 1:
        ET_sim_mixedtraff_l= [
          0.678978501641253,
          0.7748022818617738,
          0.9072886738372506,
          1.0928902616368403,
          1.43754904360929,
          2.0810587767368154,
          2.266461910378062,
          2.5977047234601125,
          3.2441553951140985,
          3.585616438620215,
          4.415600179701042,
          6.099149242270735,
          9.786138444920114,
          None, # 21.631079441147904
          None]
      elif t == 3:
        ET_sim_mixedtraff_l= [
          0.46217641274184773,
          0.5249541076176077,
          0.6065798815902482,
          0.7193352388312126,
          0.9238674360581351,
          1.363955390788439,
          1.4654931553890183,
          1.733811055160431,
          2.0493965738680795,
          2.479767271681704,
          3.065826086322138,
          4.300842192226751,
          8.05986376865404,
          None, # 35.70730644518723,
          None]
      else:
        sim_simplex_mixed_traff = True
  
  ar_l = []
  for ar in [*numpy.linspace(0.05, 0.8*ar_ub, 5, endpoint=False), *numpy.linspace(0.8*ar_ub, ar_ub, 10) ]:
  # for ar in numpy.linspace(0.05, ar_ub, 2):
    ar_l.append(ar)
    
    p_i_l = []
    if sim_simplex:
      ET_sim = test_avq(nf, ar, t, r, k, serv, servdist_m, w_sys=w_sys, p_i_l=p_i_l)
      print("*** ET_sim= {}".format(ET_sim) )
      ET_sim_l.append(ET_sim)
      # ET_sim_l.append(None)
    
    # ET_simbasedapprox_l.append(ET_simplex_approx(t, ar, servdist_m, p_i_l=p_i_l)[0] )
    # if sim_simplex_mixed_traff:
    #   ET_sim_mixedtraff_l.append(test_avq(nf, ar, t, r, k, serv, servdist_m, w_sys=w_sys, p_i_l=p_i_l, mixed_traff=True) )
    
    ET_sm_l.append(ET_simplex_sm(t, ar, servdist_m) )
    ET_lb_l.append(ET_simplex_lb(t, ar, servdist_m) )
    if serv == "Exp":
      if t == 1:
        ET_l.append(ET_reptoall_t1(ar, mu) )
        ET_matrixanalytic_l.append(ET_reptoall_t1_matrixanalytic(t, ar, mu) )
      elif t == 2:
        if w_sys:
          ET_alt_l.append(simplex_w_two_repair__E_T(ar, mu, M=2) )
          ET_l.append(simplex_w_two_repair__E_T(ar, mu, M=5) )
        else:
          ET_l.append(simplex_wo_sys_w_two_repair__E_T(ar, mu) )
    ET_naiveapprox_l.append(ET_simplex_approx(t, ar, servdist_m, naive=True)[0] )
    ET_betterapprox_l.append(ET_simplex_approx(t, ar, servdist_m)[0] )
    ET_bestapprox_l.append(ET_simplex_approx(t, ar, servdist_m, incremental=True)[0] )
    # ET_varkigauri_lb_l.append(E_T_simplex_varki_gauri_lb(t, ar, gamma, mu)[0] )
  
  ar_mixed_traff_l = []
  # for ar in numpy.linspace(0.2, 0.2, 1):
  for ar in [*numpy.linspace(0.05, 0.8*ar_ub, 5, endpoint=False), *numpy.linspace(0.8*ar_ub, 1.1*ar_ub, 10) ]:
    ar_mixed_traff_l.append(ar)
    if sim_simplex_mixed_traff:
      ET_sim_mixedtraff_l.append(test_avq(nf, ar, t, r, k, serv, servdist_m, w_sys=w_sys, mixed_traff=True) )
  
  # mew, ms = 0.1, 10
  mew, ms = 2, 5
  def plot_poster():
    # for better looking plot
    ar_approx_l = list(ar_l)
    
    ar = ar_ub + 0.03
    ar_approx_l.append(ar)
    ET_bestapprox_l.append(ET_simplex_approx(t, ar, servdist_m, incremental=True) )
    
    plot.plot(ar_l, ET_sim_l, label="FJ-FA, simulation", marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(ar_approx_l, ET_bestapprox_l, label="FJ-FA, M/G/1 approximation", zorder=2, marker=next(marker), color='black', linestyle=':', mew=mew, ms=ms)

  def get_xs_l_ys_l(_x_l, _y_l):
    x_l, y_l = [], []
    for i, y in enumerate(_y_l):
      if y is not None:
        x_l.append(_x_l[i])
        y_l.append(y)
    
    s = UnivariateSpline(x_l, y_l, s=0.001)
    xs_l = np.linspace(min(x_l), max(x_l), 20)
    ys_l = s(xs_l)
    return xs_l, ys_l
    
  def plot_():
    log(WARNING, "ET_sim_l= {}".format(pprint.pformat(ET_sim_l) ) )
    # plot.plot(ar_l, ET_simbasedapprox_l, label=r'Sim-based approximation', marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    label = 'Simulation, fixed-arrivals' if mixed_traff else 'Simulation'

    xs_l, ys_l = get_xs_l_ys_l(ar_l, ET_sim_l)
    # plot.plot(ar_l, ET_sim_l, label=label, marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    plot.plot(xs_l, ys_l, label=label, marker=next(marker), zorder=1, color=next(dark_color), linestyle=':', mew=mew, ms=ms)
    
    if mixed_traff:
      log(WARNING, "ET_sim_mixedtraff_l= {}".format(pprint.pformat(ET_sim_mixedtraff_l) ) )
      plot.plot(ar_mixed_traff_l, ET_sim_mixedtraff_l, label=r'Simulation, mixed-arrivals', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
    else:
      xs_l, ys_l = get_xs_l_ys_l(ar_l, ET_sm_l)
      # plot.plot(ar_l, ET_sm_l, label=r'Split-Merge upper bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
      plot.plot(xs_l, ys_l, label=r'Split-Merge upper bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
      
      # plot.plot(ar_l, ET_bestapprox_l, label=r'$M/G/1$ approximation', zorder=2, marker=next(marker), color='black', linestyle=':', mew=mew, ms=ms)
      xs_l, ys_l = get_xs_l_ys_l(ar_l, ET_lb_l)
      # plot.plot(ar_l, ET_lb_l, label=r'Fast-Split-Merge lower bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
      plot.plot(xs_l, ys_l, label=r'Fast-Split-Merge lower bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
      if t == 1:
        xs_l, ys_l = get_xs_l_ys_l(ar_l, ET_matrixanalytic_l)
        # plot.plot(ar_l, ET_matrixanalytic_l, label=r'Matrix-analytic upper-bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
        plot.plot(xs_l, ys_l, label=r'Matrix-analytic upper-bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)

        xs_l, ys_l = get_xs_l_ys_l(ar_l, ET_l)
        # plot.plot(ar_l, ET_l, label=r'High-traffic approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
        plot.plot(xs_l, ys_l, label=r'High-traffic approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
      # plot.plot(ar_l, ET_naiveapprox_l, label=r'Straightforward approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
      # plot.plot(ar_l, ET_betterapprox_l, label=r'Better approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew, ms=ms)
      # plot.plot(ar_l, ET_bestapprox_l, label=r'Fine-grained approximation', zorder=2, marker=next(marker), color='black', linestyle=':', mew=mew, ms=ms)
      # plot.plot(ar_l, ET_varkigauri_lb_l, label=r'$E[\hat{T}_{fast-serial}]$', color=next(dark_color), marker=next(marker), linestyle=':', mew=mew)
    # stab_lim = ET_simplex_approx(t, ar, servdist_m, incremental=True, ar_ub=True)
    # plot.axvline(stab_lim, label="Stability limit", color='black', linestyle='--')
    # plot.gca().set_xlim([0, stab_lim+0.1] )
  
  def plot_selectone():
    ar_ub = 0.9*arub_simplex_selectone(t, serv, servdist_m)
    log(WARNING, "ar_ub={}".format(ar_ub) )
    ar_l, ET_l = [], []
    for ar in numpy.linspace(0.05, ar_ub, 50):
    # for ar in numpy.linspace(0.05, ar_ub, 2):
      ar_l.append(ar)
      # if sim:
      #   ET_l.append(test_avq(nf, ar, t, r, k, serv, servdist_m, w_sys=w_sys, sching="select-one") )
      ET_l.append(ET_selectone(t, ar, mu) )
    # log(WARNING, "ET_l= {}".format(pprint.pformat(ET_l) ) )
    plot.plot(ar_l, ET_l, 'b', label=r'Select-one', linestyle='--', lw=3, mew=mew, ms=ms)
  # plot_poster()
  plot_()
  
  # plot.plot(ar_l, ET_sim_l, 'k', label=r'Replicate-to-all', linestyle='-', lw=3)
  # plot_selectone()
  fontsize = 16
  plot.yscale('log')
  plot.legend(loc='upper left', fontsize=13, framealpha=0.25)
  plot.xlabel(r'Arrival rate $\lambda$', fontsize=fontsize)
  plot.ylabel(r'Average download time', fontsize=fontsize)
  serv_in_latex = None
  if serv == "Exp":
    serv_in_latex = '\mathrm{Exp}' + r'(\mu={})'.format(mu)
  elif serv == "Pareto":
    serv_in_latex = r'Pareto(s={}, \alpha={})'.format(loc, a)
  elif serv == "Bern":
    serv_in_latex = r'Bernoulli(U={}, L={}, p={})'.format(U, L, p)
  elif serv == "Dolly":
    serv_in_latex = r'Dolly'
  plot.title(r'FJ-FA with $r= {}$, $t= {}$, $\mu= {}$'.format(r, t, mu), fontsize=fontsize)
  # plot.title(r'$t={}$, Servers $\sim {}$'.format(t, serv_in_latex) )
  fig = plot.gcf()
  fig.set_size_inches(6, 4)
  fig.tight_layout()
  plot.savefig("plot_FJFA_r{}_t{}.pdf".format(r, t) )
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
  
  # plot_winning_freqs()
  plot_reptoall()
  # plot_simplex_vs_rep()
  
