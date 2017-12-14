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
  
def plot_mds():
  n, k = 4, 2
  mu = 1
  ar_ub = mds_exact_bound_on_ar(mu, n, k) - 0.05 # mds_inner_bound_on_ar
  log(WARNING, "n= {}, k= {}, mu= {}, ar_ub={}".format(n, k, mu, ar_ub) )
  
  ar_l = []
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
  for ar in [*numpy.linspace(0.05, 0.8*ar_ub, 5, endpoint=False), *numpy.linspace(0.8*ar_ub, ar_ub, 7) ]:
    ar_l.append(ar)
    
    E_T_mds_n_k_sm_l.append(E_T_mds_n_k_sm(ar, mu, n, k) )
    p_i_l = []
    if sim_mds_n_k:
      E_T_mds_n_k_sim_l.append(sim_mds_n_k(num_f_run, ar, mu, n, k, p_i_l=p_i_l) )
    E_T_mds_n_k_sim_based_approx_l.append(E_T_mds_n_k(ar, mu, n, k, p_i_l=p_i_l) )
    E_T_mds_n_k_approx_l.append(E_T_mds_n_k(ar, mu, n, k, incremental=True) )
    
    E_T_mds_n_k_varki_gauri_lb_l.append(E_T_mds_n_k_varki_gauri_lb(ar, mu, n, k) )
  plot.plot(ar_l, E_T_mds_n_k_sm_l, label=r'Split-merge upper-bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  print("E_T_mds_n_k_sim_l= {}".format(pprint.pformat(E_T_mds_n_k_sim_l) ) )
  plot.plot(ar_l, E_T_mds_n_k_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  plot.plot(ar_l, E_T_mds_n_k_sim_based_approx_l, label=r'Sim-based approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  plot.plot(ar_l, E_T_mds_n_k_approx_l, label=r'Best approximation', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  plot.plot(ar_l, E_T_mds_n_k_varki_gauri_lb_l, label=r'Fast-serial lower-bound', color=next(dark_color), marker=next(marker), linestyle=':', mew=2)
  
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