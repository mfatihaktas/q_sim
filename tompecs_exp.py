import numpy as np
from simplex_models import simplex_sym_l__sym__rgroup_l_m, tompecs_sym_l__sym__rgroup_l_m
from simplex_sim import *
from tompecs_sim import *

ms, mew = 10, 0.1

def plot_select_one_w_mixed_fixed():
  t, r, k = 1, 2, 2
  mu = 1
  servdist_m = {'dist': "Exp", 'mu': mu}
  ar_ub = (t + 1)*mu
  log(WARNING, "t= {}, r= {}, k= {}, servdist_m= {}, ar_ub= {}".format(t, r, k, servdist_m, ar_ub) )
  
  sym_l, sym__rgroup_l_m = simplex_sym_l__sym__rgroup_l_m(t)
  log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
  
  def sim(ar, fixed=False):
    env = simpy.Environment()
    pg = MT_PG(env, 'pg', ar, sym_l, fixed=fixed)
    avq = MT_AVQ('mt_avq', env, t, sym__rgroup_l_m, servdist_m, sching='select-one')
    pg.out = avq
    pg.init()
    env.run(until=1*50000)
    st_l = avq.jsink.st_l
    return np.mean(st_l)
  
  ar_l, ET_mixed_l, ET_fixed_l = [], [], []
  for ar in np.linspace(ar_ub/10, ar_ub, 10):
    print("ar= {}".format(ar) )
    ET_mixed = sim(ar)
    ET_fixed = sim(ar, fixed=True)
    print("ET_mixed= {}, ET_fixed= {}".format(ET_mixed, ET_fixed) )
    
    if ET_fixed > 100:
      break
    
    ar_l.append(ar)
    ET_mixed_l.append(ET_mixed)
    ET_fixed_l.append(ET_fixed)
    
  
  plot.plot(ar_l, ET_mixed_l, label='Mixed', color='green', marker='^', ms=ms, mew=mew, ls=':')
  plot.plot(ar_l, ET_fixed_l, label='Fixed', color='blue', marker='v', ms=ms, mew=mew, ls=':')
  
  plot.legend(fontsize=14, framealpha=0.25, numpoints=1)
  plot.xlabel(r'$\lambda$', fontsize=20)
  plot.ylabel(r'$E[T]$', fontsize=20)
  plot.title(r'$t= {}$, $r= {}$, $\mu= {}$'.format(t, r, mu), fontsize=20)
  fig = plot.gcf()
  fig.set_size_inches(5, 4)
  fig.tight_layout()
  plot.savefig("plot_select_one_w_mixed_fixed.pdf".format(r), bbox_inches='tight')
  plot.gcf().clear()

def exp_MixedFJQ():
  n = 2
  V = Exp(1)
  
  def sim(ar, K, u=1):
    X = Exp(ar)
    env = simpy.Environment()
    jg = JobGen(env, X, K, u)
    mfjq = MixedFJQ('mfjq', env, n, V)
    
    jg.out = mfjq
    jg.init()
    env.run(until=1*50000)
    
    k__T_l_m = collections.defaultdict(list)
    rt_l = []
    for _, info_m in mfjq.jid_info_m.items():
      ft_l = info_m['finish_t_l']
      k = info_m['k']
      if len(ft_l) == k:
        T = max(ft_l) - info_m['entry_t'] 
        k__T_l_m[k].append(T)
        rt_l.append(T)
    
    blog(k_nj_m={k: len(T_l) for k, T_l in k__T_l_m.items() } )
    return np.mean(rt_l), \
           {k: np.mean(T_l) for k, T_l in k__T_l_m.items() }
  
  for ar in [0.5, 0.7, 0.9]:
    print("\n>> ar= {}".format(ar) )
    
    K = ExplicitRV(v_l=[1, 2], p_l=[0.5, 0.5] )
    ET, k_ET_m = sim(ar, K)
    log(INFO, "K= {}".format(K), ET=ET, k_ET_m=k_ET_m)
    
    K = ExplicitRV(v_l=[2], p_l=[1] )
    # ET = sim(ar, K, u=3/4) # Båcelli LB
    ET, k_ET_m = sim(ar, K)
    log(INFO, "K= {}".format(K), ET=ET, k_ET_m=k_ET_m)

def tompecs_sim(sim=False):
  mu = 1
  servdist_m = {'dist': 'Exp', 'mu': mu}
  num_sim_secs = 50000
  def sim(scheme, ar):
    env = simpy.Environment()
    if scheme == '3rep':
      pg = PG(env, 'pg', ar)
      mdsq = MDSQ(scheme, env, 1, range(3), servdist_m)
      pg.out = mdsq
      pg.init()
      env.run(until=num_sim_secs)
      st_l = mdsq.jsink.st_l
    elif scheme == 'mds':
      pg = PG(env, 'pg', ar)
      mdsq = MDSQ(scheme, env, 6, range(9), servdist_m)
      pg.out = mdsq
      pg.init()
      env.run(until=num_sim_secs)
      st_l = mdsq.jsink.st_l
    elif scheme == 'avail':
      t, n = 3, 7
      sym_l, sym__rgroup_l_m = simplex_sym_l__sym__rgroup_l_m(t)
      # log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
      pg = MT_PG(env, 'pg', ar, sym_l)
      avq = MT_AVQ(scheme, env, n, sym__rgroup_l_m, servdist_m, sching='rep-to-all')
      pg.out = avq
      pg.init()
      env.run(until=num_sim_secs)
      st_l = avq.jsink.st_l
    elif scheme == 'lrc':
      n = 10
      sym_l, sym__rgroup_l_m = tompecs_sym_l__sym__rgroup_l_m(scheme)
      # log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
      pg = MT_PG(env, 'pg', ar, sym_l)
      avq = MT_AVQ(scheme, env, n, sym__rgroup_l_m, servdist_m, sching='rep-to-all')
      pg.out = avq
      pg.init()
      env.run(until=num_sim_secs)
      st_l = avq.jsink.st_l
    return np.mean(st_l)
  
  mu_ = 10*mu/3
  scheme_conf_m = {
    '3rep': {
      'fixed_mu': {'ar_ub': 3*mu},
      'fixed_cum_mu': {'ar_ub': 3*mu_} },
    'mds': {
      'fixed_mu': {'ar_ub': 1.5*mu/(H(9) - H(3) ) },
      'fixed_cum_mu': {'ar_ub': 1.5*mu_/(H(9) - H(3) ) } },
    'avail': {
      'fixed_mu': {'ar_ub': 1.3*mu},
      'fixed_cum_mu': {'ar_ub': 1.3*mu_} },
    'lrc': {
      'fixed_mu': {'ar_ub': 1.3*mu},
      'fixed_cum_mu': {'ar_ub': 1.3*mu_} }
  }
  
  ## num_ar = 15
  scheme_conf_data_m = {
    '3rep': {
      'legend': '3-Replication',
      'fixed_mu': {
        'ar_l': [0.33999999999999997, 0.5342857142857143, 0.7285714285714284, 0.9228571428571427, 1.117142857142857, 1.3114285714285714, 1.5057142857142853, 1.6999999999999997, 1.8942857142857141, 2.088571428571428, 2.2828571428571425, 2.477142857142857, 2.671428571428571, 2.865714285714285],
        'ET_l': [0.38087510119710039, 0.40142627419722948, 0.44416643101733305, 0.47710790733290581, 0.5333616757044769, 0.59394996094310215, 0.67501479549798293, 0.76528358848143996, 0.91796223965083379, 1.111062708757246, 1.3785347556191252, 2.0011637457425362, 2.8352368688931233, 8.9362648565232394]
      },
      'fixed_cum_mu': None
    },
    'mds': {
      'legend': '(9, 6)-MDS',
      'fixed_mu': {
        'ar_l': [0.15065763252291753, 0.24750896771622166, 0.3443603029095258, 0.4412116381028299, 0.538062973296134, 0.6349143084894382, 0.7317656436827422, 0.8286169788760464, 0.9254683140693505, 1.0223196492626547, 1.1191709844559588, 1.2160223196492628, 1.312873654842567, 1.4097249900358713],
        'ET_l': [1.0584343014665178, 1.1013319413402058, 1.1534956882872194, 1.2288357073700271, 1.3100227696102535, 1.3887998227289138, 1.5161028119928401, 1.6417803016527981, 1.8509519318920735, 2.1006750336059623, 2.3897852030034756, 2.925355453556334, 4.1137560616368356, 7.746239729204734]
      },
      'fixed_cum_mu': None
    },
    'avail': {
      'legend': '(7, 3, 2, 3)-Availability',
      'fixed_mu': {
        'ar_l': [0.13, 0.21357142857142858, 0.29714285714285715, 0.38071428571428567, 0.46428571428571425, 0.5478571428571428, 0.6314285714285713, 0.715, 0.7985714285714285, 0.882142857142857, 0.9657142857142856, 1.0492857142857142, 1.1328571428571426, 1.2164285714285712],
        'ET_l': [1.0056224909803788, 1.0446995774288599, 1.096113217911292, 1.1603675363015422, 1.232460676398611, 1.3333538085078303, 1.4339629447176452, 1.6296353980493312, 1.7907076296155595, 2.1400256521158005, 2.6665740981795478, 3.526287002904418, 6.2037472561774551, 37.786638952680157]
      },
      'fixed_cum_mu': None
    },
    'lrc': {
      'legend': '(10, 6, 3, 1)-LRC'
      'fixed_mu': {
        'ar_l': [0.13, 0.21357142857142858, 0.29714285714285715, 0.38071428571428567, 0.46428571428571425, 0.5478571428571428, 0.6314285714285713, 0.715, 0.7985714285714285, 0.882142857142857, 0.9657142857142856],
        'ET_l': [1.274830117513327, 1.3767787115877124, 1.471499357217972, 1.6441380257049791, 1.7913720588206556, 2.0761145968756618, 2.42166657998707, 3.0667730866185332, 4.2185288254439168, 7.6515824222759878, 31.452833014928334]
      },
      'fixed_cum_mu': None
    }
  }
  
  sim = True
  num_ar = 15
  if sim:
    def sim_(scheme):
      ar_ub = scheme_conf_m[scheme][conf]['ar_ub']
      ar_l, ET_l = [], []
      for ar in np.linspace(ar_ub/10, ar_ub, num_ar):
        ET = sim(scheme, ar)
        print("{}, ar= {}, ET= {}".format(scheme, ar, ET) )
        # ET_model = 1/(3*mu - ar)
        # print("{}, ar= {}, ET= {}, ET_model= {}".format(scheme, ar, ET, ET_model) )
        if ET > 100:
          break
        ar_l.append(ar)
        ET_l.append(ET)
      scheme_conf_data_m[scheme][conf] = {'ar_l': ar_l, 'ET_l': ET_l}
    
    
    sim_('3rep')
    sim_('mds')
    sim_('avail')
    sim_('lrc')
  log(INFO, "done; ", scheme_conf_data_m=scheme_conf_data_m)
  return scheme_conf_data_m

def tompecs_exp():
  conf = 'fixed_cum_mu'
  scheme_conf_data_m = tompecs_sim()
  
  def plot_(data, legend):
    print("len(data['ar_l'])= {}, len(data['ET_l'])= {}".format(len(data['ar_l']), len(data['ET_l']) ) )
    plot.plot(data['ar_l'], data['ET_l'], label=data['legend'], color=next(dark_color), marker=next(marker), ms=ms, mew=mew, ls=':')
    
  plot_(scheme_conf_data_m['3rep'][conf], scheme_conf_data_m['3rep']['legend'] )
  plot_(scheme_conf_data_m['mds'][conf], scheme_conf_data_m['3rep']['legend'] )
  plot_(scheme_conf_data_m['avail'][conf], scheme_conf_data_m['3rep']['legend'] )
  plot_(scheme_conf_data_m['lrc'][conf], scheme_conf_data_m['3rep']['legend'] )
  
  plot.legend(fontsize=12, framealpha=0.25, numpoints=1)
  # plot.xscale('log')
  plot.xlabel(r'$\lambda$', fontsize=20)
  plot.ylabel(r'$E[T]$', fontsize=20)
  plot.title(r'$\mu= {}$'.format(mu), fontsize=20)
  fig = plot.gcf()
  fig.set_size_inches(6, 4)
  fig.tight_layout()
  plot.savefig("tompecs_{}.pdf".format(conf), bbox_inches='tight')
  plot.gcf().clear()
  log(INFO, "done.")

if __name__ == "__main__":
  # plot_select_one_w_mixed_fixed()
  # exp_MixedFJQ()
  
  tompecs_exp()
