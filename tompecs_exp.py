import numpy as np
from simplex_models import simplex_sym_l__sym__rgroup_l_m, tompecs_sym_l__sym__rgroup_l_m
from simplex_sim import *
from tompecs_sim import *

from simplex_models import ET_simplex_approx

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
    env.run(until=50000)
    
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

def tompecs_sim(conf='fixed_mu', pop='uniform'):
  num_sim_secs = 50000 # 500
  log(INFO, "", conf=conf, pop=pop, num_sim_secs=num_sim_secs)
  
  def sim_(scheme, ar, mu):
    ## k = 6 across all the schemes given below.
    servdist_m = {'dist': 'Exp', 'mu': mu}
    env = simpy.Environment()
    if scheme == 'rep':
      # pg = PG(env, 'pg', ar)
      # mdsq = MDSQ(scheme, env, 1, range(3), servdist_m)
      # pg.out = mdsq
      # pg.init()
      # env.run(until=num_sim_secs)
      # st_l = mdsq.jsink.st_l
      ET_model = lambda ar: 1/(3*mu - ar)
      
      if pop == 'uniform':
        st_l = [ET_model(ar/6) ]
      elif pop == 'skewed':
        p_l = scheme_conf_m[scheme]['pop_l']
        st_l = [sum(p*ET_model(p*ar) for p in p_l) ]
    elif scheme == 'mds':
      # pg = PG(env, 'pg', ar)
      # mdsq = MDSQ(scheme, env, 6, range(9), servdist_m)
      # pg.out = mdsq
      # pg.init()
      # env.run(until=num_sim_secs)
      # st_l = mdsq.jsink.st_l
      t, n = 3, 7
      sym_l, sym__rgroup_l_m = tompecs_sym_l__sym__rgroup_l_m(scheme)
      # log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
      p_l = scheme_conf_m[scheme]['pop_l'] if pop == 'skewed' else None
      pg = MT_PG(env, 'pg', ar, sym_l, p_l)
      avq = MT_AVQ(scheme, env, n, sym__rgroup_l_m, servdist_m, sching='rep-to-all')
      pg.out = avq
      pg.init()
      env.run(until=num_sim_secs)
      st_l = avq.jsink.st_l
    elif scheme == 'avail':
      t, n = 3, 7
      sym_l, sym__rgroup_l_m = simplex_sym_l__sym__rgroup_l_m(t)
      # log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
      if pop == 'skewed':
        p_l = scheme_conf_m[scheme]['pop_l']
        ar = ar*0.925
      else:
        p_l = None
        ar = ar*0.5
      pg = MT_PG(env, 'pg', ar, sym_l, p_l)
      avq = MT_AVQ(scheme, env, n, sym__rgroup_l_m, servdist_m, sching='rep-to-all')
      pg.out = avq
      pg.init()
      env.run(until=num_sim_secs)
      st_l = avq.jsink.st_l
    elif scheme == 'lrc':
      n = 10
      sym_l, sym__rgroup_l_m = tompecs_sym_l__sym__rgroup_l_m(scheme)
      # log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
      p_l = scheme_conf_m[scheme]['pop_l'] if pop == 'skewed' else None
      pg = MT_PG(env, 'pg', ar, sym_l, p_l)
      avq = MT_AVQ(scheme, env, n, sym__rgroup_l_m, servdist_m, sching='rep-to-all')
      pg.out = avq
      pg.init()
      env.run(until=num_sim_secs)
      st_l = avq.jsink.st_l
    return np.mean(st_l)
  mu = 1
  scheme_conf_m = {
    'rep': {
      'fixed_mu': {'mu': mu, 'ar_ub': 3}, # 18*mu
      'fixed_cum_mu': {'mu': mu, 'ar_ub': 2.5},
      'pop_l': [0.45, 0.45, 0.025, 0.025, 0.025, 0.025] },
    'mds': {
      'fixed_mu': {'mu': mu, 'ar_ub': 1.5*mu/(H(9) - H(3) ) },
      'fixed_cum_mu': {'mu': 18*mu/9, 'ar_ub': 1.5*(18*mu/9)/(H(9) - H(3) ) },
      'pop_l': [0.45, 0.45, 0.025, 0.025, 0.025, 0.025] },
    'avail': {
      'fixed_mu': {'mu': mu, 'ar_ub': 2*1.3*mu},
      'fixed_cum_mu': {'mu': 18*mu/14, 'ar_ub': 2*1.3*(18*mu/14) },
      'pop_l': [0.45/0.925, 0.45/0.925, 0.025/0.925] },
    'lrc': {
      'fixed_mu': {'mu': mu, 'ar_ub': 1.3*mu},
      'fixed_cum_mu': {'mu': 18*mu/10, 'ar_ub': 1.3*(18*mu/10) },
      'pop_l': [0.45, 0.45, 0.025, 0.025, 0.025, 0.025] }
  }
  
  ## num_ar = 15
  scheme_conf_data_m = {
    'rep': {
      'legend': '3-Replication',
      'fixed_mu': {
        'uniform': {
          'ar_l': [0.29999999999999999, 0.49285714285714288, 0.68571428571428572, 0.87857142857142856, 1.0714285714285714, 1.2642857142857142, 1.4571428571428573, 1.6500000000000001, 1.842857142857143, 2.0357142857142856, 2.2285714285714286, 2.4214285714285713, 2.6142857142857143, 2.8071428571428569, 3.0],
          'ET_l': [0.33898305084745761, 0.34271725826193394, 0.34653465346534656, 0.35043804755944929, 0.35443037974683544, 0.35851472471190782, 0.36269430051813473, 0.36697247706422015, 0.37135278514588854, 0.37583892617449666, 0.38043478260869568, 0.38514442916093533, 0.38997214484679671, 0.39492242595204513, 0.40000000000000002] },
        'skewed': {
          'ET_l': [0.34266733215152412, 0.34897707014492735, 0.35554669262982547, 0.36239267941193748, 0.36953293413173655, 0.37698694143334172, 0.38477594542217064, 0.39292315285325069, 0.4014539651430899, 0.41039624409089553, 0.41978061716489884, 0.42964082940167114, 0.44001415043860287, 0.45094184702727308, 0.46246973365617422],
          'ar_l': [0.20000000000000001, 0.32857142857142863, 0.45714285714285718, 0.58571428571428585, 0.71428571428571441, 0.84285714285714297, 0.97142857142857153, 1.1000000000000001, 1.2285714285714286, 1.3571428571428572, 1.4857142857142858, 1.6142857142857143, 1.7428571428571431, 1.8714285714285717, 2.0] } },
      'fixed_cum_mu': {
        'uniform': {
          'ET_l': [0.33994334277620392, 0.34433285509325678, 0.34883720930232559, 0.35346097201767307, 0.35820895522388063, 0.36308623298033288, 0.36809815950920244, 0.37325038880248834, 0.37854889589905366, 0.38400000000000001, 0.38961038961038963, 0.39538714991762763, 0.40133779264214048, 0.40747028862478779, 0.41379310344827591],
          'ar_l': [0.34999999999999998, 0.57499999999999996, 0.80000000000000004, 1.0249999999999999, 1.25, 1.4750000000000001, 1.7000000000000002, 1.9249999999999998, 2.1499999999999999, 2.375, 2.6000000000000001, 2.8250000000000002, 3.0500000000000003, 3.2750000000000004, 3.5] },
        'skewed': {
          'ar_l': [0.25, 0.4107142857142857, 0.5714285714285714, 0.73214285714285721, 0.8928571428571429, 1.0535714285714286, 1.2142857142857144, 1.375, 1.5357142857142858, 1.6964285714285716, 1.8571428571428572, 2.0178571428571432, 2.1785714285714288, 2.3392857142857144, 2.5],
          'ET_l': [0.34509123444405276, 0.35314334052184349, 0.36161782296650724, 0.37054898788294821, 0.37997496097312172, 0.38993823434967445, 0.40048631005679752, 0.4116724608580935, 0.42355663400313509, 0.43620653033555534, 0.44969889975696475, 0.46412110540425155, 0.47957302389680034, 0.49616936902256559, 0.51404255319148939] } } },
    'mds': {
      'legend': '(9, 6)-MDS',
      'fixed_mu': {
        'uniform': {
          'ET_l': [1.5333904823779256, 1.7085161423962389, 1.9293611894871168, 2.2658476473995752, 2.8230338048991781, 3.7290273003029988, 6.0708663096702216, 27.781474918929675],
          'ar_l': [0.15065763252291753, 0.24750896771622166, 0.3443603029095258, 0.44121163810282987, 0.53806297329613395, 0.63491430848943819, 0.73176564368274222, 0.82861697887604635] },
        'skewed': {
          'ET_l': [1.6195893852637242, 1.7813494241049199, 2.0637209076835812, 2.5623218447503229, 3.3338191252436995, 5.0230585082540777, 12.21869650065298],
          'ar_l': [0.15065763252291753, 0.24750896771622166, 0.3443603029095258, 0.44121163810282987, 0.53806297329613395, 0.63491430848943819, 0.73176564368274222] } },
      'fixed_cum_mu': {
        'uniform': {
          'ET_l': [0.77712014987372491, 0.86385632106340349, 0.98061685336259963, 1.131033976948401, 1.3761173935515827, 1.8737189440900888, 2.9574850117434401, 9.8070402792881417],
          'ar_l': [0.30131526504583506, 0.49501793543244332, 0.68872060581905159, 0.88242327620565975, 1.0761259465922679, 1.2698286169788764, 1.4635312873654844, 1.6572339577520927] },
        'skewed': {
          'ET_l': [0.7861860847331209, 0.9075667241434624, 1.0451637467034165, 1.2716745579784223, 1.704423325698331, 2.5846700954386508, 6.6836034057193885],
          'ar_l': [0.30131526504583506, 0.49501793543244332, 0.68872060581905159, 0.88242327620565975, 1.0761259465922679, 1.2698286169788764, 1.4635312873654844] } } },
    'avail': {
      'legend': '(7, 3, 2, 3)-Availability',
      'fixed_mu': {
        'uniform': {
          'ar_l': [0.26000000000000001, 0.42714285714285716, 0.59428571428571431, 0.76142857142857134, 0.92857142857142849, 1.0957142857142856, 1.2628571428571427, 1.4299999999999999, 1.597142857142857, 1.764285714285714, 1.9314285714285713, 2.0985714285714283, 2.2657142857142851, 2.4328571428571424],
          'ET_l': [0.99151209604698776, 1.0293698624379681, 1.0880961443988266, 1.1436574391321774, 1.2355785160484958, 1.3275081265214035, 1.4508908878772917, 1.5853011064625546, 1.8022659297945116, 2.1305986899227469, 2.5758567295476942, 3.6001300716019515, 6.401997408761873, 15.445710033736498] },
        'skewed': {
          'ET_l': [0.98009193402916783, 1.0408543210452406, 1.0736673438938225, 1.1375971125938973, 1.1989604092723134, 1.281260478582974, 1.3586011742108923, 1.5289417163617531, 1.6777137276258194, 1.9197425372770547, 2.235420042476552, 2.6776825724737865, 3.6965931822764917, 5.9754152130855873, 18.219033801860782],
          'ar_l': [0.13, 0.21357142857142858, 0.29714285714285715, 0.38071428571428567, 0.46428571428571425, 0.54785714285714282, 0.63142857142857134, 0.71499999999999997, 0.79857142857142849, 0.88214285714285701, 0.96571428571428564, 1.0492857142857142, 1.1328571428571426, 1.2164285714285712, 1.3] } },
      'fixed_cum_mu': {
        'uniform': {
          'ET_l': [0.78625370964887098, 0.81635332012828288, 0.85721288347072921, 0.88748585820053505, 0.96383318601264845, 1.0283706353657904, 1.1155905133620945, 1.2498401571747293, 1.4328234912400721, 1.6849768581145876, 2.0942931205138722, 2.6500959006144118, 4.8378123351862268, 13.052248545330935],
          'ar_l': [0.33428571428571435, 0.54918367346938779, 0.76408163265306128, 0.97897959183673477, 1.1938775510204083, 1.4087755102040818, 1.6236734693877553, 1.8385714285714287, 2.0534693877551025, 2.268367346938776, 2.4832653061224494, 2.6981632653061229, 2.9130612244897964, 3.1279591836734699] },
        'skewed': {
          'ET_l': [0.77249693560595012, 0.80199566977573478, 0.83425860516296191, 0.87847802523360619, 0.93092260647036185, 1.000168693044555, 1.070421724980273, 1.1786212388866923, 1.2712600099136284, 1.5214894794790617, 1.7437434408635732, 2.1403519110795628, 2.8259757336114446, 4.4041803581702874, 13.094126218504449],
          'ar_l': [0.16714285714285718, 0.27459183673469389, 0.38204081632653064, 0.48948979591836739, 0.59693877551020413, 0.70438775510204088, 0.81183673469387763, 0.91928571428571437, 1.0267346938775512, 1.134183673469388, 1.2416326530612247, 1.3490816326530615, 1.4565306122448982, 1.563979591836735, 1.6714285714285717] } } },
    'lrc': {
      'legend': '(10, 6, 3, 1)-LRC',
      'fixed_mu': {
        'uniform': {
          'ET_l': [1.3047182027125728, 1.3920610755049221, 1.5049300467270408, 1.6538432854533918, 1.8414571836460372, 2.0421162425415353, 2.4015501614920924, 2.9641342114200921, 4.0219863688899569, 8.9009141143498649, 29.385447222791157],
          'ar_l': [0.13, 0.21357142857142858, 0.29714285714285715, 0.38071428571428567, 0.46428571428571425, 0.54785714285714282, 0.63142857142857134, 0.71499999999999997, 0.79857142857142849, 0.88214285714285701, 0.93571428571428564] }, # 
        'skewed': {
          'ET_l': [1.2950227013879616, 1.4088161681394722, 1.5352948330270402, 1.7073792786537867, 1.9214536962973803, 2.2504472467868695, 2.7093621609544805, 3.515786152791804, 6.0505105257448699, 16.030441106489754],
          'ar_l': [0.13, 0.21357142857142858, 0.29714285714285715, 0.38071428571428567, 0.46428571428571425, 0.54785714285714282, 0.63142857142857134, 0.71499999999999997, 0.79857142857142849, 0.88214285714285701] } },
      'fixed_cum_mu': {
        'uniform': {
          'ET_l': [0.71588795709136077, 0.76956136844671541, 0.82942158895648743, 0.90620344429281396, 1.0114846741854473, 1.1491498682186312, 1.3558488947888701, 1.716060829174932, 2.4878736354564759, 4.3976876990294631], # 47.442065962192885
          'ar_l': [0.23400000000000004, 0.38442857142857145, 0.53485714285714292, 0.68528571428571428, 0.83571428571428585, 0.98614285714285721, 1.1365714285714286, 1.2870000000000001, 1.4374285714285715, 1.5878571428571429] }, # 1.7382857142857144
        'skewed': {
          'ET_l': [0.71607224651990675, 0.77523012516793088, 0.84561143816536422, 0.95261256845382813, 1.0687576291179439, 1.2745689414623611, 1.5433748476769666, 2.1034635366769341, 3.342497915826228, 13.374641615387274],
          'ar_l': [0.23400000000000004, 0.38442857142857145, 0.53485714285714292, 0.68528571428571428, 0.83571428571428585, 0.98614285714285721, 1.1365714285714286, 1.2870000000000001, 1.4374285714285715, 1.5878571428571429] } } }
  }
  
  sim = False # True
  num_ar = 15
  if sim:
    def sim__(scheme):
      mu = scheme_conf_m[scheme][conf]['mu']
      ar_ub = scheme_conf_m[scheme][conf]['ar_ub']
      ar_l, ET_l = [], []
      for ar in np.linspace(ar_ub/10, ar_ub, num_ar):
        ET = sim_(scheme, ar, mu)
        print("{}, ar= {}, ET= {}".format(scheme, ar, ET) )
        # ET_model = 1/(3*mu - ar)
        # print("{}, ar= {}, ET= {}, ET_model= {}".format(scheme, ar, ET, ET_model) )
        if ET > 100 or ET < 0:
          break
        ar_l.append(ar)
        ET_l.append(ET)
      scheme_conf_data_m[scheme][conf][pop] = {'ar_l': ar_l, 'ET_l': ET_l}
      log(INFO, "scheme= {}".format(scheme), ar_l=ar_l, ET_l=ET_l)
    
    sim__('rep')
    # sim__('mds')
    # sim__('avail')
    # sim__('lrc')
  # log(INFO, "done; ", scheme_conf_data_m=scheme_conf_data_m)
  return scheme_conf_data_m

def tompecs_exp():
  conf = 'fixed_mu' # 'fixed_cum_mu'
  pop = 'uniform' # 'skewed'
  scheme_conf_data_m = tompecs_sim(conf, pop)
  
  def plot_(data, legend):
    print("len(data['ar_l'])= {}, len(data['ET_l'])= {}".format(len(data['ar_l']), len(data['ET_l']) ) )
    plot.plot(data['ar_l'], data['ET_l'], label=legend, color=next(dark_color), marker=next(marker), ms=ms, mew=mew, ls=':')
    
  plot_(scheme_conf_data_m['rep'][conf][pop], scheme_conf_data_m['rep']['legend'] )
  plot_(scheme_conf_data_m['mds'][conf][pop], scheme_conf_data_m['mds']['legend'] )
  plot_(scheme_conf_data_m['avail'][conf][pop], scheme_conf_data_m['avail']['legend'] )
  plot_(scheme_conf_data_m['lrc'][conf][pop], scheme_conf_data_m['lrc']['legend'] )
  
  plot.legend(loc='best', fontsize=12, framealpha=0.25, numpoints=1) # loc='upper left'
  # plot.xscale('log')
  plot.xlabel(r'Arrival rate $\lambda$', fontsize=18)
  plot.ylabel('Average download time', fontsize=18)
  title = 'Cumulative service rate is fixed to 10' if conf == 'fixed_cum_mu' else "Each server's service rate is fixed to 1"
  title += '\nUniform object popularity' if pop == 'uniform' else '\nSkewed object popularity'
  plot.title(title, fontsize=16)
  fig = plot.gcf()
  fig.set_size_inches(6, 4)
  fig.tight_layout()
  plot.savefig("tompecs_{}_{}.pdf".format(conf, pop), bbox_inches='tight')
  plot.gcf().clear()
  log(INFO, "done.")

def plot_servtype_prob():
  t, r, k = 1, 2, 2
  mu = 1
  servdist_m = {'dist': 'Exp', 'mu': mu}
  log(INFO, "", t=t, r=r, k=k, servdist_m=servdist_m)
  
  def sim_(ar):
    env = simpy.Environment()
    pg = PG(env, 'pg', ar)
    avq = AVQ('avq', env, t, r, k, servdist_m)
    monitor = AVQMonitor(env, avq, poll_rate=10)
    avq.join_q.out_m = monitor
    pg.out = avq
    pg.init()
    env.run(until=50000) # 50000
    ET = np.mean(avq.jsink.st_l)
    if ET > 100:
      log(ERROR, "ET= {} > 100!".format(ET) )
      return
    
    load = np.mean(monitor.polled_busystate_l)
    total_n_types = sum(avq.servtype_num_l)
    return [n/total_n_types for n in avq.servtype_num_l], load
  
  if t == 1: ar_ub = 1.6
  elif t == 3: ar_ub = 2.4
  elif t == 7: ar_ub = float(1.1*reptoall_innerbound_on_ar(mu, t, r, w_sys=True) )
  
  def bins_labels(bins, **kwargs):
    bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
    plot.xticks(np.arange(min(bins)+bin_w/2, max(bins), bin_w), bins, **kwargs)
    plot.xlim(bins[0], bins[-1])
  
  ar_l = [ar_ub*0.3, ar_ub*0.6, ar_ub*0.9]
  barwidth = 0.1
  fig, ax = plot.subplots()
  ax.set_xticks(np.arange(t+1) + barwidth)
  # ax.set_xticks(np.arange(t+1))
  ax.set_xticklabels([str(i) for i in range(t+1) ] )
  
  plot.plot(np.NaN, np.NaN, label=r'Approximate', marker='D', markerfacecolor='None', markeredgecolor='black', ms=6, mew=1, ls='None')
  for i, ar in enumerate(ar_l):
    fi_l, load = sim_(ar)
    
    _, pi_l = ET_simplex_approx(t, ar, servdist_m, incremental=True)
    log(INFO, "ar= {}".format(ar), load=load, fi_l=fi_l, pi_l=pi_l)
    
    c = next(dark_color)
    plot.bar(np.arange(t+1)+i*barwidth, fi_l, label=r'$\rho= {}$'.format(round(load, 2) ), color=c, width=barwidth, alpha=0.5)
    # bins_labels(range(t+1), fontsize=16)
    # plot.bar(np.arange(t+1)+i*barwidth, pi_l, label=r'$\rho= {}$, approx'.format(round(load, 2) ), color=c, width=barwidth, alpha=1)
    plot.plot(np.arange(t+1)+i*barwidth, pi_l, color=c, marker='D', markeredgecolor='black', ms=6, mew=1, ls='None')

  fontsize = 16
  plot.legend(fontsize=13, framealpha=0.25, numpoints=1)
  plot.xlabel('Request service type', fontsize=fontsize)
  plot.ylabel(r'$f_i$', fontsize=fontsize)
  k = int(math.log2(t+1))+1
  plot.title(r'$n= {}$, $k= {}$, $r= {}$, $t= {}$, $\mu= {}$'.format(2**k-1, k, r, t, mu), fontsize=fontsize)
  fig = plot.gcf()
  fig.set_size_inches(6, 4)
  fig.tight_layout()
  plot.savefig("plot_servtype_prob_t{}.pdf".format(t), bbox_inches='tight')
  plot.gcf().clear()
  log(INFO, "done.")

if __name__ == "__main__":
  # plot_select_one_w_mixed_fixed()
  # exp_MixedFJQ()
  
  # tompecs_exp()
  plot_servtype_prob()
