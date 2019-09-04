import numpy as np
from simplex_models import simplex_sym_l__sym__rgroup_l_m
from simplex_sim import *

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
    env.run(until=1*5000) # 50000
    st_l = avq.jsink.st_l
    return np.mean(st_l)
  
  ar_l, ET_mixed_l, ET_fixed_l = [], [], []
  for ar in np.linspace(ar_ub/10, ar_ub, 10):
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

if __name__ == "__main__":
  plot_select_one_w_mixed_fixed()
