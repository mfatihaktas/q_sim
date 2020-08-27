from simplex_sim import *
from simplex_models import *
from plot_utils import *
from log_utils import *

import numpy as np

def sim_ET(t, ar):
  sdist_m = {'dist': 'Exp', 'mu': 1}
  # n = 1 + 2*t
  sym_l, sym__rgroup_l_m = simplex_sym_l__sym__rgroup_l_m(t)
  # log(WARNING, "sym__rgroup_l_m=\n {}".format(pprint.pformat(sym__rgroup_l_m) ) )
  env = simpy.Environment()
  # pg = MT_PG(env, 'pg', ar, sym_l, fixed=True)
  # avq = MT_AVQ('mt_avq', env, n, sym__rgroup_l_m, sdist_m, sching='rep-to-all')
  pg = PG(env, "pg", ar)
  r, k = 2, len(sym_l)
  print("k= {}".format(k) )
  avq = AVQ("avq", env, t, r, k, sdist_m)
  pg.out = avq
  pg.init()
  env.run(until=5000)
  st_l = avq.jsink.st_l
  return np.mean(st_l)

def ar_upperBound(t, sdist_m):
  return float(1/ES_typei(t, t, sdist_m))

def plot_ET_bounds():
  sdist_m = {'dist': 'Exp', 'mu': 1}
  mew, ms = 1, 4
  def plot_(type_, t, c, m):
    ar_l, ET_l = [], []
    for ar in np.linspace(0.1, ar_upperBound(t, sdist_m), 30):
      ET = None
      if type_ == 'UB':
        ET = ET_simplex_sm(t, ar, sdist_m)
      elif type_ == 'Approx':
        ET = ET_simplex_approx(t, ar, sdist_m, incremental=True)[0]
      elif type_ == 'LB':
        ET = ET_simplex_lb(t, ar, sdist_m)
      elif type_ == 'Sim':
        ET = sim_ET(t, ar)
      print("type= {}, ar= {}, ET= {}".format(type_, ar, ET) )
      if ET is None or ET > MAX_ET:
        break
      ar_l.append(ar)
      ET_l.append(ET)
    ls = '-'
    if type_ == 'UB':
      ls = '-.'
    elif type_ == 'Approx':
      ls = '--'
    elif type_ == 'LB':
      ls = ':'
    plot.plot(ar_l, ET_l, label=r'{}, $t= {}$'.format(type_, t), color=c, marker=m, mew=mew, ms=ms, linestyle=ls)

  def plot_t(t):
    print("\n>> t= {}".format(t) )
    c, m = next(dark_color_c), next(marker_c)
    plot_('LB', t, c, m)
    plot_('Approx', t, c, m)
    plot_('UB', t, c, m)
    plot_('Sim', t, c, m)

  # plot_t(1)
  plot_t(3)
  # plot_t(7)
  
  plot.legend(prop={'size':12} )
  plot.xlabel(r'Arrival rate', fontsize=12)
  plot.ylabel(r'Average download time', fontsize=12)
  plot.title(r'Locality $r=2$, Server service rate $\mu = 1$', fontsize=12)
  fig = plot.gcf()
  fig.set_size_inches(2*5, 2*4)
  plot.savefig("plot_ET_bounds.pdf", bbox_inches='tight')
  fig.clear()
  log(INFO, "done.")

if __name__ == "__main__":
  plot_ET_bounds()
