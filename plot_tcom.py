from simplex_sim import *
from simplex_models import *
from plot_utils import *
from log_utils import *

import numpy as np

def sim_ET(t, ar):
  sdist_m = {'dist': 'Exp', 'mu': 1}
  env = simpy.Environment()
  pg = PG(env, "pg", ar)
  avq = AVQ("avq", env, t, 2, 2, sdist_m)
  pg.out = avq
  pg.init()
  env.run(until=50000)
  st_l = avq.jsink.st_l
  return np.mean(st_l)
  # if t = 1:
  #   k = 3
  # test_avq(nf, ar, t, r, k, serv, servdist_m, w_sys=w_sys, p_i_l=p_i_l)

def ar_upperBound(t, sdist_m):
  return float(1/ES_typei(t, t, sdist_m))

t__ar_ETsim_m_l_m = {
  1 : [
    {'ar': 0.1, 'ET': 0.7102228730216715},
    {'ar': 0.16551724137931034, 'ET': 0.7295542761645167},
    {'ar': 0.23103448275862068, 'ET': 0.7584656670519933},
    {'ar': 0.296551724137931, 'ET': 0.8086295000436953},
    {'ar': 0.3620689655172413, 'ET': 0.8278738423288505},
    {'ar': 0.4275862068965517, 'ET': 0.8582962220339132},
    {'ar': 0.49310344827586206, 'ET': 0.9083779256689641},
    {'ar': 0.5586206896551723, 'ET': 0.9495163696128596},
    {'ar': 0.6241379310344827, 'ET': 1.0120631619743001},
    {'ar': 0.689655172413793, 'ET': 1.0631885006518174},
    {'ar': 0.7551724137931034, 'ET': 1.135267336099861},
    {'ar': 0.8206896551724137, 'ET': 1.2255770024332167},
    {'ar': 0.886206896551724, 'ET': 1.3060730608720186},
    {'ar': 0.9517241379310344, 'ET': 1.3749550740790109},
    {'ar': 1.0172413793103448, 'ET': 1.5283061313932638},
    {'ar': 1.0827586206896551, 'ET': 1.7199315573212008},
    {'ar': 1.1482758620689655, 'ET': 1.91994669410001},
    {'ar': 1.2137931034482758, 'ET': 2.149631325370201},
    {'ar': 1.2793103448275862, 'ET': 2.4114674798294393},
    {'ar': 1.3448275862068966, 'ET': 2.999086698296884},
    {'ar': 1.410344827586207, 'ET': 3.539488548614376},
    {'ar': 1.475862068965517, 'ET': 4.73643544753308},
    {'ar': 1.5413793103448274, 'ET': 7.00546384177934},
    {'ar': 1.6068965517241378, 'ET': 14.03306774829778},
    {'ar': 1.6724137931034482, 'ET': 119.30288364751723}
  ],
  3 : [
    {'ar': 0.1, 'ET': 0.47605063351567284},
    {'ar': 0.23448275862068965, 'ET': 0.4947770386657556},
    {'ar': 0.36896551724137927, 'ET': 0.5248109935350759},
    {'ar': 0.503448275862069, 'ET': 0.5502518506418188},
    {'ar': 0.6379310344827586, 'ET': 0.5778066304583614},
    {'ar': 0.7724137931034482, 'ET': 0.6171552772859458},
    {'ar': 0.9068965517241379, 'ET': 0.6579595723125362},
    {'ar': 1.0413793103448277, 'ET': 0.7089548114956226},
    {'ar': 1.1758620689655173, 'ET': 0.766297648548337},
    {'ar': 1.3103448275862069, 'ET': 0.8606130813424128},
    {'ar': 1.4448275862068964, 'ET': 0.9336376992524446},
    {'ar': 1.5793103448275863, 'ET': 1.0720855439293964},
    {'ar': 1.7137931034482758, 'ET': 1.1991951232816207},
    {'ar': 1.8482758620689654, 'ET': 1.4360654118997662},
    {'ar': 1.9827586206896552, 'ET': 1.7301464883276272},
    {'ar': 2.117241379310345, 'ET': 2.5910808935248357},
    {'ar': 2.2517241379310344, 'ET': 3.9680826663437503},
    {'ar': 2.386206896551724, 'ET': 14.116445371638155},
    {'ar': 2.5206896551724136, 'ET': 770.1011042526911} ]
}

def plot_ET_bounds():
  sdist_m = {'dist': 'Exp', 'mu': 1}
  mew, ms = 2, 5
  def plot_(type_, t, c, m):
    ar_l, ET_l = [], []
    if (type_ == 'Sim') and (t in t__ar_ETsim_m_l_m):
      for ar_ETsim_m in t__ar_ETsim_m_l_m[t]:
        ar, ET = ar_ETsim_m['ar'], ar_ETsim_m['ET']
        print("type= {}, ar= {}, ET= {}".format(type_, ar, ET) )
        if ET is None or ET > MAX_ET:
          break
        ar_l.append(ar)
        ET_l.append(ET)
    else:
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
    ls = (0, (1, 10))
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

  plot_t(1)
  plot_t(3)
  # plot_t(7)

  fontsize = 16
  plot.yscale('log')
  plot.legend(prop={'size':12} )
  plot.xlabel(r'Arrival rate', fontsize=fontsize)
  plot.ylabel(r'Average download time', fontsize=fontsize)
  plot.title(r'Server service rate $\mu = 1$, Locality $r=2$, Availability $t$', fontsize=fontsize)
  fig = plot.gcf()
  fig.set_size_inches(7, 5) # 2*5, 2*4
  plot.savefig("plot_FJFA_r2_t.pdf", bbox_inches='tight')
  fig.clear()
  log(INFO, "done.")

if __name__ == "__main__":
  plot_ET_bounds()
