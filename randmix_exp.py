import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

from mixed_sim import MixedPG, MixedNet
from randmix_sim import *
from randmix_model import *
from patch import *
from commonly_used import *

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def sim_nkmix(num_frun, n, k, dist_m):
  ET_sum = 0
  for f in range(num_frun):
    log(WARNING, "n= {}, k= {}".format(n, k) )
    env = simpy.Environment()
    
    pg = MixedPG(env, "pg", [dist_m for i in range(n) ] )
    m = MixedNet(env, n, k)
    pg.out = m
    env.run(until=50000*1)
    
    ET_sum += m.ET_ET2()[0]
  return ET_sum/num_frun

def sim_samplekmix(num_frun, n, k, dist_m):
  ET_sum = 0
  for f in range(num_frun):
    log(WARNING, "n= {}, k= {}".format(n, k) )
    env = simpy.Environment()
    
    pg = MixedPG(env, "pg", [dist_m for i in range(n) ] )
    m = SamplekMix(env, n, k)
    pg.out = m
    env.run(until=50000*1)
    
    ET_sum += m.ET()
  return ET_sum/num_frun

def plot_mix():
  num_frun = 1
  
  def plot_ET_vs_k(n, dist_m):
    k_l = []
    ET_samplekmix_sim_l, ET_samplekmix_l = [], []
    ET_nkmix_sim_l = []
    
    sim = False # True
    # n = 40, mu = 1
    ET_nkmix_sim_l= [0.01314522595292312,
      0.027029970757678906,
      0.041651515422629494,
      0.057184546391097144,
      0.07342264685219294,
      0.09086822482832366,
      0.10947453365352897,
      0.1289673736066132,
      0.15013226305337055,
      0.1724385992827565,
      0.19657959655041332,
      0.22156000755747338,
      0.24995448806724285,
      0.27977123605586557,
      0.31226272433600655,
      0.34774058734582525,
      0.3857361720597433,
      0.4279659888153143,
      0.4760837075111895,
      0.5270868075403595,
      0.5833746252899382,
      0.645674404228407,
      0.7179074631413703,
      0.798376642092802,
      0.8949136414452482,
      1.0009853249422354,
      1.1211971246177983,
      1.2718412178713563,
      1.449197633658389,
      1.6710895578579272,
      1.9381493113271806,
      2.281011251664213,
      2.766378510057035,
      3.398239547959631,
      4.400850644625362,
      6.084681257925894,
      9.273220333737331,
      18.563958286272367]
    
    for k in range(2, n, 1):
      k_l.append(k)
      print(">> k= {}".format(k) )
      
      if sim:
        ET_sim = sim_nkmix(num_frun, n, k, dist_m)
        print("ET_nkmix_sim= {}".format(ET_sim) )
        ET_nkmix_sim_l.append(ET_sim)
      
      # ET_sim = sim_samplekmix(num_frun, n, k, dist_m)
      # print("ET_samplekmix_sim= {}".format(ET_sim) )
      # ET_samplekmix_sim_l.append(ET_sim)
      
      ET = ET_randmix(n, k, dist_m['mu'] )
      print("ET_samplekmix= {}".format(ET) )
      ET_samplekmix_l.append(ET)
    print("ET_nkmix_sim_l= {}".format(pprint.pformat(ET_nkmix_sim_l) ) )
    plot.plot(k_l, ET_nkmix_sim_l, label=r'Batch mix', color='g', marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.plot(k_l, ET_samplekmix_sim_l, label=r'Sampling Mix', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.plot(k_l, ET_samplekmix_l, label=r'Sampling mix', color='r', marker=next(marker), mew=mew, ms=ms, linestyle=':')
    
    plot.legend(fontsize=13, loc=2)
    plot.xlabel(r'$k$', fontsize=14)
    plot.ylabel(r'Average delay', fontsize=14)
    plot.title(r'$n= {}$, $\lambda= {}$'.format(n, dist_m['mu'] ) )
    
    ax = plot.gca()
    axins = zoomed_inset_axes(ax, 2, loc=6)
    axins.plot(k_l, ET_nkmix_sim_l, label=r'Batch mix', color='g', marker=next(marker), mew=mew, ms=ms, linestyle=':')
    axins.plot(k_l, ET_samplekmix_l, label=r'Sampling mix', color='r', marker=next(marker), mew=mew, ms=ms, linestyle=':')
    x1, x2, y1, y2 = 1, 18.5, -0.5, 2
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    plot.yticks(visible=False)
    plot.xticks(visible=False)
    axins.xaxis.set_visible('False')
    axins.yaxis.set_visible('False')
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
  n = 40
  dist_m = {'dist': 'Exp', 'mu': 1}
  print("n= {}, X ~ {}".format(n, dist_m) )
  plot_ET_vs_k(n, dist_m)
    
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_randmix_n{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_mix()
