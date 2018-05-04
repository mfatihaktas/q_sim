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
    
    for k in range(2, n, 1):
      k_l.append(k)
      print(">> k= {}".format(k) )
      
      # ET_sim = sim_samplekmix(num_frun, n, k, dist_m)
      # print("ET_samplekmix_sim= {}".format(ET_sim) )
      # ET_samplekmix_sim_l.append(ET_sim)
      
      ET = ET_randmix(n, k, dist_m['mu'] )
      print("ET_samplekmix= {}".format(ET) )
      ET_samplekmix_l.append(ET)
      
      # ET_approx = ET_newmg1approx_(n, k, dist_m)
      # print("ET_newapprox_= {}".format(ET_approx) )
      # ET_approx_l.append(ET_approx)
      
      # ET_sim = sim_nkmix(num_frun, n, k, dist_m)
      # print("ET_nkmix_sim= {}".format(ET_sim) )
      # ET_nkmix_sim_l.append(ET_sim)
      
      plot.xlabel(r'$E[D]$', fontsize=13)
    plot.plot(k_l, ET_samplekmix_sim_l, label=r'Sample-k Mix Sim', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.plot(k_l, ET_approx_l, label=r'M/G/1 Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    # plot.plot(k_l, ET_nkmix_sim_l, label=r'(n, k) Mix Sim', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
    plot.xlabel(r'$k$', fontsize=12)
    plot.ylabel(r'$E[D]$', fontsize=13)
    plot.title(r'$n= {}$, $X \sim {}$'.format(n, dist_to_latex(dist_m) ) )
  
  n = 10
  dist_m = {'dist': 'Exp', 'mu': 1}
  print("n= {}, X ~ {}".format(n, dist_m) )
  plot_ET_vs_k(n, dist_m)
    
  plot.legend()
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mix_n{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_mix()
