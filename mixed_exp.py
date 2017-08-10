import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt, itertools

from mixed_sim import *
from mixed_models import *

# ******************************  E[T]  ****************************** #
def sim_E_T_mixednet(num_f_run, n, k, qlambda_l):
  E_T_sum = 0
  for f in range(num_f_run):
    log(WARNING, "n= {}, k= {}, qlambda_l= {}".format(n, k, qlambda_l) )
    env = simpy.Environment()
    
    pg = MixedPG(env, _id="pg", qlambda_l=qlambda_l)
    mn = MixedNet(env, n, k)
    pg.out = mn
    env.run(until=50000)
    
    E_T_sum += mn.E_T()
  return E_T_sum/num_f_run

def sim_mixednet_throughput(num_f_run, n, k, qlambda_l):
  throughput_sum = 0
  for f in range(num_f_run):
    log(WARNING, "n= {}, k= {}, qlambda_l= {}".format(n, k, qlambda_l) )
    env = simpy.Environment()
    
    pg = MixedPG(env, _id="pg", qlambda_l=qlambda_l)
    mn = MixedNet(env, n, k)
    pg.out = mn
    env.run(until=50000)
    
    throughput_sum += mn.throughput()
  throughput = throughput_sum/num_f_run
  log(WARNING, "throughput= {}".format(throughput) )
  return throughput
  
def plot_mixednet_E_T():
  n, k = 10, 3 # 100, 10 #
  
  # Assuming each q is identical
  x_l = []
  E_T_sim_l, E_T_approx_l, E_T_ub_l, E_T_lb_l = [], [], [], []
  throughput_sim_l = []
  
  # pe__E_T_l_m = {}
  # for pe in numpy.linspace(0.1, 1, 10):
  #   pe__E_T_l_m[pe] = []
  E_T_mg1_approx_l = []
  
  sim = False
  # by default: for l in numpy.linspace(0.1, 1, 10)
  if n == 10 and k == 3:
    E_T_sim_l= [
      2.1629958454966474,
      1.0718169256987349,
      0.71192203775192842,
      0.53520743064466814,
      0.42740244796207716,
      0.35797970005238161,
      0.30530781328068618,
      0.26729799440082369,
      0.23885349550723575,
      0.21400923474330563]
    # for l in numpy.linspace(1, 10, 10):
    # E_T_sim_l= [
    #   0.21370991075291515,
    #   0.10688086540466046,
    #   0.071415195972370477,
    #   0.053543987548646997,
    #   0.042804095961577408,
    #   0.035765186080486476,
    #   0.030588360059672991,
    #   0.026801612201940089,
    #   0.023791433498213648,
    #   0.021447199148346113]
  elif n == 10 and k == 6:
    E_T_sim_l= [
      7.3981595931117825,
      3.7629979661159312,
      2.5149118894148765,
      1.8633248032001355,
      1.4943548006806542,
      1.2521735927560544,
      1.0742916736200598,
      0.93855937443735571,
      0.83455156153035315,
      0.75032428941029505]
  elif n == 10 and k == 9:
    E_T_sim_l= [
      46.168055155656702,
      24.720469359143856,
      15.17157115845826,
      10.885192749303764,
      8.8291490916279507,
      7.363940075114817,
      6.2862493758335249,
      5.5168866334520201,
      4.9871904538739669,
      4.5303155260079304]
  else:
    sim = True
  
  num_f_run = 1
  for l in numpy.linspace(0.1, 1, 10):
  # for l in numpy.linspace(1, 10, 10):
    x_l.append(l)
    
    if sim:
      E_T_sim = sim_E_T_mixednet(num_f_run, n, k, qlambda_l=n*[l] )
      E_T_sim_l.append(E_T_sim)
      print("E_T_sim= {}".format(E_T_sim) )
    # E_T_ub_l.append(E_T_mixednet_ub(n, k, l) )
    # E_T_lb_l.append(E_T_mixednet_lb(n, k, l) )
    
    # E_T_approx = E_T_mixednet_approx(n, k, l)
    # E_T_approx_l.append(E_T_approx)
    # print("E_T_approx= {}".format(E_T_approx) )
    
    # for pe,E_T_l in pe__E_T_l_m.items():
    #   E_T_l.append(E_T_mg1_approx(n, k, l, pe) )
    
    pe = p_empty_iteratively(n, k, l)
    E_T_approx = E_T_mg1_approx(n, k, l, pe)
    E_T_mg1_approx_l.append(E_T_approx)
    print("l= {}, pe= {}, E_T_approx= {}".format(l, pe, E_T_approx) )
    
    # throughput_sim_l.append(sim_mixednet_throughput(num_f_run, n, k, qlambda_l=n*[l]) )
  mew, ms = 3, 5
  # plot.plot(x_l, E_T_lb_l, label=r'lower-bound', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # plot.plot(x_l, E_T_ub_l, label=r'upper-bound', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  print("E_T_sim_l=\n {}".format(pprint.pformat(E_T_sim_l) ) )
  plot.plot(x_l, E_T_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # print("E_T_approx_l=\n {}".format(pprint.pformat(E_T_approx_l) ) )
  # plot.plot(x_l, E_T_approx_l, label=r'Approximation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # plot.plot(x_l, throughput_sim_l, label=r'sim', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  # for pe,E_T_l in pe__E_T_l_m.items():
  #   plot.plot(x_l, E_T_l, label=r'Approx,$p_0={0:.2f}$'.format(pe), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.plot(x_l, E_T_mg1_approx_l, label=r'M/G/1 Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.title(r'$n= {}$, $k= {}$'.format(n, k) )
  plot.xlabel(r'$\lambda$', fontsize=12)
  #
  plot.legend()
  plot.ylabel(r'$E[T]$ (s)', fontsize=13)
  # plot.ylabel(r'Throughput', fontsize=13)
  fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mixednet_E_T_n_{}_k_{}.png".format(n, k) )
  log(WARNING, "done; n= {}, k= {}".format(n, k) )

# ****************************** Steady state prob  ******************************* #
def sim_steadystate_prob_mixednet(n, k, qlambda_l):
  log(WARNING, "n= {}, k= {}, qlambda_l= {}".format(n, k, qlambda_l) )
  env = simpy.Environment()
  
  pg = MixedPG(env, _id="pg", qlambda_l=qlambda_l)
  mn = MixedNet(env, n, k)
  monitor = MixedNetMonitor(env, mn, 0.05)
  pg.out = mn
  env.run(until=50000)
  
  # for qid, state_counter_map in monitor.qid__state_counter_map_map.items():
  #   total_c = sum([c for s,c in state_counter_map.items() ] )
  #   monitor.qid__state_counter_map_map[qid] = {s:c/total_c for s,c in state_counter_map.items() }
  # print("qid__state_counter_map_map=\n {}".format(pprint.pformat(monitor.qid__state_counter_map_map) ) )
  
  sstate_prob_map = monitor.steadystate_prob_map()
  # log(WARNING, "sstate_prob_map= {}".format(pprint.pformat(sstate_prob_map) ) )
  return sstate_prob_map

def plot_mixednet_steadystate():
  n, k = 100, 3
  
  k_l, Pr_busy_sim_l, Pr_busy_l = [], [], []
  
  # for k in range(1, n+1, 2):
  for k in range(1, n+1, 20):
    k_l.append(k)
    for l in numpy.linspace(1, 1, 1):
      qlambda_l=n*[l]
      steadystate_prob_map = sim_steadystate_prob_mixednet(n, k, qlambda_l)
      print("k= {}, l= {}, steadystate_prob_map= \n{}".format(k, l, pprint.pformat(steadystate_prob_map) ) )
      Pr_busy_sim = 1 - steadystate_prob_map[0]
      print("Pr_busy_sim= {}, Pr_busy= {}".format(Pr_busy_sim, Pr_busy(n, k) ) )
      Pr_busy_sim_l.append(Pr_busy_sim)
      Pr_busy_l.append(Pr_busy(n, k) )
  mew, ms = 3, 5
  plot.plot(k_l, Pr_busy_sim_l, label=r'Simulation', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  plot.plot(k_l, Pr_busy_l, label=r'Approx', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.legend()
  plot.xlabel(r'$k$', fontsize=13)
  plot.ylabel(r'Steady-state probability of busy', fontsize=13)
  fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mixednet_steadystate_n_{}.pdf".format(n) )
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_mixednet_E_T()
  # plot_mixednet_steadystate()
  # Pr_busy_mixednet_approx()
  