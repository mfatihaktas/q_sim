import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot
import sys, pprint, math, numpy, simpy, getopt, itertools

from mixed_sim_components import *
from mixed_models import *

def E_T_mixed_net(num_f_run, n, k, qlambda_l):
  qid__E_T_sum_l = len(qlambda_l)*[0]
  for f in range(num_f_run):
    log(WARNING, "n= {}, k= {}, qlambda_l= {}".format(n, k, qlambda_l) )
    env = simpy.Environment()
    
    pg = MixedPG(env, _id="pg", qlambda_l=qlambda_l)
    mn = MixedNet(env, n, k)
    pg.out = mn
    env.run(until=50000)
    
    for i,q in enumerate(mn.id_q_map):
      qid__E_T_sum_l[i] += q.avg_qtime()
  qid__E_T_l = [e/num_f_run for e in qid__E_T_sum_l]
  log(WARNING, "qid__E_T_sum_l= {}".format(pprint.pformat(qid__E_T_sum_l) ) )
  return qid__E_T_l[0]

def throughput_mixed_net(num_f_run, n, k, qlambda_l):
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
  
def plot_mixed_net():
  n = 4 # 3
  k = 3 # 4
  
  num_f_run = 1
  x_l = []
  """
  # l_0 != l_1 = ... = l_{n-1}
  qid__E_T_sim_l = {i:[] for i in range(n) }
  qid__E_T_l = {i:[] for i in range(n) }
  
  qlambda_l = n*[1]
  for highest_l in numpy.linspace(2, 2.9, 10):
    x_l.append(highest_l)
    
    qlambda_l[0] = highest_l
    i__E_T_l = E_T_mixed_net(n, k, qlambda_l)
    i__E_T_sim_l = E_T_mixed_net(num_f_run, n, k, qlambda_l)
    for i in range(n):
      qid__E_T_l[i].append(i__E_T_l[i] )
      qid__E_T_sim_l[i].append(i__E_T_sim_l[i] )
  
  for i in range(n):
    plot.plot(x_l, qid__E_T_sim_l[i], label=r'sim, $\lambda_{}$'.format(i), color=next(dark_color), marker=next(marker), mew=2, ms=5, linestyle=':')
    plot.plot(x_l, qid__E_T_l[i], label=r'$\lambda_{}$'.format(i), color=next(dark_color), marker=next(marker), mew=2, ms=5, linestyle=':')
  plot.title(r'Other $\lambda=1$')
  plot.xlabel(r'Highest arrival rate $\lambda_0$', fontsize=12)
  """
  # l_0 = l_1 = ... = l_{n-1}
  E_T_ub_l, E_T_sim_l, E_T_lb_l = [], [], []
  
  
  throughput_sim_l = []
  for l in numpy.linspace(1, 10, 10):
    x_l.append(l)
    
    qlambda_l=n*[l]
    qlambda_l[0] += 1
    
    # E_T_ub_l.append(E_T_mixed_net_ub(n, k, l) )
    # E_T_lb_l.append(E_T_mixed_net_lb(n, k, l) )
    # E_T_sim_l.append(E_T_mixed_net(num_f_run, n, k, qlambda_l=n*[l] )[0] )
    # throughput_sim_l.append(throughput_mixed_net(num_f_run, n, k, qlambda_l=n*[l]) )
    throughput_sim_l.append(throughput_mixed_net(num_f_run, n, k, qlambda_l) )
  mew, ms = 3, 5
  # plot.plot(x_l, E_T_ub_l, label=r'upper-bound', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # plot.plot(x_l, E_T_sim_l, label=r'sim', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  # plot.plot(x_l, E_T_lb_l, label=r'lower-bound', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  plot.plot(x_l, throughput_sim_l, label=r'sim', color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.title(r'$\lambda_i=\lambda$')
  plot.xlabel(r'$\lambda$', fontsize=12)
  # 
  plot.legend()
  # plot.ylabel(r'Average queue time $E[T]$ (s)', fontsize=12)
  plot.ylabel(r'Throughput', fontsize=12)
  fig = plot.gcf()
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.tight_layout()
  plot.savefig("plot_mixed_net_n_{}_k_{}.pdf".format(n, k) )
  log(WARNING, "done; n= {}, k= {}".format(n, k) )

if __name__ == "__main__":
  plot_mixed_net()
  