import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import simpy, numpy

from mds_sim import *
from rvs import *
from patch import *

# ########################  Models  ######################## #
def ar_ub_fj(n, X):
  return float(1/moment_i_n_k(1, n, n, X) )

def E_T_fj(ar, n, X):
  # def max_cdf(x):
  #   return X.cdf(x)**n
  # def max_moment(i):
  #   return mpmath.quad(lambda x: i*x**(i-1) * (1 - max_cdf(x) ), [0, mpmath.inf] )
  # return PK(max_moment(1), max_moment(2), ar)
  return PK(moment_i_n_k(1, n, n, X), moment_i_n_k(2, n, n, X), ar)

# ########################  Sim  ######################## #
def test_fj(num_f_run, ar, n, serv, serv_dist_m):
  E_T_f_sum = 0
  for f in range(num_f_run):
    log(WARNING, "ar= {}, n= {}, serv= {}, serv_dist_m= {}".format(ar, n, serv, serv_dist_m) )
    
    env = simpy.Environment()
    pg = PG(env, "pg", ar)
    q = MDSQ("mdsq", env, n, range(n), serv, serv_dist_m)
    pg.out = q
    pg.init()
    env.run(until=10*10*50000)
    
    l = q.jsink.st_l
    if len(l): E_T_f_sum += float(sum(l) )/len(l)
    
    total_n_wins = sum([n for i, n in q.jsink.qid__num_win_map.items() ] )
    print("pg.n_sent= {}, total_n_wins= {}".format(pg.n_sent, total_n_wins) )
    qid__win_freq_map = {i:float(n)/total_n_wins for i, n in q.jsink.qid__num_win_map.items() }
    print("qid__win_freq_map= {}".format(pprint.pformat(qid__win_freq_map) ) )
  E_T = E_T_f_sum/num_f_run
  print(">> E_T= {}".format(E_T) )
  return E_T

def plot_fj():
  n = 2
  serv = "Pareto" # "TPareto"
  l, u, a = 1, 10**6, 2
  if serv == "TPareto":
    X = TPareto(l, u, a)
    serv_dist_m = {'l':l, 'u':u, 'a':a}
  elif serv == "Pareto":
    X = Pareto(l, a)
    serv_dist_m = {'loc':l, 'a':a}
  ar_ub = ar_ub_fj(n, X)
  log(WARNING, "n= {}, serv= {}, serv_dist_m= {}, ar_ub= {}".format(n, serv, serv_dist_m, ar_ub) )
  
  E_T_l, E_T_sim_l = [], []
  num_f_run = 1
  sim = False
  if serv == "TPareto":
    if n == 22:
      pass
    else:
      sim = True
  elif serv == "Pareto":
    if n == 22:
      E_T_sim_l= [
        3.7875159802925884,
        3.6594505295950768,
        4.223943206950012,
        4.589334674521958,
        6.524796278389641,
        5.64633614293259,
        7.252958280015537,
        8.035109860019876,
        8.463351261567757,
        39.12300569764332,
        11.573032446153153,
        13.929789522860153,
        14.965936063862987,
        20.40743954754556,
        27.105625093446594]
    else:
      sim = True
  
  ar_l = []
  for ar in numpy.linspace(0.05, ar_ub, 15):
    ar_l.append(ar)
    if sim:
      E_T_sim_l.append(test_fj(num_f_run, ar, n, serv, serv_dist_m) )
    E_T_l.append(E_T_fj(ar, n, X) )
  log(WARNING, "E_T_sim_l= {}".format(pprint.pformat(E_T_sim_l) ) )
  plot.plot(ar_l, E_T_sim_l, label=r'sim, n={}'.format(n), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  plot.plot(ar_l, E_T_l, label=r'n={}'.format(n), color=next(dark_color), marker=next(marker), mew=mew, ms=ms, linestyle=':')
  
  plot.legend(prop={'size':12})
  plot.xlabel(r'Arrival rate $\lambda$ (Request/s)', fontsize=12)
  plot.ylabel(r'Average download time (s)', fontsize=12)
  if serv == "TPareto":
    serv_in_latex = r'TPareto(l={}, u={}, a={})'.format(l, u, a)
  elif serv == "Pareto":
    serv_in_latex = r'Pareto(l={}, a={})'.format(l, a)
  plot.title(r'$X \sim {}$, $n= {}$'.format(serv_in_latex, n) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_fj_n_{}.pdf".format(n) )
  fig.clear()
  log(WARNING, "done; n= {}".format(n) )

if __name__ == "__main__":
  plot_fj()
