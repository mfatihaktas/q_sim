import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')
import matplotlib.pyplot as plot

from sim import *
import sys, pprint, math, numpy, simpy, getopt, itertools

class JS(object):
  def __init__(self, env):
    self.env = env
    
    self.qt_l = []
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    self.qt_l.append(self.env.now - p.entrance_time)

class MultiQ(object):
  def __init__(self, env, n, sching_m, serv, serv_dist_m):
    self.env = env
    self.n = n
    self.sching_m = sching_m
    
    self.sching_t = self.sching_m["t"]
    
    self.js = JS(env)
    self.q_l = [FCFS(i, env, serv, serv_dist_m, out=self.js) for i in range(self.n) ]
    
    self.jid_counter = 0
    
    self.store = simpy.Store(env)
    self.action = env.process(self.run() )
    
    self.roundrobin_c = 0
    
  def __repr__(self):
    return "MultiQ[n= {}]".format(self.n)
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      
      to_i = None
      if self.sching_t == "power-of-x":
        x = self.sching_m["x"]
        i_l = [i for i in range(self.n) ]
        random.shuffle(i_l)
        i_l = i_l[0:x]
        
        l_min = float("Inf")
        for i in i_l:
          l = self.q_l[i].length()
          if l_min > l:
            l_min = l
            to_i = i
      elif self.sching_t == "rand":
        to_i = random.randint(0, self.n-1)
      elif self.sching_t == "idle-then-random":
        l_l = [self.q_l[i].length() for i in range(self.n) ]
        if 0 in l_l:
          to_i = next(i for i, l in enumerate(l_l) if l == 0)
        else:
          to_i = random.randint(0, self.n-1)
      elif self.sching_t == "idle-then-roundrobin":
        l_l = [self.q_l[i].length() for i in range(self.n) ]
        if 0 in l_l:
          to_i = next(i for i, l in enumerate(l_l) if l == 0)
        else:
          to_i = self.roundrobin_c % self.n
          self.roundrobin_c += 1
      self.q_l[to_i].put(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.jid_counter += 1
    p.job_id = self.jid_counter
    return self.store.put(p)

# *********************************  Sim  *********************************** #
def sim_multiq(num_f_run, ar, n, sching_m, serv, serv_dist_m):
  E_T_sum = 0
  for f in range(num_f_run):
    log(WARNING, "ar= {}, n= {}, sching_m= {}, serv_dist_m= {}".format(ar, n, sching_m, serv_dist_m) )
    
    env = simpy.Environment()
    pg = PG(env, "pg", ar)
    multiq = MultiQ(env, n, sching_m, serv, serv_dist_m)
    pg.out = multiq
    pg.init()
    env.run(until=50000)
    
    l = multiq.js.qt_l
    if len(l): E_T_sum += float(sum(l) )/len(l)
  E_T = E_T_sum/num_f_run
  print(">> E_T= {}".format(E_T) )
  return E_T

def plot_multiq():
  n = 10
  serv = "Exp"
  mu = 1
  if serv == "Exp":
    E_rate = mu
    serv_dist_m = {'mu': mu}
    serv_in_latex = "Exp(\mu={})".format(mu)
  ar_ub = n*E_rate # + 1
  log(WARNING, "n= {}, ar_ub= {}, serv= {}, serv_dist_m= {}".format(n, ar_ub, serv, serv_dist_m) )
  
  E_T_random_l, E_T_powerof2_l, E_T_powerofn_l = [], [], []
  E_T_idlethenrand_l, E_T_idlethenrr_l = [], []
  num_f_run = 1
  sim = False
  if serv == "Exp":
    if n == 100:
      pass
    else:
      sim = True
  
  ar_l = []
  for ar in numpy.arange(1, ar_ub, 1):
  # for ar in numpy.linspace(9, 9.9, 6):
    ar_l.append(ar)
    
    # sching_m = {"t": "rand"}
    # E_T = sim_multiq(num_f_run, ar, n, sching_m, serv, serv_dist_m)
    # E_T_random_l.append(E_T)
    
    sching_m = {"t": "power-of-x", "x": 2}
    E_T = sim_multiq(num_f_run, ar, n, sching_m, serv, serv_dist_m)
    E_T_powerof2_l.append(E_T)
    
    sching_m = {"t": "power-of-x", "x": n}
    E_T = sim_multiq(num_f_run, ar, n, sching_m, serv, serv_dist_m)
    E_T_powerofn_l.append(E_T)
    
    sching_m = {"t": "idle-then-random"}
    E_T = sim_multiq(num_f_run, ar, n, sching_m, serv, serv_dist_m)
    E_T_idlethenrand_l.append(E_T)
    
    sching_m = {"t": "idle-then-roundrobin"}
    E_T = sim_multiq(num_f_run, ar, n, sching_m, serv, serv_dist_m)
    E_T_idlethenrr_l.append(E_T)
    
  # print("E_T_random_l= \n{}".format(pprint.pformat(E_T_random_l) ) )
  # plot.plot(ar_l, E_T_random_l, label="Random", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  print("E_T_powerof2_l= \n{}".format(pprint.pformat(E_T_powerof2_l) ) )
  plot.plot(ar_l, E_T_powerof2_l, label="Power-of-2", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  print("E_T_powerofn_l= \n{}".format(pprint.pformat(E_T_powerofn_l) ) )
  plot.plot(ar_l, E_T_powerofn_l, label="Power-of-n", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  print("E_T_idlethenrand_l= \n{}".format(pprint.pformat(E_T_idlethenrand_l) ) )
  plot.plot(ar_l, E_T_idlethenrand_l, label="Idle-then-random", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  print("E_T_idlethenrr_l= \n{}".format(pprint.pformat(E_T_idlethenrr_l) ) )
  plot.plot(ar_l, E_T_idlethenrr_l, label="Idle-then-roundrobin", marker=next(marker), color=next(dark_color), linestyle=':', mew=mew, ms=ms)
  
  plot.legend(prop={'size':11} )
  plot.xlabel(r'$\lambda$', fontsize=12)
  plot.ylabel(r'$E[T]$', fontsize=12)
  plot.title(r'$S \sim {}$, $n={}$'.format(serv_in_latex, n) )
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]/1.2, def_size[1]/1.2)
  fig.tight_layout()
  plot.savefig("plot_multiq_n_{}.pdf".format(n) )
  log(WARNING, "done; n= {}, serv_dist_m= {}".format(n, serv_dist_m) )
  
if __name__ == "__main__":
  plot_multiq()
  