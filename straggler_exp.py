import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot

from arepeat_models import *

def E_C_workload(task_t, task_dist_m, T, rv_k, J, c=0, n_k=0, r=0, red_top_percent=0):
  w_cancel = True
  d = 0
  
  red_k_min = 0
  if red_top_percent:
    red_k_min = rv_k.inv_cdf(1-red_top_percent)
    print("red_k_min= {}".format(red_k_min) )
  
  E_C = 0
  for k in range(1, T+1):
    l = k
    num_jobs = J*rv_k.pdf(k)
    if (c or n_k or r) and k >= red_k_min:
      if c:
        E_C += num_jobs*E_C_k_c(task_t, task_dist_m, d, k, c, w_cancel)
      elif n_k:
        n = k + n_k
        E_C += num_jobs*E_C_k_l_n(task_t, task_dist_m, d, k, l, n, w_cancel)
      elif r:
        n = math.floor(k*(1+r) )
        # print("n= {}".format(n) )
        E_C += num_jobs*E_C_k_l_n(task_t, task_dist_m, d, k, l, n, w_cancel)
    else:
      E_C += num_jobs*E_C_k_c(task_t, task_dist_m, d, k, c, w_cancel)
  return E_C
  
def E_eload_workload(T, rv_k, c=0, n_k=0):
  E_eload = 0
  for k in range(1, T+1):
    if c:
      E_eload += rv_k.pdf(k) * c
    elif n_k:
      E_eload += rv_k.pdf(k) * ((k+n_k)/k - 1)
  return E_eload

def plot_workload():
  J = 1000 # Total # of jobs
  T = 100 # Max # of tasks
  rv_k = DiscreteUniform(1, T) # BoundedZipf(1, T)
  
  task_t = "Pareto"
  loc, a = 3, 2
  # task_dist_m = {"loc": loc, "a": a}
  task_t_in_latex = r'X \sim Pareto(\lambda={}, \alpha)'.format(loc)
  
  def plot_cost(c=0, n_k=0, r=0, red_top_percent=0):
    a_l, E_C_l = [], []
    for a in numpy.linspace(1.05, 2.5, 30):
      a_l.append(a)
      
      task_dist_m = {"loc": loc, "a": a}
      E_C = E_C_workload(task_t, task_dist_m, T, rv_k, J, c, n_k, r, red_top_percent)
      E_C_l.append(E_C)
    if c: label = "c= {}".format(c)
    elif n_k: label = "n-k= {}".format(n_k)
    elif r: label = "r= {}".format(r)
    else: label = "No red"
    if red_top_percent: label += ", top %{} red".format(100*red_top_percent)
    plot.plot(a_l, E_C_l, label=label, color=next(dark_color), marker=next(marker), ms=10, mew=2, linestyle=':')
  def plot_extra_load(c=0, n_k=0, red_top_percent=0):
    a_l, load_l = [], []
    for a in numpy.linspace(1.05, 2.5, 30):
      a_l.append(a)
      load_l.append(E_eload_workload(T, rv_k, c, n_k) )
      
    if c: label = "c= {}".format(c)
    elif n_k: label = "n-k= {}".format(n_k)
    if red_top_percent: label += ", top %{} red".format(100*red_top_percent)
    plot.plot(a_l, load_l, label=label, color=next(dark_color), marker=next(marker), ms=10, mew=2, linestyle=':')
  plot_cost()
  plot_cost(c=1)
  plot_cost(n_k=1)
  plot_cost(r=0.1)
  # plot_cost(n_k=1, red_top_percent=0.1)
  
  # plot_extra_load(c=1)
  # plot_extra_load(n_k=1)
  
  plot.legend()
  plot.title(r'${}$'.format(task_t_in_latex) + '\n' + r'$J= {}$, $K \sim {}$'.format(J, rv_k), fontsize=12)
  plot.xlabel(r'$\alpha$', fontsize=12)
  plot.ylabel(r'Expected Workload Cost', fontsize=12)
  fig = plot.gcf()
  plot.savefig("plot_workload.pdf", bbox_inches='tight', dpi=fig.dpi)
  plot.gcf().clear()
  log(WARNING, "done; J= {}, T= {}".format(J, T) )

def plot_tail():
  """
  lb, ub = 1, 100
  rv = BoundedZipf(lb, ub)
  v_l, p_l = [], []
  for v in range(lb, ub+1):
    v_l.append(v)
    p_l.append(rv.pdf(v) )
  plot.plot(v_l, p_l, label='{}'.format(rv), color=next(dark_color), marker=next(marker), ms=10, mew=2, linestyle=':')
  plot.legend()
  plot.xlabel(r'$x$', fontsize=12)
  plot.ylabel(r'$Pr\{X = x\}$', fontsize=12)
  """
  rv_l = []
  # rv_l.append(NegBinomial(10, 0.5) )
  # rv_l.append(NegBinomial(10, 0.1) )
  # rv_l.append(NegBinomial(10, 0.01) )
  # rv_l.append(TPareto(l=3, u=1000, a=0.1) )
  # rv_l.append(Exp(1) )
  rv_l.append(Gamma(num_exp=100, rate=0.1) )
  # rv = Pareto(3, 2)
  for rv in rv_l:
    x_l, tail_l = [], []
    for x in range(1000):
      x_l.append(x)
      tail_l.append(rv.tail(x) )
    plot.plot(x_l, tail_l, label='{}'.format(rv), color=next(dark_color), marker=next(marker), linestyle=':')
  
  plot.legend()
  # plot.xscale('log')
  plot.yscale('log')
  plot.xlabel(r'$x$', fontsize=12)
  plot.ylabel(r'$Pr\{X > x\}$', fontsize=12)
  fig = plot.gcf()
  plot.savefig("plot_tail.png", bbox_inches='tight', dpi=fig.dpi)
  plot.gcf().clear()
  log(WARNING, "done; rv= {}".format(rv) )

if __name__ == "__main__":
  # plot_workload()
  # plot_zipf()
  plot_tail()
  
  