import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.cm as cm # cm.rainbow
import numpy, math, collections
from simplex_sim import *

# def plot_split_merge_approx_error(arr_rate, mu, n):
#   gamma = arr_rate
#   k__i_list_map, k__p_i_list_map = {}, {}
#   k_list, k__avg_active_num_server_list = [], []
#   for k in range(2, n+1):
#     if k not in k__p_i_list_map:
#       k__i_list_map[k], k__p_i_list_map[k] = [], []
    
#     p_k_1 = 1/((n-k+1)*mu/gamma * (harmonic_sum(n) - harmonic_sum(n-k+1) ) + 1)
#     def p_i(i):
#       return (n-k+1)*mu/((n-i)*gamma) * p_k_1 if i < k-1 else p_k_1
#     for i in range(0, k):
#       k__i_list_map[k].append(i)
#       k__p_i_list_map[k].append(p_i(i) )
#     #
#     k_list.append(k)
#     k__avg_active_num_server_list.append((n-k+1)*(k-1)*mu/gamma * p_k_1)
#   # print("k__i_list_map=\n{}".format(pprint.pformat(k__i_list_map) ) )
#   # print("k__p_i_list_map=\n{}".format(pprint.pformat(k__p_i_list_map) ) )
#   plot.figure(1)
#   plot.subplot(211)
#   for k in k__i_list_map:
#     plot.plot(k__i_list_map[k], k__p_i_list_map[k], '-o', label="k={}".format(k), color=plot_color_list[k%len(plot_color_list) ] )
#     # plot.legend()
#   plot.xlabel('i')
#   plot.ylabel(r'$p_i$')
#   plot.title(r'n= {}, $\lambda$= {}, $\mu$= {}'.format(n, arr_rate, mu) )
#   plot.subplot(212)
#   plot.plot(k_list, k__avg_active_num_server_list, 'b-o')
#   # plot.legend()
#   plot.xlabel('k')
#   plot.ylabel('Average active number of servers')
#   # plot.title(r'n= {}, $\lambda$= {}, $\mu$= {}'.format(n, arr_rate, mu) )
#   plot.tight_layout()
#   # ax = plot.gca()
#   # ax.grid(True)
#   plot.savefig("plot_split_merge_approx_error__n_{}.png".format(n) )

def plot_dist(dist=None):
  N = 10000 * 100
  if dist == None:
    data = [random.expovariate(1) for i in range(N) ]
  else:
    data = [dist() for i in range(N) ]
  plot.clf()
  
  # plot.hist(data, bins='auto', normed=True)
  
  # w_l = numpy.ones_like(data)/float(len(data) )
  # plot.hist(data, weights=w_l)
  
  plot.hist(data, bins=numpy.arange(min(data)-0.5, max(data)-0.5, 1), normed=True)
  d_f_m = {d:f/len(data) for d,f in collections.Counter(data).items() }
  print("d_f_m= {}".format(d_f_m) )
  
  plot.savefig("plot_dist.png")

def plot_binomial_dist__approx():
  n = 15
  p = 0.4
  def comp_dist(k):
    if k > n:
      return 0
    sum_ = 0
    for i in range(k, n+1):
      sum_ += binomial(n, i) * p**i * (1-p)**(n-i)
    return sum_
  def dist(k):
    if k > n:
      return 0
    sum_ = 0
    for i in range(0, k+1):
      sum_ += binomial(n, i) * p**i * (1-p)**(n-i)
    return sum_
  def chernoff_bound_on_upper_tail(k):
    p_ = k/n
    print("p_= {}".format(p_) )
    return math.exp(n*((p_*math.log(p/p_) ) + (1-p_)*math.log((1-p)/(1-p_) ) ) )
  
  k_l, dist_l, approx_dist_l = [], [], []
  # for k in range(0, n+2):
  for k in range(int(p*n), n):
    k_l.append(k)
    dist_l.append(comp_dist(k) )
    approx_dist_l.append(chernoff_bound_on_upper_tail(k) )
  marker = itertools.cycle(('^', 'p', 'x', '+', '*', 'v', 'o') )
  # print("k_l= {}".format(pprint.pformat(k_l) ) )
  plot.plot(k_l, approx_dist_l, color='red', label=r'$Pr\{N \leq k\}_{UB}$', marker=next(marker), linestyle='', mew=2)
  plot.plot(k_l, dist_l, color='black', label=r'$Pr\{N \leq k\}$', marker=next(marker), linestyle='', mew=2)
  plot.xlabel(r'$k$')
  plot.ylabel("E[T] (s)")
  plot.title(r'n= {}, p= {}'.format(n, p) )
  plot.savefig("plot_binomial_dist__approx_n_{}.png".format(n) )
  plot.gcf().clear()
  log(WARNING, "done; n= {}, p= {}".format(n, p) )

def sum_of_harmonics():
  k = 10
  mu = 1
  
  def E_C(n):
    sum_ = sum([H(n-i) for i in range(1, k+1) ] )
    return 1/mu*(k*H(n) - sum_ + (n-k)*(H(n)-H(n-k) ) )
  
  for n in range(k, 4*k):
    print("n= {}, E_C= {}".format(n, E_C(n) ) )

def plot_deneme():
  # # Checking G(n-r+1)/G(n-r+1 - 1/a) ~= (n-r)**(1/a)
  # n = 100
  # r = 90
  # a_l = []
  # actual_l, approx_l = [], []
  # def exact(a):
  #   return G(n-r+1)/G(n-r+1 - 1/a)
  # def approx(a):
  #   return (n-r)**(1/a)
  # plot.title(r'$n= {}$, $r= {}$'.format(n, r) )
  
  # # Checking R ~ Bin(k,q), E[(n-R)^(1/a)] =~ (n-kq)^(1/a)
  # n = 100
  # k = 5
  # q = 0.9
  # def exact(a):
  #   sum_ = 0
  #   for r in range(k+1):
  #     sum_ += (n-r)**(1/a) * binomial(k, r) * q**r * (1-q)**(k-r)
  #   return sum_
  # def approx(a):
  #   return (n - k*q)**(1/a)
  # plot.title(r'$n= {}$, $k= {}$, $q= {}$'.format(n, k, q) )
  
  # # Checking R ~ Bin(k,q), E[G(k-R+1-1/a)/G(k-R)] =~ G(k-kq+1-1/a)/G(k-kq)
  # k = 10
  # q = 0.2
  # def exact(a):
  #   sum_ = 0
  #   for r in range(k+1):
  #     sum_ += G(k-r+1-1/a)/G(k-r) * binomial(k, r) * q**r * (1-q)**(k-r)
  #   return sum_
  # def approx(a):
  #   return G(k-k*q+1-1/a)/G(k-k*q)
  # plot.title(r'$k= {}$, $q= {}$'.format(k, q) )
  
  # a_l = []
  # actual_l, approx_l = [], []
  # for a in numpy.linspace(1, 100, 2*100):
  #   a_l.append(a)
  #   actual_l.append(exact(a) )
  #   approx_l.append(approx(a) )
  
  # plot.plot(a_l, actual_l, 'bx', label='Actual', zorder=1, ms=8)
  # plot.plot(a_l, approx_l, 'ro', label='Approx', zorder=2, ms=5)
  # plot.legend()
  # plot.xlabel(r'$a$')
  
  # Plotting CDF of Exp or Pareto
  loc, a = 3, 0.3
  X = Pareto(loc, a)
  # D, mu = 3, 0.1
  # X = Exp(mu, D)
  x_l, cdf_l = [], []
  for x in numpy.logspace(0, 6, 100):
    x_l.append(x)
    cdf_l.append(X.cdf(x) )
  plot.xscale('log')
  plot.plot(x_l, cdf_l, 'bx', label='{}'.format(X) )
  plot.legend()
  plot.xlabel(r'$x$')
  plot.ylabel(r'CDF')
  
  plot.savefig("plot_deneme.png")
  plot.gcf().clear()
  log(WARNING, "done.")

if __name__ == "__main__":
  # plot_dist()
  # plot_dist(dist=dolly_slowdown_dist)
  # sum_of_harmonics()
  plot_deneme()
  