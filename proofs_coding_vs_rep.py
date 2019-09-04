import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy as np

from rvs import *
from patch import *

def test_Twcoding_leq_Twrep():
  def Pr_Twcoding_g_t(k, n, q): # q: Pr{X > x}
    return sum([q**i*(1-q)**(n-i)*binom(n, i) for i in range(n-k+1, n+1) ] )
  
  def Pr_Twrep_g_t(k, c, q):
    return 1 - (1 - q**(c+1) )**k
  
  k = 200
  c = 1
  q_l = []
  Pr_Twcoding_g_t_l, Pr_Twrep_g_t_l = [], []
  for q in np.linspace(0.001, 0.999, 100):
    q_l.append(q)
    Pr_Twcoding_g_t_l.append(Pr_Twcoding_g_t(k, k*(c+1), q) )
    Pr_Twrep_g_t_l.append(Pr_Twrep_g_t(k, c, q) )
  
  plot.plot(q_l, Pr_Twcoding_g_t_l, label='Coding', color='red', marker='x', ls=':', lw=2)
  plot.plot(q_l, Pr_Twrep_g_t_l, label='Replication', color='green', marker='o', ls=':', lw=2)
  
  plot.legend(loc='upper right', fontsize=14, framealpha=0.25)
  plot.xlabel(r'$Pr\{X > t\}$', fontsize=20)
  plot.ylabel(r'$Pr\{T > t\}$', fontsize=20)
  plot.title('k= {}, c= {}'.format(k, c) )
  prettify(plot.gca() )
  plot.gcf().set_size_inches(4, 3)
  plot.savefig("plot_Twcoding_leq_Twrep.pdf", bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done.")

if __name__ == "__main__":
  test_Twcoding_leq_Twrep()
