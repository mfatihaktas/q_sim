from plot_utils import *
from patch import *

def E_S_r_t(mu, r, t):
  return 1/mu * sum([binom(t, i)*(-1)**(i-1) * H(r*i) for i in range(1, t+1) ] )

def E_S_r_t_lb(mu, r, t):
  return 1/(mu*t) * H(r)

def E_S_r_t_ub(mu, r, t):
  return B(t, 1/r)/(mu*r) * H(r)

def Pr_S_r_t_g_S(r, t):
  return 1/r*B(t+1, 1/r)

ms, mew = 10, 0.1
def plot_E_S_r_t():
  mu = 1
  
  def plot_wrt_t(r):
    print("r= {}".format(r) )
    t_l = []
    E_S_r_t_l, E_S_r_t_lb_l, E_S_r_t_ub_l = [], [], []
    
    for t in range(1, 30):
      t_l.append(t)
      E_S_r_t_l.append(E_S_r_t(mu, r, t) )
      E_S_r_t_lb_l.append(E_S_r_t_lb(mu, r, t) )
      E_S_r_t_ub_l.append(E_S_r_t_ub(mu, r, t) )
    
    plot.plot(t_l, E_S_r_t_ub_l, label=r'Upper-bound', color="red", marker='^', ms=ms, mew=mew, ls=':')
    plot.plot(t_l, E_S_r_t_l, label=r'Exact', color="black", marker='x', ms=9, mew=3, ls=':')
    plot.plot(t_l, E_S_r_t_lb_l, label=r'Lower-bound', color="green", marker='v', ms=ms, mew=mew, ls=':')
    
    plot.legend(fontsize=14, framealpha=0.25, numpoints=1) # loc='upper left'
    plot.xlabel(r'$t$', fontsize=20)
    plot.ylabel(r'$E[S_{r, t}]$', fontsize=20)
    plot.title(r'$r= {}$'.format(r), fontsize=20)
    fig = plot.gcf()
    fig.set_size_inches(5, 4)
    fig.tight_layout()
    plot.savefig("plot_E_S_r_t_r{}.pdf".format(r), bbox_inches='tight')
    plot.gcf().clear()
  
  plot_wrt_t(r=2)
  plot_wrt_t(r=4)
  
  log(WARNING, "done.")

def plot_Pr_S_r_t_g_S():
  def plot_(r):
    t_l, Pr_S_r_t_g_S_l = [], []
    for t in range(1, 30):
      t_l.append(t)
      Pr_S_r_t_g_S_l.append(Pr_S_r_t_g_S(r, t) )
    
    plot.plot(t_l, Pr_S_r_t_g_S_l, label=r'$r= {}$'.format(r), color=next(dark_color), marker=next(marker), ms=ms, mew=mew, ls=':')
  
  plot_(r=1)
  plot_(r=2)
  plot_(r=3)
  plot_(r=4)
  
  plot.legend(fontsize=14, framealpha=0.25, numpoints=1) # loc='upper left'
  plot.xlabel(r'$t$', fontsize=20)
  plot.ylabel(r'$\Pr\{S_{r, t} > S_0\}$', fontsize=20)
  fig = plot.gcf()
  fig.set_size_inches(5, 4)
  fig.tight_layout()
  plot.savefig("plot_Pr_S_r_t_g_S.pdf", bbox_inches='tight')
  plot.gcf().clear()

if __name__ == "__main__":
  # plot_E_S_r_t()
  plot_Pr_S_r_t_g_S()
