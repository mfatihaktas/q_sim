from simplex_models import ES
from plot_utils import *
from log_utils import *

def plot_ES_wrt_r_t():
  mew, ms = 2, 6
  def plot_(r):
    t_l, ES_l = [], []
    for t in range(1, 11):
      t_l.append(t)
      ES_l.append(ES(r, t, mu=1) )
    plot.plot(t_l, ES_l, label=r'$r = {}$'.format(r), color=next(dark_color_c), marker=next(marker_c), mew=mew, ms=ms, linestyle=':')
  
  plot_(r=1)
  plot_(r=2)
  plot_(r=3)
  plot_(r=4)
  
  plot.legend(prop={'size':12} )
  plot.xlabel(r'Availability $t$', fontsize=12)
  plot.ylabel(r'Average download time', fontsize=12)
  plot.title(r'Locality $r$,  Server service rate $\mu = 1$', fontsize=12)
  fig = plot.gcf()
  fig.set_size_inches(5, 4)
  plot.savefig("plot_ES_wrt_r_t.pdf", bbox_inches='tight')
  fig.clear()
  log(INFO, "done.")

if __name__ == "__main__":
  plot_ES_wrt_r_t()
