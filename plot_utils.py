import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import itertools

NICE_BLUE = '#66b3ff'
NICE_RED = '#ff9999'
NICE_GREEN = '#99ff99'
NICE_ORANGE = '#ffcc99'

nice_color_c = itertools.cycle((NICE_BLUE, NICE_RED, NICE_ORANGE, NICE_GREEN))
dark_color_c = itertools.cycle(('green', 'purple', 'blue', 'magenta', 'purple', 'gray', 'brown', 'turquoise', 'gold', 'olive', 'silver', 'rosybrown', 'plum', 'goldenrod', 'lightsteelblue', 'lightpink', 'orange', 'darkgray', 'orangered'))
# dark_color_c = itertools.cycle(('magenta', 'purple', 'gray', 'brown', 'turquoise', 'gold', 'olive', 'silver', 'rosybrown', 'plum', 'goldenrod', 'lightsteelblue', 'lightpink', 'orange', 'darkgray', 'orangered'))
light_color_c = itertools.cycle(('silver', 'rosybrown', 'plum', 'lightsteelblue', 'lightpink', 'orange', 'turquoise'))
linestyle_c = itertools.cycle(('-', '--', ':', '-.') )
marker_c = itertools.cycle(('o', 'v', '^', 'p', 'd', '<', '>', 'h', 'H', '*', 's', '1' , '2', '3', '4') )
skinny_marker_l = ['x', '+', '1', '2', '3', '4']

mew, ms = 1, 2 # 3, 5

def prettify(ax):
  plot.tick_params(top='off', right='off', which='both')
  ax.patch.set_alpha(0.2)
  ax.spines['right'].set_visible(False)
  ax.spines['top'].set_visible(False)

def plot_points(x_y_l, fname):
  x_l, y_l = [], []
  for x_y in x_y_l:
    x_l.append(x_y[0] )
    y_l.append(x_y[1] )
  
  plot.plot(x_l, y_l, color=NICE_BLUE, marker='o', ls='None')
  plot.savefig('{}.png'.format(fname), bbox_inches='tight')
  plot.gcf().clear()
  log(INFO, "done.")