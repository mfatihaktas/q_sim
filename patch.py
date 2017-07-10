import inspect, math, mpmath, scipy, itertools

dark_color = itertools.cycle(('green', 'red', 'goldenrod', 'blue', 'purple', 'gray', 'brown', 'magenta', 'goldenrod', 'gold', 'olive', 'orangered', 'silver', 'rosybrown', 'plum', 'lightsteelblue', 'lightpink', 'orange', 'turquoise', 'darkgray'))
light_color = itertools.cycle(('silver', 'rosybrown', 'plum', 'lightsteelblue', 'lightpink', 'orange', 'turquoise'))
linestyle = itertools.cycle(('-', '--', '-.', ':') )
marker = itertools.cycle(('^', 'p', '+', 'x', '*', 'v', '<', '>', 'd', '1' , '2', '3', '4') )
skinny_marker_l = ['x', '+', '1', '2', '3', '4']

INFO = 0
DEBUG = 1
WARNING = 2
ERROR = 3

# DEBUG_LEVEL = INFO
DEBUG_LEVEL = WARNING
# DEBUG_LEVEL = ERROR

debug_level__string_map = {INFO: "INFO", DEBUG: "DEBUG", WARNING: "WARNING", ERROR: "ERROR"}

"""
*log: To have a unified logging which can be refactored easily
"""
def sim_log(dlevel, env, caller, action, affected):
  """
  Parameters
  ----------
  dlevel= int -- debug level
  env= simpy.Environment
  caller= string -- name of the sim component acting
  action= string
  affected= any -- whatever component being acted on/with e.g., packet
  """
  if DEBUG_LEVEL <= dlevel:
    print("{} t: {:.2f}] {} {}\n\t{}".format(debug_level__string_map[dlevel], env.now, caller, action, affected) )

def log(dlevel, log):
  """
  Parameters
  ----------
  dlevel= int -- debug level
  log= string to be logged
  """
  if DEBUG_LEVEL <= dlevel:
    print("{}] {}:: {}".format(debug_level__string_map[dlevel], inspect.stack()[1][3], log) )

def list_to_str(l):
  return ",".join("%s" % e for e in l)

def H_cont(n):
  return mpmath.quad(lambda x: (1-x**n)/(1-x), [0, 1] )

def H(n):
  if n == 0:
    return 0
  sum_ = 0
  for i in range(1, n+1):
    sum_ += float(1/i)
  return sum_

def H_2(n):
  sum_ = 0
  for i in range(1, n+1):
    sum_ += float(1/(i**2) )
  return sum_

def gen_H(n, k):
  sum_ = 0
  for i in range(1, n+1):
    if (i - k) == 0:
      continue
    sum_ += float(1/(i*(i - k) ) )
  return sum_

def binomial(n, k):
  if n == k:
    return 1
  elif k == 1:
    return n
  elif k == 0:
    return 1
  elif k > n:
    return 0
  else:
    return math.factorial(n)/math.factorial(k)/math.factorial(n-k)

# def binomial(x, y):
#   try:
#     binom = factorial(x) // factorial(y) // factorial(x - y)
#   except ValueError:
#     binom = 0
#   return binom

def I(u_l, m, n):
  # return B(m, n, u_l=u_l)/B(m, n)
  return scipy.special.betainc(m, n, u_l)

def B(m, n, u_l=1):
  if u_l == 1:
    return scipy.special.beta(m, n)
  return mpmath.quad(lambda x: x**(m-1) * (1-x)**(n-1), [0, u_l] )
  # else:
  #   return I(u_l, m, n)*B(m, n)

def G(z):
  return scipy.special.gamma(z)
