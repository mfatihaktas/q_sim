import inspect

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

def harmonic_sum(n):
  sum_ = 0
  for i in range(1, n+1):
    sum_ += float(1/i)
  return sum_

def harmonic_2_sum(n):
  sum_ = 0
  for i in range(1, n+1):
    sum_ += float(1/(i**2) )
  return sum_

def gen_harmonic_sum(n, k):
  sum_ = 0
  for i in range(1, n+1):
    if (i - k) == 0:
      continue
    sum_ += float(1/(i*(i - k) ) )
  return sum_
