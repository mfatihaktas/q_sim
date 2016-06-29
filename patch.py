
INFO = 0
DEBUG = 1
WARNING = 2
ERROR = 3

# DEBUG_LEVEL = INFO
DEBUG_LEVEL = ERROR

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
    print("{}] {}".format(debug_level__string_map[dlevel], log) )
  