import random, numpy, time, getopt

def random_dot(num):
  t_s = time.time()
  a_l = [random.randrange(1,101,1) for i in range(num) ]
  b_l = [random.randrange(1,101,1) for i in range(num) ]
  numpy.dot(a_l, b_l)
  t_e = time.time()
  return t_e - t_s

def get_opts(argv):
  opt_map = {}
  try:
    opts, args = getopt.getopt(argv, '', ['num_q='] )
  except getopt.GetoptError:
    log(ERROR, "Unexpected command line arg, expecting: exp.py --num_q=<>")
    sys.exit(1)
  
  for opt, arg in opts:
    opt_map[opt] = arg
  return opt_map

if __name__ == "__main__":
  # for num in numpy.linspace(100000, 100000*1000, 10):
  #   num = int(num)
  #   t = random_dot(num)
  #   print("num= {}, t= {}".format(num, t) )
  
  num = int(0.8*500*1000)
  t = random_dot(num)
  # print("num= {}, t= {}".format(num, t) )
  print("{}".format(t) )
