import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy, scipy, cvxpy, pprint, itertools
from scipy.spatial import ConvexHull
from patch import *

def random_2d_convex_hull():
  ps = numpy.random.rand(30, 2)
  hull = ConvexHull(ps)
  
  print("ps= {}".format(ps) )
  
  plot.plot(ps[:,0], ps[:,1], 'o')
  for simplex in hull.simplices:
    plot.plot(ps[simplex, 0], ps[simplex, 1], 'k-')
  plot.plot(ps[hull.vertices,0], ps[hull.vertices,1], 'r--', lw=2)
  plot.plot(ps[hull.vertices[0],0], ps[hull.vertices[0],1], 'ro')
  # plot.show()
  plot.savefig("random_2d_convex_hull.png")
  plot.gcf().clear()

def stability_region_mds_4_2():
  ps = numpy.array([ [0, 0],
    [2.5, 0], [0, 2.5], [2, 0], [0, 2], [2, 1], [1, 2] ] )
  hull = ConvexHull(ps)
  
  for simplex in hull.simplices:
    plot.plot(ps[simplex, 0], ps[simplex, 1], 'k-')
  plot.plot(ps[hull.vertices,0], ps[hull.vertices,1], 'r--', lw=2)
  plot.plot(ps[hull.vertices[0],0], ps[hull.vertices[0],1], 'ro')
  # plot.show()
  plot.savefig("stability_region_mds_4_2.png")
  plot.gcf().clear()

def opt():
  # Problem data.
  m = 30
  n = 20
  numpy.random.seed(1)
  A = numpy.random.randn(m, n)
  b = numpy.random.randn(m)
  
  # Construct the problem.
  x = Variable(n)
  obj = Minimize(sum_squares(A*x - b))
  constraints = [0 <= x, x <= 1]
  prob = Problem(obj, constraints)
  
  # The optimal obj is returned by prob.solve().
  result = prob.solve()
  # The optimal value for x is stored in x.value.
  print("x= {}".format(x.value) )
  # The optimal Lagrange multiplier for a constraint is stored in constraint.dual_value.
  print("constraints[0].dual_value= {}".format(constraints[0].dual_value) )

def plot_hull_of_ps(p_l, fname, title):
  ps = numpy.empty((len(p_l), 2))
  for i,p in enumerate(p_l):
    # ps[i,:] = [[p[0], p[1]]]
    ps[i,0] = p[0]
    ps[i,1] = p[1]
  print("ps= {}".format(ps) )
  plot.plot(ps[:,0], ps[:,1], 'o')
  """
  hull = ConvexHull(ps)
  for simplex in hull.simplices:
    plot.plot(ps[simplex, 0], ps[simplex, 1], 'k-')
  plot.plot(ps[hull.vertices,0], ps[hull.vertices,1], 'r--', lw=2)
  plot.plot(ps[hull.vertices[0],0], ps[hull.vertices[0],1], 'ro')
  """
  
  # plot.show()
  # axes = plot.gca()
  # axes.set_xlim([0, 2] )
  plot.title(title)
  plot.savefig(fname)
  plot.gcf().clear()

"""
# Matrix factorization
alpha = 0.002 # 0.0002
D = numpy.zeros((n,2))
step = 0
while step < 100000:
  step += 1
  for i in range(n):
    for j in range(2):
      if M[i,j] == 1:
        e_ij = M[i,j] - numpy.dot(D[i,:], C[:,j])
        for k in range(2):
          D[i][k] = D[i][k] + alpha*(2*e_ij * C[k][j] )
  # print("step= {}, D=\n{}".format(step, D) )
  e = 0
  for i in range(n):
    for j in range(2):
      e = e + pow(M[i][j] - numpy.dot(D[i,:], C[:,j]), 2)
  if e < 0.1:
    break
print("step= {}, D=\n{}".format(step, D) )

M_ = numpy.dot(D, C)
print("M=\n{}".format(M_) )
"""

def generator_matrix_to_M_C(G):
  n = G.shape[1]
  s_rg_l = [[], []]
  
  for s in range(0, 2):
    for c in range(n):
      m = numpy.column_stack((G[:,s], G[:,c]))
      if numpy.linalg.det(m) == 0: # c is a systematic node for s
        s_rg_l[s].append((c,))
    
    for subset in itertools.combinations(range(n), 2):
      if s in subset:
        continue
      m = numpy.column_stack((G[:,subset[0]], G[:,subset[1]]))
      # print("m= {}".format(m) )
      if numpy.linalg.det(m):
        # print("columns {}, {} are LI".format(os, c) )
        s_rg_l[s].append((subset[0], subset[1]))
  print("s_rg_l= {}".format(pprint.pformat(s_rg_l) ) )
  
  r_0 = len(s_rg_l[0] )
  r_1 = len(s_rg_l[1] )
  r = r_0 + r_1
  # if r != len(s_rg_l[1] ):
  #   log(ERROR, "Code was supposed to be symmetric, but it is not.")
  #   return 1
  x = s_rg_l[0] + s_rg_l[1]
  M = numpy.zeros((n,r))
  for i in range(n):
    for j in range(r):
      if i in x[j]:
        M[i,j] = 1
  print("M= {}".format(M) )
  
  C = numpy.zeros((2 ,r))
  C[0,0:r_0] = 1
  C[1,r_0:r] = 1
  print("C= {}".format(C) )
  
  return r, M, C

def generator_matrix(code, n, k=2):
  if k != 2:
    log(ERROR, "Only for k=2")
    return 1
  if code == 'MDS':
    G = numpy.zeros((2, n))
    for j in range(n):
      if j == 0:
        G[0,j] = 1
        G[1,j] = 0
      elif j == 1:
        G[0,j] = 0
        G[1,j] = 1
      else:
        G[0,j] = j-1
        G[1,j] = 1
    return G
  else:
    log(ERROR, "unexpected code= {}".format(code) )
    return 1

def plot_xy_stability_region(n):
  # Rep(2)
  # G = numpy.matrix([[1,0], [0,1]])
  # G = numpy.matrix([[1,0], [0,1], [1,0], [0,1] ])
  # MDS(3,2)
  # G = numpy.matrix([[1,0], [0,1], [1,1]])
  # MDS(4,2)
  # G = numpy.matrix([[1,0], [0,1], [1,1], [2,1]])
  # MDS(5,2)
  # G = numpy.matrix([[1,0], [0,1], [1,1], [2,1], [3,1]])
  # MDS(6,2)
  # G = numpy.matrix([[1,0], [0,1], [1,1], [2,1], [3,1], [4,1]])
  # MDS(7,2)
  # G = numpy.matrix([[1,0], [0,1], [1,1], [2,1], [3,1], [4,1], [5,1]])
  # Mixed
  # G = numpy.matrix([[1,0], [0,1], [1,1], [2,1], [3,1], [1,0] ])
  # G = numpy.matrix([[1,0], [0,1], [1,1], [2,1], [3,1], [1,1] ])
  # G = numpy.matrix([[1,0], [0,1], [1,1], [1,0], [0,1], [1,1] ])
  # G = G.transpose()
  
  code = 'MDS' # 'Rep' # 'MDS' # 'Mixed'
  G = generator_matrix(code, n)
  
  print("G= {}".format(G) )
  n = G.shape[1]
  r, M, C = generator_matrix_to_M_C(G)
  p_l = []
  # 
  x = cvxpy.Variable(r, 1, name='x')
  for b in numpy.linspace(0, 1, 20):
    # print("b= {}".format(b) )
    
    length = math.sqrt((1-b)**2 + b**2)
    w = numpy.matrix([[(1-b)/length, b/length]] )
    # print("w.shape= {}, w= {}".format(w.shape, w) )
    w_ = w*C
    # print("w_= {}".format(w_) )
    # obj = cvxpy.Maximize(w*(C*x) )
    obj = cvxpy.Maximize(w_*x)
    # print("obj= {}".format(obj) )
    constraints = [M*x == 1, x >= 0] # [M*x <= 1, x >= 0]
    prob = cvxpy.Problem(obj, constraints)
    # print("prob= {}".format(prob) )
    prob.solve()
    print("status= {}".format(prob.status) )
    # print("optimal value= {}".format(prob.value) )
    y = C*(x.value)
    # print("optimal y= {}".format(y) )
    p_l.append((y[0], y[1]) )
  plot_hull_of_ps(p_l, "plot_xy_stability_region_{}_n_{}.png".format(code, n),
                  title='{}, n= {}, k= 2'.format(code, n) )
  log(WARNING, "done, code= {}, n= {}".format(code, n) )

if __name__ == "__main__":
  # random_2d_convex_hull()
  # stability_region_mds_4_2()
  # opt()
  
  plot_xy_stability_region(n=4)
  # for n in range(3, 10):
  #   plot_xy_stability_region(n)
  # plot_xy_stability_region(n=100)
