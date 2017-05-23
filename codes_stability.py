import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy, scipy, cvxpy
# from cvxpy import *
from scipy.spatial import ConvexHull

def random_2d_convex_hull():
  points = numpy.random.rand(30, 2)
  hull = ConvexHull(points)
  
  plot.plot(points[:,0], points[:,1], 'o')
  for simplex in hull.simplices:
    plot.plot(points[simplex, 0], points[simplex, 1], 'k-')
  plot.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
  plot.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
  # plot.show()
  plot.savefig("random_2d_convex_hull.png")
  plot.gcf().clear()

def stability_region_mds_4_2():
  points = numpy.array([ [0, 0],
    [2.5, 0], [0, 2.5], [2, 0], [0, 2], [2, 1], [1, 2] ] )
  hull = ConvexHull(points)
  
  for simplex in hull.simplices:
    plot.plot(points[simplex, 0], points[simplex, 1], 'k-')
  plot.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
  plot.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
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
  objective = Minimize(sum_squares(A*x - b))
  constraints = [0 <= x, x <= 1]
  prob = Problem(objective, constraints)
  
  # The optimal objective is returned by prob.solve().
  result = prob.solve()
  # The optimal value for x is stored in x.value.
  print("x= {}".format(x.value) )
  # The optimal Lagrange multiplier for a constraint is stored in constraint.dual_value.
  print("constraints[0].dual_value= {}".format(constraints[0].dual_value) )

if __name__ == "__main__":
  # random_2d_convex_hull()
  # stability_region_mds_4_2()
  opt()