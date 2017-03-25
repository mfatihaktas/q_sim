import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import numpy, scipy
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
  points = numpy.array([
    [2.5, 0], [0, 2.5], [2, 0], [0, 2], [2, 1], [1, 2] ] )
  hull = ConvexHull(points)
  
  for simplex in hull.simplices:
    plot.plot(points[simplex, 0], points[simplex, 1], 'k-')
  plot.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
  plot.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
  # plot.show()
  plot.savefig("stability_region_mds_4_2.png")
  plot.gcf().clear()

if __name__ == "__main__":
  # random_2d_convex_hull()
  stability_region_mds_4_2()