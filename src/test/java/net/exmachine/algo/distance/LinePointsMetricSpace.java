// This algorithm calculates all distances which do not exceed a fixed value epsilon,
// for a given set of points and a distance function r(x, y).
// The algorithm works faster than the complete search in the case when the distance
// function r(x, y) complies with the relaxed triangle inequality:
//                          r(x, y) <= K(r(x, z) + r(y, z)),
// for any three data points x, y, z, and a constant K.
//
// Author: Danila Milanov (danila.milanov@gmail.com)
// Year:   2026

package net.exmachine.algo.distance;

/**
 * Test double: a {@link MetricSpace} of points on the real line, with the plain (exact,
 * {@code K = 1}) triangle inequality distance {@code |a - b|}.
 */
class LinePointsMetricSpace implements MetricSpace {
  private final double[] points;

  LinePointsMetricSpace(double... points) {
    this.points = points;
  }

  @Override
  public int getElementsNumber() {
    return points.length;
  }

  @Override
  public double getTriangleCoefficient() {
    return 1.0;
  }

  @Override
  public double distanceEps2(int elemA, int elemB, double epsSquared) {
    double d = Math.abs(points[elemA] - points[elemB]);
    return d * d < epsSquared ? d : Double.POSITIVE_INFINITY;
  }
}