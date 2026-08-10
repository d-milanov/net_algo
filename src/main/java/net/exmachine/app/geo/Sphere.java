// This algorithm calculates all distances which do not exceed a fixed value epsilon,
// for a given set of points and a distance function r(x, y).
// The algorithm works faster than the complete search in the case when the distance
// function r(x, y) complies with the relaxed triangle inequality:
//                          r(x, y) <= K(r(x, z) + r(y, z)),
// for any three data points x, y, z, and a constant K.
//
// Author: Danila Milanov (danila.milanov@gmail.com)
// Year:   2026

package net.exmachine.app.geo;

import net.exmachine.algo.distance.MetricSpace;

/**
 * A {@link MetricSpace} whose elements are points on the unit sphere, represented as rows of a
 * {@code double[][]} array: {@code row[LAT]} is the latitude in radians, range
 * {@code [-pi/2, pi/2]}, and {@code row[LON]} is the longitude in radians, range
 * {@code [-pi, pi]}. The distance between two points is their great-circle (angular) distance,
 * computed with the haversine formula.
 */
public class Sphere implements MetricSpace {
  /** Index of the latitude, in radians, range {@code [-pi/2, pi/2]}. */
  public static final int LAT = 0;
  /** Index of the longitude, in radians, range {@code [-pi, pi]}. */
  public static final int LON = 1;

  private final double[][] elements;

  /**
   * Creates a sphere metric space for the given points.
   *
   * @param elements the points, one row per element, as {@code {latitude, longitude}}
   */
  public Sphere(double[][] elements) {
    this.elements = elements;
  }

  @Override
  public int getElementsNumber() {
    return elements.length;
  }

  /**
   * {@inheritDoc}
   *
   * @return {@code 1.0}, since the great-circle distance is a proper metric
   */
  @Override
  public double getTriangleCoefficient() {
    return 1.0;
  }

  @Override
  public double distanceEps2(int elemA, int elemB, double epsSquared) {
    double[] p1 = elements[elemA];
    double[] p2 = elements[elemB];

    double dLat = 0.5 * (p1[LAT] - p2[LAT]),
           dLon = 0.5 * (p1[LON] - p2[LON]),
           sinDLat = Math.sin(dLat),
           sinDLon = Math.sin(dLon),
           a = sinDLat * sinDLat + Math.cos(p1[LAT]) * Math.cos(p2[LAT]) * sinDLon * sinDLon;

    // a = sin^2(d/2), which is monotonically increasing in the angular distance d over its
    // whole range [0, pi], so a can be compared against the threshold without computing d.
    if (epsSquared < Math.PI * Math.PI) {
      double sinHalfEps = Math.sin(0.5 * Math.sqrt(epsSquared));

      if (a >= sinHalfEps * sinHalfEps)
        return Double.POSITIVE_INFINITY;
    }

    double sqrtA = Math.sqrt(a),
           d = 2.0 * Math.asin(Math.min(sqrtA, 1.0));

    return d * d < epsSquared ? d : Double.POSITIVE_INFINITY;
  }
}