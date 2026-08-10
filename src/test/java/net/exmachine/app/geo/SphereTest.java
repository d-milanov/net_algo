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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class SphereTest {
  private static double rad(double deg) {
    return Math.toRadians(deg);
  }

  private static double[] point(double latDeg, double lonDeg) {
    return new double[] { rad(latDeg), rad(lonDeg) };
  }

  private static double[][] samplePoints() {
    return new double[][] {
      point(90, 0),      // north pole
      point(-90, 0),     // south pole
      point(0, 0),       // equator, prime meridian
      point(0, 90),      // equator, quarter turn east
      point(0, 180),     // equator, antipodal to (0, 0)
      point(45, -120),   // an arbitrary mid-latitude point
      point(-30, 60),    // an arbitrary southern point
    };
  }

  @Test
  void distanceIsZeroForIdenticalPoints() {
    Sphere space = new Sphere(samplePoints());

    assertEquals(0.0, space.distance(0, 0), 1e-9);
    assertEquals(0.0, space.distance(5, 5), 1e-9);
  }

  @Test
  void distanceIsZeroForCoincidentPointsWithDifferentLongitudeAtThePole() {
    // Longitude is degenerate at the poles: any two rows with latitude = +-pi/2 name the
    // same point regardless of their longitude value.
    double[][] elements = { point(90, 0), point(90, 137) };
    Sphere space = new Sphere(elements);

    assertEquals(0.0, space.distance(0, 1), 1e-9);
  }

  @Test
  void distanceIsSymmetric() {
    double[][] elements = samplePoints();
    Sphere space = new Sphere(elements);

    for (int x = 0; x < elements.length; x++)
      for (int y = 0; y < elements.length; y++)
        assertEquals(space.distance(x, y), space.distance(y, x), 1e-12);
  }

  @Test
  void distanceBetweenPolesIsPi() {
    Sphere space = new Sphere(samplePoints());

    assertEquals(Math.PI, space.distance(0, 1), 1e-9);
  }

  @Test
  void distanceFromPoleToEquatorIsHalfPi() {
    Sphere space = new Sphere(samplePoints());

    assertEquals(Math.PI / 2.0, space.distance(0, 2), 1e-9);
  }

  @Test
  void distanceBetweenQuarterTurnEquatorPointsIsHalfPi() {
    Sphere space = new Sphere(samplePoints());

    assertEquals(Math.PI / 2.0, space.distance(2, 3), 1e-9);
  }

  @Test
  void distanceBetweenAntipodalEquatorPointsIsPi() {
    Sphere space = new Sphere(samplePoints());

    assertEquals(Math.PI, space.distance(2, 4), 1e-9);
  }

  @Test
  void distanceEps2ReturnsInfinityBelowThresholdAndExactValueAtOrAboveIt() {
    double[][] elements = samplePoints();
    Sphere space = new Sphere(elements);

    double d = space.distance(2, 5);

    assertEquals(Double.POSITIVE_INFINITY, space.distanceEps2(2, 5, d * d * 0.99));
    assertEquals(d, space.distanceEps2(2, 5, d * d * 1.01), 1e-9);
  }

  @Test
  void distanceEps2AgreesWithDistanceForAThresholdLargerThanTheMaximumPossibleDistance() {
    double[][] elements = samplePoints();
    Sphere space = new Sphere(elements);

    for (int x = 0; x < elements.length; x++)
      for (int y = 0; y < elements.length; y++)
        assertEquals(space.distance(x, y), space.distanceEps2(x, y, 1e6), 1e-9);
  }

  @Test
  void triangleCoefficientIsOne() {
    Sphere space = new Sphere(samplePoints());

    assertEquals(1.0, space.getTriangleCoefficient());
  }

  @Test
  void satisfiesTriangleInequality() {
    double[][] elements = samplePoints();
    Sphere space = new Sphere(elements);
    double k = space.getTriangleCoefficient();

    for (int x = 0; x < elements.length; x++)
      for (int y = 0; y < elements.length; y++)
        for (int z = 0; z < elements.length; z++) {
          if (x == y || x == z || y == z)
            continue;

          double lhs = space.distance(x, y);
          double rhs = k * (space.distance(x, z) + space.distance(y, z));

          assertTrue(lhs <= rhs + 1e-9,
              "triangle inequality violated for (" + x + "," + y + "," + z + "): " + lhs + " > " + rhs);
        }
  }
}