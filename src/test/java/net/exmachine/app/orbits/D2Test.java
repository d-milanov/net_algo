// This algorithm calculates all distances which do not exceed a fixed value epsilon,
// for a given set of points and a distance function r(x, y).
// The algorithm works faster than the complete search in the case when the distance
// function r(x, y) complies with the relaxed triangle inequality:
//                          r(x, y) <= K(r(x, z) + r(y, z)),
// for any three data points x, y, z, and a constant K.
//
// Author: Danila Milanov (danila.milanov@gmail.com)
// Year:   2026

package net.exmachine.app.orbits;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class D2Test {
  private static double rad(double deg) {
    return Math.toRadians(deg);
  }

  private static double[][] sampleOrbits() {
    return new double[][] {
      OrbitTestSupport.orbit(0.5, 0.3, rad(5), rad(20), rad(30)),
      OrbitTestSupport.orbit(0.8, 0.6, rad(8), rad(120), rad(200)),
      OrbitTestSupport.orbit(1.0, 0.9, rad(2), rad(300), rad(10)),
      OrbitTestSupport.orbit(0.3, 0.1, rad(6), rad(350), rad(100)),
    };
  }

  @Test
  void distanceIsZeroForIdenticalOrbits() {
    D2 space = new D2(sampleOrbits());

    assertEquals(0.0, space.distance(0, 0), 1e-6);
    assertEquals(0.0, space.distance(2, 2), 1e-6);
  }

  @Test
  void distanceIsSymmetric() {
    double[][] elements = sampleOrbits();
    D2 space = new D2(elements);

    for (int x = 0; x < elements.length; x++)
      for (int y = 0; y < elements.length; y++)
        assertEquals(space.distance(x, y), space.distance(y, x), 1e-12);
  }

  @Test
  void distanceEqualsAbsolutePerihelionDeltaWhenOnlyPerihelionDistanceDiffers() {
    double[][] elements = {
      OrbitTestSupport.orbit(0.5, 0.3, rad(5), rad(20), rad(30)),
      OrbitTestSupport.orbit(0.8, 0.3, rad(5), rad(20), rad(30)),
    };
    D2 space = new D2(elements);

    assertEquals(0.3, space.distance(0, 1), 1e-9);
  }

  @Test
  void distanceEps2ReturnsInfinityBelowThresholdAndExactValueAtOrAboveIt() {
    double[][] elements = sampleOrbits();
    D2 space = new D2(elements);

    double d = space.distance(0, 1);

    assertEquals(Double.POSITIVE_INFINITY, space.distanceEps2(0, 1, d * d * 0.99));
    assertEquals(d, space.distanceEps2(0, 1, d * d * 1.01), 1e-9);
  }

  @Test
  void triangleCoefficientIsSqrtTwo() {
    D2 space = new D2(sampleOrbits());

    assertEquals(Math.sqrt(2), space.getTriangleCoefficient(), 1e-12);
  }

  @Test
  void satisfiesRelaxedTriangleInequality() {
    double[][] elements = sampleOrbits();
    D2 space = new D2(elements);
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
