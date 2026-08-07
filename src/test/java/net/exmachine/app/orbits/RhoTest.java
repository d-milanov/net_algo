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

class RhoTest {
  private static double rad(double deg) {
    return Math.toRadians(deg);
  }

  private static double[][] sampleOrbits() {
    return new double[][] {
      OrbitTestSupport.orbit(0.5, 0.3, rad(10), rad(20), rad(30)),
      OrbitTestSupport.orbit(0.8, 0.6, rad(45), rad(120), rad(200)),
      OrbitTestSupport.orbit(1.0, 0.9, rad(80), rad(300), rad(10)),
      OrbitTestSupport.orbit(0.3, 0.1, rad(5), rad(350), rad(100)),
    };
  }

  @Test
  void distanceIsZeroForIdenticalOrbits() {
    Rho space = new Rho(sampleOrbits());

    assertEquals(0.0, space.distance(0, 0), 1e-6);
    assertEquals(0.0, space.distance(2, 2), 1e-6);
  }

  @Test
  void distanceIsSymmetric() {
    double[][] elements = sampleOrbits();
    Rho space = new Rho(elements);

    for (int x = 0; x < elements.length; x++)
      for (int y = 0; y < elements.length; y++)
        assertEquals(space.distance(x, y), space.distance(y, x), 1e-12);
  }

  @Test
  void distanceEqualsAbsoluteSqrtSemiLatusRectumDeltaForCircularOrbitsWithSameOrientation() {
    // With e = 0 the eccentricity term vanishes, and with a shared inclination/node the
    // inclination term collapses too, leaving S^2 = (sqrt(p1) - sqrt(p2))^2.
    double[][] elements = {
      OrbitTestSupport.orbit(0.5, 0.0, rad(10), rad(20), rad(30)),
      OrbitTestSupport.orbit(0.8, 0.0, rad(10), rad(20), rad(30)),
    };
    Rho space = new Rho(elements);

    double expected = Math.abs(Math.sqrt(0.5) - Math.sqrt(0.8));

    assertEquals(expected, space.distance(0, 1), 1e-9);
  }

  @Test
  void distanceEps2ReturnsInfinityBelowThresholdAndExactValueAtOrAboveIt() {
    double[][] elements = sampleOrbits();
    Rho space = new Rho(elements);

    double d = space.distance(0, 1);

    assertEquals(Double.POSITIVE_INFINITY, space.distanceEps2(0, 1, d * d * 0.99));
    assertEquals(d, space.distanceEps2(0, 1, d * d * 1.01), 1e-9);
  }

  @Test
  void triangleCoefficientIsOne() {
    Rho space = new Rho(sampleOrbits());

    assertEquals(1.0, space.getTriangleCoefficient());
  }

  @Test
  void satisfiesTriangleInequality() {
    double[][] elements = sampleOrbits();
    Rho space = new Rho(elements);
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
