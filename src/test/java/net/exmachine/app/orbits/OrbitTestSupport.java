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

/**
 * Test helper for building {@link OrbitSpace} rows from the classic orbital elements
 * (perihelion distance, eccentricity, inclination, longitude of the ascending node,
 * argument of perihelion), filling in the derived columns (semi-latus rectum, its square
 * root, and the cached sines/cosines) that the {@code net.exmachine.app.orbits} distance
 * criteria read directly out of the row.
 */
final class OrbitTestSupport {
  private OrbitTestSupport() {
  }

  /**
   * Builds a single orbit row.
   *
   * @param q  perihelion distance
   * @param e  eccentricity
   * @param i  inclination, in radians
   * @param OM longitude of the ascending node, in radians
   * @param om argument of perihelion, in radians
   */
  static double[] orbit(double q, double e, double i, double OM, double om) {
    double[] row = new double[14];

    double p = q * (1.0 + e);

    row[OrbitSpace.a] = e == 1.0 ? Double.POSITIVE_INFINITY : q / (1.0 - e);
    row[OrbitSpace.q] = q;
    row[OrbitSpace.p] = p;
    row[OrbitSpace.e] = e;
    row[OrbitSpace.i] = i;
    row[OrbitSpace.OM] = OM;
    row[OrbitSpace.om] = om;
    row[OrbitSpace.sinI] = Math.sin(i);
    row[OrbitSpace.cosI] = Math.cos(i);
    row[OrbitSpace.sinom] = Math.sin(om);
    row[OrbitSpace.cosom] = Math.cos(om);
    row[OrbitSpace.sinOM] = Math.sin(OM);
    row[OrbitSpace.cosOM] = Math.cos(OM);
    row[OrbitSpace.sp] = Math.sqrt(p);

    return row;
  }
}