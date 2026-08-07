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

import static net.exmachine.app.orbits.OrbitSpace.cosI;
import static net.exmachine.app.orbits.OrbitSpace.cosOM;
import static net.exmachine.app.orbits.OrbitSpace.cosom;
import static net.exmachine.app.orbits.OrbitSpace.sinI;
import static net.exmachine.app.orbits.OrbitSpace.sinOM;
import static net.exmachine.app.orbits.OrbitSpace.sinom;

/**
 * Trigonometric identities on pairs of {@link OrbitSpace} elements (orbit rows), shared by
 * orbit distance metrics such as {@link D1}.
 */
class DistUtil {
  static double cosI(double[] b1, double[] b2) {
    return b1[cosI] * b2[cosI] + b1[sinI] * b2[sinI] * cos_d_OM(b1, b2);
  }

  static double cosP(double[] b1, double[] b2) {
    double s1 = b1[sinI],
           s2 = b2[sinI],
           c1 = b1[cosI],
           c2 = b2[cosI],
           so1 = b1[sinom],
           so2 = b2[sinom],
           co1 = b1[cosom],
           co2 = b2[cosom];

    return s1 * s2 * so1 * so2 +
           (co1 * co2 + c1 * c2 * so1 * so2) * cos_d_OM(b1, b2) +
           (c2 * co1 * so2 - c1 * so1 * co2) * sin_d_OM(b1, b2);
  }

  static double cos_d_OM(double[] b1, double[] b2) {
    return b1[cosOM] * b2[cosOM] + b1[sinOM] * b2[sinOM];
  }

  static double sin_d_OM(double[] b1, double[] b2) {
    return b1[sinOM] * b2[cosOM] - b1[cosOM] * b2[sinOM];
  }

  static double cos_d_om(double[] b1, double[] b2) {
    return b1[cosom] * b2[cosom] + b1[sinom] * b2[sinom];
  }

  static double sin_d_om(double[] b1, double[] b2) {
    return b1[sinom] * b2[cosom] - b1[cosom] * b2[sinom];
  }

  static double cos_sum_i(double[] b1, double[] b2) {
    return b1[cosI] * b2[cosI] - b1[sinI] * b2[sinI];
  }
}