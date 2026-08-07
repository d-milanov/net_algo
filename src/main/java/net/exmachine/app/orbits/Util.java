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
 * Numeric helper functions that tolerate small out-of-domain rounding errors, common when
 * working with values derived from trigonometric identities on orbital elements.
 */
class Util {
  public static final double PRECISION = 1e-6;

  /**
   * Same as {@link Math#acos(double)}, but clamps arguments that are out of {@code [-1, 1]} by
   * no more than {@link #PRECISION} instead of returning {@code NaN}.
   *
   * @throws RuntimeException if {@code x} is out of {@code [-1, 1]} by more than {@link #PRECISION}
   */
  static double acos(double x) {
    if (x >= 1) {
      x -= 1.0;

      if (x > PRECISION)
        throw new RuntimeException("acos: " + x);

      return 0.0;
    } else if (x <= -1.0) {
      x += 1.0;

      if (x < -PRECISION)
        throw new RuntimeException("acos: " + x);

      return Math.PI;
    }

    return Math.acos(x);
  }

  /**
   * Same as {@link Math#asin(double)}, but clamps arguments that are out of {@code [-1, 1]} by
   * no more than {@link #PRECISION} instead of returning {@code NaN}.
   *
   * @throws RuntimeException if {@code x} is out of {@code [-1, 1]} by more than {@link #PRECISION}
   */
  static double asin(double x) {
    if (x >= 1) {
      x -= 1.0;

      if (x > PRECISION)
        throw new RuntimeException("asin: " + x);

      return Math.PI / 2.0;
    } else if (x <= -1.0) {
      x += 1.0;

      if (x < -PRECISION)
        throw new RuntimeException("asin: " + x);

      return -Math.PI / 2.0;
    }

    return Math.asin(x);
  }

  /**
   * Same as {@link Math#sqrt(double)}, but clamps negative arguments that are within
   * {@link #PRECISION} of zero to {@code 0.0} instead of returning {@code NaN}.
   *
   * @throws RuntimeException if {@code x} is negative by more than {@link #PRECISION}
   */
  static double sqrt(double x) {
    if (x <= 0.0) {
      if (x < -PRECISION)
        throw new RuntimeException("sqrt: " + x);

      return 0.0;
    }

    return Math.sqrt(x);
  }
}