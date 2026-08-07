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
 * DH-criterion by Jopek.
 * <p>
 * Jopek T.: Remarks on the meteor orbital similarity D-criterion. Icarus 106(2), 603-607 (1993)
 */
public class D_H extends OrbitSpace {
  private final double triangleCoefficient;

  /**
   * Creates a Jopek criterion space for the given orbits and maximal value of an
   * eccentricity {@code E}.
   *
   * @param elements the orbits, one row per element
   * @param E        maximal expected eccentricity, used to derive the relaxed triangle
   *                 inequality constant
   */
  public D_H(double[][] elements, double E) {
    super(elements);
    this.triangleCoefficient = Util.sqrt(2 + 4 * E * E);
  }

  @Override
  public double getTriangleCoefficient() {
    return triangleCoefficient;
  }

  @Override
  public double distanceEps2(int elemA, int elemB, double epsSquared) {
    double[] b1 = elements[elemA];
    double[] b2 = elements[elemB];

    double de = b1[e] - b2[e],
           S2 = de * de;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double rq = (b1[q] - b2[q]) / (b1[q] + b2[q]);

    S2 += rq * rq;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double cosI = DistUtil.cosI(b1, b2),
           sin_2_I_2 = (1.0 - cosI) / 2.0;

    S2 += 4.0 * sin_2_I_2;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double se = b1[e] + b2[e],
           dOM = b1[OM] - b2[OM],
           cos_delta = DistUtil.cos_d_om(b1, b2),
           sin_delta = DistUtil.sin_d_om(b1, b2),
           ksi;

    if (sin_2_I_2 < 1.0 - 1e-5) {
      double t = (1.0 + DistUtil.cos_sum_i(b1, b2)) *
                  (1.0 - DistUtil.cos_d_OM(b1, b2)) / (1.0 - sin_2_I_2);

      if (t < -1e-4)
        throw new RuntimeException("t < 0: " + t + " elements: " + elemA + "; " + elemB);

      ksi = t > 0 ? Util.sqrt(t) / 2.0 : 0.0;
    } else
      ksi = Math.abs(DistUtil.sin_d_OM(b1, b2));

    double ksi2 = ksi * ksi;

    if (b1[i] + b2[i] > Math.PI)
      ksi = -ksi;

    if (dOM < 0.0)
      ksi = -ksi;

    if (Math.abs(dOM) > Math.PI)
      ksi = -ksi;

    double sin_2_P_2 = 0.5 *
            (1.0 - cos_delta * (1.0 - 2.0 * ksi2) + 2.0 * sin_delta * ksi * Util.sqrt(1.0 - ksi2));

    S2 += se * se * sin_2_P_2;

    return S2 < epsSquared ? Util.sqrt(S2) : Double.POSITIVE_INFINITY;
  }
}