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
 * Approximation of the D-criterion by Southworth &amp; Hawkins for the case of small
 * inclinations {@code i1}, {@code i2}.
 * <p>
 * Southworth R., Hawkins G.: Statistics of meteor streams. Smithsonian Contributions
 * to Astrophysics 7, 261 (1963)
 */
public class D2 extends OrbitSpace {
  /**
   * Creates a Southworth &amp; Hawkins approximation criterion space for the given orbits.
   *
   * @param elements the orbits, one row per element
   */
  public D2(double[][] elements) {
    super(elements);
  }

  @Override
  public double getTriangleCoefficient() {
    return Util.sqrt(2);
  }

  @Override
  public double distanceEps2(int elemA, int elemB, double epsSquared) {
    double[] b1 = elements[elemA];
    double[] b2 = elements[elemB];

    double dq = b1[q] - b2[q],
           S2 = dq * dq;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double de = b1[e] - b2[e];

    S2 += de * de;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double se = b1[e] + b2[e],
           cosP = DistUtil.cos_d_OM(b1, b2) * DistUtil.cos_d_om(b1, b2)
                - DistUtil.sin_d_OM(b1, b2) * DistUtil.sin_d_om(b1, b2),
           sin_2_P_2 = (1.0 - cosP) / 2.0;

    S2 += se * se * sin_2_P_2;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double cosI = DistUtil.cosI(b1, b2),
           sin_2_I_2 = (1.0 - cosI) / 2.0;

    S2 += 4.0 * sin_2_I_2;

    return S2 < epsSquared ? Util.sqrt(S2) : Double.POSITIVE_INFINITY;
  }
}