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
 * The rho_2 orbital distance by Kholshevnikov.
 * <p>
 * Kholshevnikov K., Kokhirova G., Babadzhanov P., Khamroev U.: Metrics in the space of orbits
 * and their application to searching for celestial objects of common origin. Monthly Notices
 * of the Royal Astronomical Society 462(2), 2275-2283 (2016), formula (15)
 */
public class Rho extends OrbitSpace {
  /**
   * Creates a Kholshevnikov rho_2 criterion space for the given orbits.
   *
   * @param elements the orbits, one row per element
   */
  public Rho(double[][] elements) {
    super(elements);
  }

  @Override
  public double getTriangleCoefficient() {
    return 1.0;
  }

  @Override
  public double distanceEps2(int elemA, int elemB, double epsSquared) {
    double[] b1 = elements[elemA];
    double[] b2 = elements[elemB];

    double sp1p2 = 2.0 * b1[sp] * b2[sp],
           S2 = b1[p] + b2[p] - sp1p2 * DistUtil.cosI(b1, b2);

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    S2 += b1[e] * b1[e] * b1[p] + b2[e] * b2[e] * b2[p]
            - sp1p2 * b1[e] * b2[e] * DistUtil.cosP(b1, b2);

    return S2 < epsSquared ? Util.sqrt(S2) : Double.POSITIVE_INFINITY;
  }
}