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
 * D-criterion by Drummond.
 * <p>
 * Drummond J.: A test of comet and meteor shower associations. Icarus 45(3), 545-553 (1981)
 */
public class D_D extends OrbitSpace {
  private final double triangleCoefficient;

  /**
   * Creates a Drummond criterion space for the given orbits and maximal value of an
   * eccentricity {@code E}.
   *
   * @param elements the orbits, one row per element
   * @param E        maximal expected eccentricity, used to derive the relaxed triangle
   *                 inequality constant
   */
  public D_D(double[][] elements, double E) {
    super(elements);
    this.triangleCoefficient = Util.sqrt(1 + Math.max(1.0, E * E));
  }

  @Override
  public double getTriangleCoefficient() {
    return triangleCoefficient;
  }

  @Override
  public double distanceEps2(int elemA, int elemB, double epsSquared) {
    double[] b1 = elements[elemA];
    double[] b2 = elements[elemB];

    double rq = (b1[q] - b2[q]) / (b1[q] + b2[q]),
           S2 = rq * rq;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double re = (b1[e] - b2[e]) / (b1[e] + b2[e]);

    S2 += re * re;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double I = Util.acos(DistUtil.cosI(b1, b2)) / Math.PI;

    S2 += I * I;

    if (S2 >= epsSquared)
      return Double.POSITIVE_INFINITY;

    double se_2 = (b1[e] + b2[e]) / 2.0,
           P = Util.acos(DistUtil.cosP(b1, b2)) / Math.PI;

    S2 += se_2 * se_2 * P * P;

    return S2 < epsSquared ? Util.sqrt(S2) : Double.POSITIVE_INFINITY;
  }
}