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


import net.exmachine.algo.distance.MetricSpace;

/**
 * A {@link MetricSpace} whose elements are orbits, represented as rows of a {@code double[][]}
 * array. Each row holds a fixed set of orbital elements, addressed via the named index
 * constants declared on this class (e.g. {@code row[a]} for the semi-major axis).
 */
public abstract class OrbitSpace implements MetricSpace {
  /** Index of the semi-major axis. */
  protected static final int a = 0;
  /** Index of the perihelion distance. */
  protected static final int q = 1;
  /** Index of the semi-latus rectum. */
  protected static final int p = 2;
  /** Index of the eccentricity. */
  protected static final int e = 3;
  /** Index of the inclination. */
  protected static final int i = 4;
  /** Index of the longitude of the ascending node. */
  protected static final int OM = 5;
  /** Index of the argument of perihelion. */
  protected static final int om = 6;
  /** Index of the sine of the inclination. */
  protected static final int sinI = 7;
  /** Index of the cosine of the inclination. */
  protected static final int cosI = 8;
  /** Index of the sine of the argument of perihelion. */
  protected static final int sinom = 9;
  /** Index of the cosine of the argument of perihelion. */
  protected static final int cosom = 10;
  /** Index of the sine of the longitude of the ascending node. */
  protected static final int sinOM = 11;
  /** Index of the cosine of the longitude of the ascending node. */
  protected static final int cosOM = 12;
  /** Index of sqrt(p) of the orbit. */
  protected static final int sp = 13;

  /**
   * Backing storage for the orbits, one row per element, columns addressed by the index
   * constants declared on this class.
   */
  protected final double[][] elements;

  /**
   * Creates an orbit space backed by the given elements.
   *
   * @param elements the orbits, one row per element
   */
  protected OrbitSpace(double[][] elements) {
    this.elements = elements;
  }

  /**
   * {@inheritDoc}
   *
   * @return the number of orbits, i.e. the number of rows in the backing array
   */
  @Override
  public int getElementsNumber() {
    return elements.length;
  }
}