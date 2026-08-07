// This algorithm calculates all distances which do not exceed a fixed value epsilon,
// for a given set of points and a distance function r(x, y).
// The algorithm works faster than the complete search in the case when the distance
// function r(x, y) complies with the relaxed triangle inequality:
//                          r(x, y) <= K(r(x, z) + r(y, z)),
// for any three data points x, y, z, and a constant K.
//
// Author: Danila Milanov (danila.milanov@gmail.com)
// Year:   2026

package net.exmachine.algo.distance;


/**
 * Represents a finite set of points together with a distance function {@code r} that
 * satisfies:
 * <ul>
 *   <li>non-negativity and identity of indiscernibles:
 *       {@code r(x, y) >= 0}, with {@code r(x, y) = 0} if and only if {@code x = y};</li>
 *   <li>symmetry: {@code r(x, y) = r(y, x)};</li>
 *   <li>the relaxed triangle inequality:
 *       {@code r(x, y) <= K(r(x, z) + r(y, z))} for any three points {@code x},
 *       {@code y}, {@code z}, and a constant {@code K}.</li>
 * </ul>
 * The elements of the space are enumerated with contiguous {@code int} indices
 * starting from {@code 0} up to (exclusive) {@link #getElementsNumber()}.
 */
public interface MetricSpace {
  /** Sentinel value denoting the absence of an element. */
  int NULL_ELEMENT = -1;

  /**
   * Returns the number of elements (points) in this metric space.
   */
  int getElementsNumber();

  /**
   * Returns the constant {@code K} used in the relaxed triangle inequality
   * {@code r(x, y) <= K(r(x, z) + r(y, z))}.
   */
  double getTriangleCoefficient();

  /**
   * Returns the distance between two points if its square is less than
   * {@code epsSquared}, otherwise returns {@link Double#POSITIVE_INFINITY}.
   * Implementations of this method may terminate the computation early as
   * soon as it becomes clear that the squared distance will not be less
   * than {@code epsSquared}, without computing the exact distance.
   *
   * @param elemA index of the first point
   * @param elemB index of the second point
   * @param epsSquared threshold on the squared distance
   * @return the distance between {@code elemA} and {@code elemB}, or
   *         {@link Double#POSITIVE_INFINITY} if its square is not less than
   *         {@code epsSquared}
   */
  double distanceEps2(int elemA, int elemB, double epsSquared);

  /**
   * Returns the distance between two points, with no threshold applied.
   *
   * @param elemA index of the first point
   * @param elemB index of the second point
   * @return the distance between {@code elemA} and {@code elemB}
   */
  default double distance(int elemA, int elemB) {
    return distanceEps2(elemA, elemB, Double.POSITIVE_INFINITY);
  }


}