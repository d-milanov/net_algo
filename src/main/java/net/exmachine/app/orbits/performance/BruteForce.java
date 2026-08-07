// This algorithm calculates all distances which do not exceed a fixed value epsilon,
// for a given set of points and a distance function r(x, y).
// The algorithm works faster than the complete search in the case when the distance
// function r(x, y) complies with the relaxed triangle inequality:
//                          r(x, y) <= K(r(x, z) + r(y, z)),
// for any three data points x, y, z, and a constant K.
//
// Author: Danila Milanov (danila.milanov@gmail.com)
// Year:   2026

package net.exmachine.app.orbits.performance;

import net.exmachine.algo.distance.DistanceMatrix;
import net.exmachine.algo.distance.MetricSpace;
import net.exmachine.algo.distance.NetAlgoListener;

/**
 * Complete-search construction of a distance matrix: compares every pair of elements directly,
 * without using the relaxed triangle inequality to prune comparisons. Serves as the correctness
 * and performance baseline for {@link net.exmachine.algo.distance.NetAlgo}.
 */
public final class BruteForce {
  private BruteForce() {
  }

  /**
   * Computes the distance matrix of all element pairs whose distance is strictly less than
   * {@code eps}, by comparing every pair of elements in {@code space}.
   *
   * @param space    the metric space to search
   * @param eps      the distance threshold; only pairs with distance {@code < eps} are included
   *                 in the result
   * @param listener an optional listener notified of algorithm progress, or {@code null} if no
   *                 progress reporting is needed
   * @return a {@link DistanceMatrix} containing all element pairs found within {@code eps}
   */
  public static DistanceMatrix calculate(MetricSpace space, double eps, NetAlgoListener listener) {
    return calculate(space, eps, new DistanceMatrix(), listener);
  }

  /**
   * Same as {@link #calculate(MetricSpace, double, NetAlgoListener)}, but records the found
   * pairs into the given {@code matrix} instead of a freshly allocated one (e.g. a
   * {@link DistanceMatrix.Stub} to measure the algorithm without matrix storage overhead).
   *
   * @param matrix the matrix to record found pairs into
   * @return {@code matrix}, for convenience
   */
  public static DistanceMatrix calculate(MetricSpace space, double eps, DistanceMatrix matrix, NetAlgoListener listener) {
    if (listener != null) listener.onEvent("calculation started");

    double eps2 = eps * eps;
    int n = space.getElementsNumber();

    for (int elemA = 0; elemA < n; elemA++) {
      for (int elemB = elemA + 1; elemB < n; elemB++) {
        double d = space.distanceEps2(elemA, elemB, eps2);

        if (d < eps)
          matrix.add(elemA, elemB, d);
      }

      if (listener != null) listener.onEvent("element processed");
    }

    if (listener != null) listener.onEvent("calculation ended");

    return matrix;
  }
}
