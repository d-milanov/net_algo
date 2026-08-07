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


import java.util.*;
import java.util.stream.IntStream;

/**
 * The algorithm for the fast construction of a distance matrix.
 * <p>
 * Given a {@link MetricSpace} of elements and a distance threshold {@code eps}, computes all
 * pairwise distances strictly less than {@code eps}. Instead of comparing every pair of
 * elements (an {@code O(n^2)} operation), the algorithm selects a set of representative
 * "pivot" elements and uses the relaxed triangle inequality guaranteed by the metric space
 * to prune the search: an element only needs to be compared against elements that share a
 * pivot neighbourhood, plus a small fallback set of elements not covered by any pivot.
 */
public class NetAlgo {
  private final MetricSpace space;
  private final double eps;
  private final double delta;
  private final int maxPivots;
  private final int pivotStep;
  private final NetAlgoListener listener;

  private final int[][] adjacent;
  private final ArrayList<NavigableSet<Integer>> neighbourhoods = new ArrayList<>();
  private final BitSet heap;

  /**
   * Computes the distance matrix of all element pairs whose distance is strictly less than
   * {@code eps}.
   *
   * @param space     the metric space to search; must comply with the relaxed triangle
   *                  inequality reported by {@link MetricSpace#getTriangleCoefficient()}
   * @param eps       the distance threshold; only pairs with distance {@code < eps} are
   *                  included in the result
   * @param delta     the neighbourhood radius used to group elements around pivots; smaller
   *                  values yield smaller, more localized neighbourhoods
   * @param maxPivots the maximum number of pivot elements to select
   * @param pivotStep the sampling step used when scanning the space for candidate pivots;
   *                  {@code 1} scans every element, larger values sample sparsely for speed
   * @param listener  an optional listener notified of algorithm progress, or {@code null} if
   *                  no progress reporting is needed
   * @return a {@link DistanceMatrix} containing all element pairs found within {@code eps}
   */
  public static DistanceMatrix calculate(MetricSpace space, double eps, double delta, int maxPivots, int pivotStep, NetAlgoListener listener) {
    return calculate(space, eps, delta, maxPivots, pivotStep, new DistanceMatrix(), listener);
  }

  /**
   * Same as {@link #calculate(MetricSpace, double, double, int, int, NetAlgoListener)}, but
   * records the found pairs into the given {@code matrix} instead of a freshly allocated one
   * (e.g. a {@link DistanceMatrix.Stub} to measure the algorithm without matrix storage
   * overhead).
   *
   * @param matrix the matrix to record found pairs into
   * @return {@code matrix}, for convenience
   */
  public static DistanceMatrix calculate(MetricSpace space, double eps, double delta, int maxPivots, int pivotStep, DistanceMatrix matrix, NetAlgoListener listener) {
    return new NetAlgo(space, eps, delta, maxPivots, pivotStep, listener).calculate(matrix);
  }

  private NetAlgo(MetricSpace space, double eps, double delta, int maxPivots, int pivotStep, NetAlgoListener listener) {
    this.space = space;
    this.eps = eps;
    this.delta = delta;
    this.maxPivots = maxPivots;
    this.pivotStep = pivotStep;
    this.adjacent = new int[space.getElementsNumber()][];
    this.heap = new BitSet(space.getElementsNumber());
    this.listener = listener;
  }


  /**
   * Runs the full algorithm: selects pivots, builds the pivot adjacency structure, then
   * iterates every element to collect its neighbours within {@code eps}.
   */
  private DistanceMatrix calculate(DistanceMatrix result) {
    if (listener != null) listener.onEvent("calculation started");

    int[] pivots = Pivots.selectPivots(space, maxPivots, pivotStep, delta);
    fillAdjacent(pivots);

    streamAllElements().forEach(elem -> {
      iterateNeighbours(elem, result);
      if (listener != null) listener.onEvent("element processed");
    });

    if (listener != null) listener.onEvent("calculation ended");

    return result;
  }

  /**
   * For every element, records pivots which s are "adjacent" to it (within
   * {@code K * (delta + eps)}, close enough to potentially matter for the {@code eps} search)
   * and builds the pivot's {@code delta}-neighbourhood of closely bound elements. Elements
   * that fall in a pivot's neighbourhood are marked in {@link #heap}; elements never covered
   * by any neighbourhood are left clear and handled by direct comparison in
   * {@link #iterateNeighbours}.
   */
  private void fillAdjacent(int[] pivots) {
    if (listener != null) listener.onEvent("fillAdjacent started");

    final int len1 = pivots.length + 1;
    final int len2 = pivots.length + 2;
    final double th = space.getTriangleCoefficient() * (delta + eps);
    final double th2 = th * th;

    int nCnt = 0;

    for (int p : pivots) {
      TreeSet<Integer> neighbourhood = new TreeSet<>();

      for (int elem = 0, len = space.getElementsNumber(); elem < len; elem++) {
        int[] adj = adjacent[elem];
        int ind = 0;

        if (adj == null) {
          adj = adjacent[elem] = new int[len2];
          adj[0] = -1;
        } else while (adj[ind] >= 0) ind++;

        double d = space.distanceEps2(p, elem, th2);

        if (d < th) {
          adj[ind] = nCnt;
          adj[ind + 1] = -1;

          if (d < delta && adj[len1] == 0) {
            neighbourhood.add(elem);
            heap.set(elem);
            adj[len1] = 1;
          }
        }
      }

      nCnt++;
      neighbourhoods.add(neighbourhood);

      if (listener != null) listener.info("Pivot " + p + ": " + neighbourhood.size());
    }

    if (listener != null) {
      listener.onEvent("fillAdjacent ended");
      listener.info("Heap: " + (space.getElementsNumber() - heap.cardinality()));
    }
  }


  /**
   * Finds all not-yet-visited elements within {@code eps} of {@code elem} and adds them to
   * {@code matrix}. Candidates come from two sources: elements sharing a pivot neighbourhood
   * with {@code elem} (via {@link #adjacent} / {@link #neighbourhoods}), and elements not
   * covered by any pivot neighbourhood (via {@link #heap}), which are compared directly.
   */
  private void iterateNeighbours(int elem, DistanceMatrix matrix) {
    final double eps2 = eps * eps;
    int[] adj = adjacent[elem];

    if(adj != null) { // adj is null when no pivots were constructed. Fall back to brute force in this case 
      for (int i : adj) {
        if (i < 0) break;

        neighbourhoods.get(i).tailSet(elem, false).forEach(neighbor -> addToMatrix(elem, neighbor, eps2, matrix));
      }
    }

    for (int hElem = heap.nextClearBit(elem + 1), len = space.getElementsNumber(); hElem < len; hElem = heap.nextClearBit(hElem + 1)) {
      addToMatrix(elem, hElem, eps2, matrix);
    }

  }

  private IntStream streamAllElements() {
    return IntStream.range(0, space.getElementsNumber());
  }

  private void addToMatrix(int elemA, int elemB, double eps2, DistanceMatrix matrix) {
    double d = space.distanceEps2(elemA, elemB, eps2);
    if (d < eps) matrix.add(elemA, elemB, d);
  }
}