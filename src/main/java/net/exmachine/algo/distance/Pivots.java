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
import java.util.stream.Collectors;

import static net.exmachine.algo.distance.MetricSpace.NULL_ELEMENT;

/**
 * Various pivots selection algorithms.
 * <p>
 * Pivots are elements chosen so that their {@code delta}-neighbourhoods together cover as much
 * of the space as possible with as little overlap as possible, which lets {@link NetAlgo} prune
 * most pairwise comparisons when searching for close elements.
 */
public class Pivots {
  /**
   * Selects a set of pivot elements from {@code space} using a greedy set-cover-like strategy:
   * elements are sampled every {@code pivotStep} and grouped by shared {@code delta}-closeness,
   * then repeatedly the element covering the largest number of not-yet-covered elements is
   * picked as a pivot until either {@code maxPivots} pivots have been selected or no elements
   * remain to cover.
   *
   * @param space     the metric space to select pivots from
   * @param maxPivots the maximum number of pivots to select
   * @param pivotStep the sampling step used when scanning the space for candidate pivots;
   *                  {@code 1} scans every element, larger values sample sparsely for speed
   * @param delta     the neighbourhood radius used to decide whether two elements are close
   *                  enough to be considered neighbours
   * @return the indices of the selected pivot elements, in selection order; may contain fewer
   *         than {@code maxPivots} elements if the space is fully covered earlier, or be empty
   *         if no pair of elements is found within {@code delta}, making pivot selection
   *         impossible; {@link NetAlgo} falls back to comparing every pair directly in that case
   */
  static int[] selectPivots(MetricSpace space, int maxPivots, int pivotStep, double delta) {
    Map<Integer, List<Integer>> neighbours = new HashMap<>();
    Map<Integer, Integer> elemSz = new HashMap<>();

    double delta2 = delta * delta;

    for (int elem1 = 0, len = space.getElementsNumber(); elem1 < len; elem1 += pivotStep) {
      for (int elem2 = elem1 + pivotStep; elem2 < len; elem2 += pivotStep) {
        double d = space.distanceEps2(elem1, elem2, delta2);

        if (d < delta)
          addElements(elem1, elem2, neighbours, elemSz);
      }
    }
    
    NavigableMap<Integer, Set<Integer>> neighbSz = makeNeighbSz(elemSz);
    int[] res = new int[maxPivots];
    Set<Integer> removed = new HashSet<>();
    int cnt = 0;

    while (!neighbSz.isEmpty()) {
      Map.Entry<Integer, Set<Integer>> nSzElems = neighbSz.lastEntry();
      Set<Integer> elems = nSzElems.getValue();
      int elemsSz = nSzElems.getKey();

      int b = NULL_ELEMENT;
      for (Iterator<Integer> i = elems.iterator(); i.hasNext(); ) {
        int b1 = i.next();

        i.remove();

        if (!removed.contains(b1)) {
          b = b1;
          break;
        }
      }

      if (elems.isEmpty()) {
        neighbSz.remove(elemsSz);

        if (b == NULL_ELEMENT)
          continue;
      }

      List<Integer> elemNs = neighbours.get(b);

      elemNs.removeAll(removed);

      int bNsz = elemNs.size();

      if (bNsz < elemsSz) {
        if (bNsz > 0)
          putSzB(bNsz, b, neighbSz);
      } else {
        res[cnt++] = b;

        if (cnt == maxPivots)
          break;

        removed.addAll(elemNs);
      }
    }

    if (cnt < maxPivots) {
      int[] r = new int[cnt];
      System.arraycopy(res, 0, r, 0, cnt);
      res = r;
    }

    return res;
  }


  /** Re-inserts element {@code b}, whose neighbour count changed to {@code sz}, into {@code neighbSz}. */
  private static void putSzB(int sz, int b, NavigableMap<Integer, Set<Integer>> neighbSz) {
    neighbSz.computeIfAbsent(sz, _x -> new HashSet<>()).add(b);
  }

  /** Records {@code elem1} and {@code elem2} as mutual neighbours. */
  private static void addElements(int elem1, int elem2, Map<Integer, List<Integer>> neighbours, Map<Integer, Integer> elemSz) {
    addElements0(elem1, elem2, neighbours, elemSz);
    addElements0(elem2, elem1, neighbours, elemSz);
  }

  /** Adds {@code elem1} to {@code elem}'s neighbour list and updates its recorded neighbour count. */
  private static void addElements0(int elem, int elem1, Map<Integer, List<Integer>> neighbours, Map<Integer, Integer> elemSz) {
    List<Integer> elemNeighbours = neighbours.computeIfAbsent(elem, _x -> new LinkedList<>());
    elemNeighbours.add(elem1);
    elemSz.put(elem, elemNeighbours.size());
  }

  /** Groups elements of {@code elemSz} by their neighbour count, ordered by count ascending. */
  private static NavigableMap<Integer, Set<Integer>> makeNeighbSz(Map<Integer, Integer> elemSz) {
    return elemSz.entrySet().stream().collect(
            Collectors.groupingBy(Map.Entry::getValue, TreeMap::new,
                    Collectors.mapping(Map.Entry::getKey, Collectors.toSet())));
  }
}