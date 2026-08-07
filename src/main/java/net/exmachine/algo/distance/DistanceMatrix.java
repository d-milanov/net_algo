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


import java.util.HashMap;
import java.util.Map;

/**
 * Distance matrix: the algorithm execution result storage class.
 * <p>
 * Stores pairwise distances between elements identified by integer index. Since distance
 * is assumed to be symmetric ({@code r(x, y) == r(y, x)}), each pair is stored only once,
 * keyed by the smaller index first.
 */
public class DistanceMatrix {
  private final Map<Integer, Map<Integer, Double>> matrix = new HashMap<>();

  /**
   * Records the distance between two elements.
   * <p>
   * The pair is stored in a canonical order (smaller index first), so calling this method
   * with {@code (elemA, elemB)} or {@code (elemB, elemA)} has the same effect. If a distance
   * was already recorded for the pair, it is overwritten.
   *
   * @param elemA index of the first element
   * @param elemB index of the second element
   * @param dist  distance between {@code elemA} and {@code elemB}
   */
  public void add(int elemA, int elemB, double dist) {
    if (elemA > elemB) {
      int tmp = elemB;
      elemB = elemA;
      elemA = tmp;
    }

    matrix.computeIfAbsent(elemA, _x -> new HashMap<>()).put(elemB, dist);
  }


  /**
   * Returns the number of distances stored in the matrix.
   *
   * @return the count of recorded element pairs
   */
  public int size() {
    return matrix.values().stream().mapToInt(Map::size).sum();
  }

  /**
   * Compares two distance matrices pair by pair, e.g. as a correctness check between two
   * algorithms expected to compute the same matrix (such as {@link NetAlgo} against a
   * brute-force baseline).
   *
   * @param m1 the first matrix
   * @param m2 the second matrix
   * @return a summary of the pairs found in only one of the two matrices, and the largest
   *         distance discrepancy among pairs found in both; {@link DeltaStat#matches()} is
   *         {@code true} when the matrices agree exactly
   */
  public static DeltaStat deltaStat(DistanceMatrix m1, DistanceMatrix m2) {
    long onlyInFirst = 0;
    double maxAbsoluteDelta = 0.0;

    for (Map.Entry<Integer, Map<Integer, Double>> row : m1.matrix.entrySet()) {
      int elemA = row.getKey();
      Map<Integer, Double> distances1 = row.getValue(),
                            distances2 = m2.matrix.get(elemA);

      if (distances2 == null) {
        onlyInFirst += distances1.size();
        continue;
      }

      for (Map.Entry<Integer, Double> entry : distances1.entrySet()) {
        Double d2 = distances2.get(entry.getKey());

        if (d2 == null) {
          onlyInFirst++;
          continue;
        }

        maxAbsoluteDelta = Math.max(maxAbsoluteDelta, Math.abs(entry.getValue() - d2));
      }
    }

    long size1 = m1.size(),
         size2 = m2.size(),
         onlyInSecond = size2 - (size1 - onlyInFirst);

    return new DeltaStat(size1, size2, onlyInFirst, onlyInSecond, maxAbsoluteDelta);
  }

  /**
   * A {@link DistanceMatrix} that only counts {@link #add} calls instead of storing the
   * distances, for use in timing tests where the storage overhead of the full matrix would
   * skew the measurement.
   */
  public static class Stub extends DistanceMatrix {
    private long count;

    @Override
    public void add(int elemA, int elemB, double dist) {
      count++;
    }

    @Override
    public int size() {
      return (int) count;
    }

    @Override
    public String toString() {
      return "Size: " + count;
    }
  }

  /** Summary of the differences found by {@link DistanceMatrix#deltaStat}. */
  public static final class DeltaStat {
    private final long size1;
    private final long size2;
    private final long onlyInFirst;
    private final long onlyInSecond;
    private final double maxAbsoluteDelta;

    private DeltaStat(long size1, long size2, long onlyInFirst, long onlyInSecond, double maxAbsoluteDelta) {
      this.size1 = size1;
      this.size2 = size2;
      this.onlyInFirst = onlyInFirst;
      this.onlyInSecond = onlyInSecond;
      this.maxAbsoluteDelta = maxAbsoluteDelta;
    }

    /** Number of pairs recorded in the first matrix. */
    public long size1() {
      return size1;
    }

    /** Number of pairs recorded in the second matrix. */
    public long size2() {
      return size2;
    }

    /** Number of pairs recorded in the first matrix but not the second. */
    public long onlyInFirst() {
      return onlyInFirst;
    }

    /** Number of pairs recorded in the second matrix but not the first. */
    public long onlyInSecond() {
      return onlyInSecond;
    }

    /** Largest absolute distance discrepancy among pairs recorded in both matrices. */
    public double maxAbsoluteDelta() {
      return maxAbsoluteDelta;
    }

    /** Whether the two matrices record exactly the same pairs with exactly the same distances. */
    public boolean matches() {
      return onlyInFirst == 0 && onlyInSecond == 0 && maxAbsoluteDelta == 0.0;
    }

    @Override
    public String toString() {
      return String.format("#M1 = %d  #M2 = %d  #(M1 \\ M2) = %d  #(M2 \\ M1) = %d  max|d1 - d2| = %f",
          size1, size2, onlyInFirst, onlyInSecond, maxAbsoluteDelta);
    }
  }
}
