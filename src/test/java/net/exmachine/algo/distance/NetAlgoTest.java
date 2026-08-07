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

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class NetAlgoTest {

  @Test
  void calculateFindsAllPairsWithinEpsForClusteredPoints() {
    // Three well-separated clusters of four consecutive points each.
    double[] points = {0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23};
    double eps = 1.5;

    DistanceMatrix result = NetAlgo.calculate(
            new LinePointsMetricSpace(points), eps, 2.0, 4, 1, null);

    assertEquals(bruteForceCount(points, eps), result.size());
  }

  @Test
  void fallsBackToBruteForceWhenDeltaIsTooSmallForAnyPivot() {
    // No pair is within delta = 1.0 of each other, so no pivot can be formed; the algorithm
    // must fall back to comparing every pair directly instead of failing.
    double[] points = {0, 1};
    double eps = 1.5;

    DistanceMatrix result = NetAlgo.calculate(
            new LinePointsMetricSpace(points), eps, 1.0, 4, 1, null);

    assertEquals(bruteForceCount(points, eps), result.size());
  }

  @Test
  void calculateFindsAllPairsWhenEpsCoversEntireSpace() {
    double[] points = {0, 1, 2, 3, 4, 5};
    double eps = 100.0;

    DistanceMatrix result = NetAlgo.calculate(
            new LinePointsMetricSpace(points), eps, 1.1, 3, 1, null);

    int totalPairs = points.length * (points.length - 1) / 2;
    assertEquals(totalPairs, result.size());
    assertEquals(bruteForceCount(points, eps), result.size());
  }

  @Test
  void calculateReturnsEmptyMatrixWhenNoPairIsWithinEps() {
    // A pivot can still be formed from the close pair (0, 1), but eps is too small
    // for even that pair to qualify.
    double[] points = {0, 1, 100};
    double eps = 0.5;

    DistanceMatrix result = NetAlgo.calculate(
            new LinePointsMetricSpace(points), eps, 2.0, 2, 1, null);

    assertEquals(0, result.size());
    assertEquals(bruteForceCount(points, eps), result.size());
  }

  @Test
  void fallsBackToBruteForceWhenSparsePivotSamplingFindsNoPivots() {
    // pivotStep = 3 only samples points 0, 3, 6, 9 as pivot candidates, none of which are
    // within delta = 1.5 of each other, so no pivot can be formed from the sample even though
    // the full point set has plenty of close pairs.
    double[] points = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    double eps = 1.5;

    DistanceMatrix result = NetAlgo.calculate(
            new LinePointsMetricSpace(points), eps, 1.5, 3, 3, null);

    assertEquals(bruteForceCount(points, eps), result.size());
  }

  @Test
  void calculateNotifiesListenerOfProgress() {
    double[] points = {0, 1, 2, 3, 10, 11, 12, 13};
    RecordingListener listener = new RecordingListener();

    NetAlgo.calculate(new LinePointsMetricSpace(points), 1.5, 10, 3, 1, listener);

    assertEquals(List.of("calculation started", "fillAdjacent started", "fillAdjacent ended"),
            listener.events.subList(0, 3));
    assertEquals("calculation ended", listener.events.get(listener.events.size() - 1));
    assertEquals(points.length,
            listener.events.stream().filter("element processed"::equals).count());
  }

  @Test
  void calculateReturnsEmptyMatrixForEmptySpace() {
    DistanceMatrix result = NetAlgo.calculate(
            new LinePointsMetricSpace(), 1.5, 1.0, 3, 1, null);

    assertEquals(0, result.size());
  }


  private static int bruteForceCount(double[] points, double eps) {
    int count = 0;

    for (int i = 0; i < points.length; i++) {
      for (int j = i + 1; j < points.length; j++) {
        if (Math.abs(points[i] - points[j]) < eps) count++;
      }
    }

    return count;
  }

  private static class RecordingListener implements NetAlgoListener {
    final List<String> events = new ArrayList<>();

    @Override
    public void onEvent(String name) {
      events.add(name);
    }

    @Override
    public void info(String message) {
      // Not asserted on; info messages are diagnostic only.
    }
  }
}