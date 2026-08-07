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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

/** Collects wall-clock timings across repeated runs of a timed action, for later reporting. */
final class Stat {
  private final List<Long> timesNanos = new ArrayList<>();

  /** Runs {@code action}, records its wall-clock duration, and returns its result. */
  <T> T record(Supplier<T> action) {
    long start = System.nanoTime();
    T result = action.get();

    timesNanos.add(System.nanoTime() - start);

    return result;
  }

  /** Average duration across all recorded runs, in milliseconds; {@code 0} if none were recorded. */
  double avgTimeMs() {
    return timesNanos.stream().mapToLong(Long::longValue).average().orElse(0.0) / 1e6;
  }

  void print(PrintStream out) {
    if (timesNanos.isEmpty()) {
      out.println("  no runs recorded");
      return;
    }

    long minNanos = timesNanos.stream().mapToLong(Long::longValue).min().getAsLong();
    long maxNanos = timesNanos.stream().mapToLong(Long::longValue).max().getAsLong();

    out.printf("  runs: %d, avg: %.3f ms, min: %.3f ms, max: %.3f ms%n",
        timesNanos.size(), avgTimeMs(), minNanos / 1e6, maxNanos / 1e6);
  }
}
