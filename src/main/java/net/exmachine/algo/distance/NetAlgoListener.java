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
 * Callback interface notified of {@link NetAlgo} execution progress.
 * <p>
 * Implementations may be used for logging, progress reporting, or diagnostics; the algorithm
 * itself does not depend on the listener's behavior.
 */
public interface NetAlgoListener {
  /**
   * Called when the algorithm reaches a named stage of its execution
   * (e.g. {@code "calculation started"}, {@code "element processed"}).
   *
   * @param name identifier of the stage or event that occurred
   */
  void onEvent(String name);

  /**
   * Called to report a free-form diagnostic or informational message.
   *
   * @param message the message to report
   */
  void info(String message);
}
