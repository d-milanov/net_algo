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

import net.exmachine.algo.distance.MetricSpace;

/**
 * A {@link MetricSpace} decorator that counts how many times {@link #distanceEps2} is called,
 * so the number of distance evaluations an algorithm performs can be measured directly instead
 * of only its wall-clock time.
 */
final class CountingMetricSpace implements MetricSpace {
  private final MetricSpace delegate;
  private long callCount;

  CountingMetricSpace(MetricSpace delegate) {
    this.delegate = delegate;
  }

  long callCount() {
    return callCount;
  }

  @Override
  public int getElementsNumber() {
    return delegate.getElementsNumber();
  }

  @Override
  public double getTriangleCoefficient() {
    return delegate.getTriangleCoefficient();
  }

  @Override
  public double distanceEps2(int elemA, int elemB, double epsSquared) {
    callCount++;
    return delegate.distanceEps2(elemA, elemB, epsSquared);
  }
}
