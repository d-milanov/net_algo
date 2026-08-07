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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class DistanceMatrixTest {

  @Test
  void sizeIsZeroForEmptyMatrix() {
    DistanceMatrix matrix = new DistanceMatrix();

    assertEquals(0, matrix.size());
  }

  @Test
  void sizeCountsEachAddedPairOnce() {
    DistanceMatrix matrix = new DistanceMatrix();

    matrix.add(1, 2, 3.5);
    matrix.add(1, 3, 4.5);
    matrix.add(2, 3, 5.5);

    assertEquals(3, matrix.size());
  }

  @Test
  void addIsSymmetricRegardlessOfArgumentOrder() {
    DistanceMatrix matrix = new DistanceMatrix();

    matrix.add(2, 1, 3.5);

    assertEquals(1, matrix.size());
  }

  @Test
  void addingSamePairTwiceOverwritesAndDoesNotDuplicate() {
    DistanceMatrix matrix = new DistanceMatrix();

    matrix.add(1, 2, 3.5);
    matrix.add(2, 1, 9.9);

    assertEquals(1, matrix.size());
  }

  @Test
  void addAllowsSelfDistance() {
    DistanceMatrix matrix = new DistanceMatrix();

    matrix.add(1, 1, 0.0);

    assertEquals(1, matrix.size());
  }

  // --- deltaStat ---

  @Test
  void deltaStatMatchesForTwoEmptyMatrices() {
    DistanceMatrix.DeltaStat delta = DistanceMatrix.deltaStat(new DistanceMatrix(), new DistanceMatrix());

    assertTrue(delta.matches());
    assertEquals(0, delta.size1());
    assertEquals(0, delta.size2());
    assertEquals(0, delta.onlyInFirst());
    assertEquals(0, delta.onlyInSecond());
    assertEquals(0.0, delta.maxAbsoluteDelta());
  }

  @Test
  void deltaStatMatchesForIdenticalMatricesRegardlessOfPairArgumentOrder() {
    DistanceMatrix m1 = new DistanceMatrix();
    m1.add(1, 2, 1.0);
    m1.add(2, 3, 2.0);

    DistanceMatrix m2 = new DistanceMatrix();
    m2.add(2, 1, 1.0);
    m2.add(3, 2, 2.0);

    DistanceMatrix.DeltaStat delta = DistanceMatrix.deltaStat(m1, m2);

    assertTrue(delta.matches());
    assertEquals(2, delta.size1());
    assertEquals(2, delta.size2());
    assertEquals(0, delta.onlyInFirst());
    assertEquals(0, delta.onlyInSecond());
    assertEquals(0.0, delta.maxAbsoluteDelta());
  }

  @Test
  void deltaStatCountsPairsOnlyInFirstMatrix() {
    DistanceMatrix m1 = new DistanceMatrix();
    m1.add(1, 2, 1.0);
    m1.add(1, 3, 2.0);

    DistanceMatrix m2 = new DistanceMatrix();
    m2.add(1, 2, 1.0);

    DistanceMatrix.DeltaStat delta = DistanceMatrix.deltaStat(m1, m2);

    assertFalse(delta.matches());
    assertEquals(2, delta.size1());
    assertEquals(1, delta.size2());
    assertEquals(1, delta.onlyInFirst());
    assertEquals(0, delta.onlyInSecond());
  }

  @Test
  void deltaStatCountsPairsOnlyInSecondMatrix() {
    DistanceMatrix m1 = new DistanceMatrix();
    m1.add(1, 2, 1.0);

    DistanceMatrix m2 = new DistanceMatrix();
    m2.add(1, 2, 1.0);
    m2.add(1, 3, 2.0);

    DistanceMatrix.DeltaStat delta = DistanceMatrix.deltaStat(m1, m2);

    assertFalse(delta.matches());
    assertEquals(1, delta.size1());
    assertEquals(2, delta.size2());
    assertEquals(0, delta.onlyInFirst());
    assertEquals(1, delta.onlyInSecond());
  }

  @Test
  void deltaStatReportsMaxAbsoluteDeltaForDivergingCommonPairs() {
    DistanceMatrix m1 = new DistanceMatrix();
    m1.add(1, 2, 1.0);
    m1.add(1, 3, 5.0);

    DistanceMatrix m2 = new DistanceMatrix();
    m2.add(1, 2, 1.25);
    m2.add(1, 3, 4.9);

    DistanceMatrix.DeltaStat delta = DistanceMatrix.deltaStat(m1, m2);

    assertFalse(delta.matches());
    assertEquals(0, delta.onlyInFirst());
    assertEquals(0, delta.onlyInSecond());
    assertEquals(0.25, delta.maxAbsoluteDelta(), 1e-9);
  }

  // --- Stub ---

  @Test
  void stubCountsAddsWithoutDeduplicatingPairs() {
    DistanceMatrix.Stub stub = new DistanceMatrix.Stub();

    stub.add(1, 2, 3.5);
    stub.add(1, 2, 9.9); // same pair again: a real DistanceMatrix would dedupe, Stub just counts
    stub.add(2, 3, 1.0);

    assertEquals(3, stub.size());
  }

  @Test
  void stubToStringReportsCount() {
    DistanceMatrix.Stub stub = new DistanceMatrix.Stub();

    stub.add(1, 2, 3.5);
    stub.add(2, 3, 1.0);

    assertEquals("Size: 2", stub.toString());
  }
}