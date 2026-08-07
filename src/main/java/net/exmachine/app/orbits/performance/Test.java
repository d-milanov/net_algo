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

import java.io.File;
import java.io.IOException;
import java.util.Locale;

import net.exmachine.algo.distance.DistanceMatrix;
import net.exmachine.algo.distance.MetricSpace;
import net.exmachine.algo.distance.NetAlgo;

/**
 * Net algorithm performance test program.
 * <p>
 * Measures {@link NetAlgo} against {@link BruteForce}, a plain complete-search baseline: a
 * correctness check confirms both find the same set of close pairs, and a timing test compares
 * their running time. Test and algorithm parameters are read from a YAML configuration file; see
 * {@code conf/sample_conf.yaml}.
 * <p>
 * Usage: {@code Test <config-file.yaml>}.
 */
public class Test {
  public static void main(String[] args) throws IOException {
    Locale.setDefault(new Locale("us"));

    new Test().perform(new File(args[0]));
  }

  public void perform(File confFile) throws IOException {
    Config conf = Config.load(confFile);

    conf.out.printf("Testing Net algorithm (distance=%s, eps=%s, delta=%s, maxPivots=%d, pivotStep=%d) "
            + "against brute force%n",
        conf.distanceName, conf.eps, conf.delta, conf.maxPivots, conf.pivotStep);

    if (conf.checkEnabled)
      doCheck(conf);
    else
      conf.out.println("\nCorrectness check is skipped");

    if (conf.testEnabled)
      doTest(conf);
    else
      conf.out.println("\nTiming test is skipped");
  }

  private void doCheck(Config conf) throws IOException {
    conf.out.println("\nCorrectness check");

    MetricSpace space = conf.makeSpace(conf.checkData.loader.load(conf.checkData.file, conf.checkData.filter));

    DistanceMatrix netResult = NetAlgo.calculate(space, conf.eps, conf.delta, conf.maxPivots, conf.pivotStep, null);
    DistanceMatrix bfResult = BruteForce.calculate(space, conf.eps, null);

    DistanceMatrix.DeltaStat delta = DistanceMatrix.deltaStat(netResult, bfResult);

    conf.out.println(delta);
    conf.out.println(delta.matches() ? "OK" : "MISMATCH");
  }

  private void doTest(Config conf) throws IOException {
    conf.out.println("\nTiming test");

    MetricSpace space = conf.makeSpace(conf.testData.loader.load(conf.testData.file, conf.testData.filter));

    CountingMetricSpace netSpace = new CountingMetricSpace(space);
    CountingMetricSpace bfSpace = new CountingMetricSpace(space);

    Stat netStat = new Stat(),
         bfStat = new Stat();

    for (int i = 0; i < conf.nTests; i++) {
      netStat.record(() -> NetAlgo.calculate(
          netSpace, conf.eps, conf.delta, conf.maxPivots, conf.pivotStep, new DistanceMatrix.Stub(), null));
      bfStat.record(() -> BruteForce.calculate(bfSpace, conf.eps, new DistanceMatrix.Stub(), null));
    }

    conf.out.println("Net algo stat:");
    netStat.print(conf.out);
    conf.out.println("  distanceEps2 calls: " + netSpace.callCount());

    conf.out.println("Brute force stat:");
    bfStat.print(conf.out);
    conf.out.println("  distanceEps2 calls: " + bfSpace.callCount());

    conf.out.printf("Avg brute force / avg net: %.3f%n", bfStat.avgTimeMs() / netStat.avgTimeMs());
  }
}
