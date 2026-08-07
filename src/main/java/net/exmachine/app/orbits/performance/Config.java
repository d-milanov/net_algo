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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Map;

import org.yaml.snakeyaml.Yaml;

import net.exmachine.algo.distance.MetricSpace;
import net.exmachine.app.orbits.D1;
import net.exmachine.app.orbits.D2;
import net.exmachine.app.orbits.D_D;
import net.exmachine.app.orbits.D_H;
import net.exmachine.app.orbits.DataLoader;
import net.exmachine.app.orbits.Rho;

/**
 * Parameters for a {@link Test} run, loaded from a human-readable YAML file. See
 * {@code conf/sample_conf.yaml} for the expected structure and a description of each field.
 */
public final class Config {
  public final PrintStream out;

  public final boolean checkEnabled;
  public final DataSet checkData;

  public final boolean testEnabled;
  public final DataSet testData;

  public final double eps;
  public final int nTests;

  public final String distanceName;
  public final double eMax;

  public final double delta;
  public final int maxPivots;
  public final int pivotStep;

  private Config(PrintStream out, boolean checkEnabled, DataSet checkData, boolean testEnabled, DataSet testData,
                 double eps, int nTests, String distanceName, double eMax,
                 double delta, int maxPivots, int pivotStep) {
    this.out = out;
    this.checkEnabled = checkEnabled;
    this.checkData = checkData;
    this.testEnabled = testEnabled;
    this.testData = testData;
    this.eps = eps;
    this.nTests = nTests;
    this.distanceName = distanceName;
    this.eMax = eMax;
    this.delta = delta;
    this.maxPivots = maxPivots;
    this.pivotStep = pivotStep;
  }

  /**
   * Loads and validates a configuration from the given YAML file.
   *
   * @throws IllegalArgumentException if a required key or section is missing, or a value has
   *                                   the wrong type
   */
  public static Config load(File file) throws IOException {
    Map<String, Object> root;

    try (InputStream in = new FileInputStream(file)) {
      root = new Yaml().load(in);
    }

    Map<String, Object> check = section(root, "check");
    Map<String, Object> test = section(root, "test");
    Map<String, Object> distance = requireSection(root, "distance");
    Map<String, Object> net = requireSection(root, "net");

    return new Config(
        openOut(root.get("out")),
        check != null && bool(check, "enabled", false),
        check == null ? null : DataSet.of(check),
        test != null && bool(test, "enabled", false),
        test == null ? null : DataSet.of(test),
        number(root, "eps").doubleValue(),
        number(root, "nTests").intValue(),
        string(distance, "name"),
        distance.containsKey("eMax") ? number(distance, "eMax").doubleValue() : 1.0,
        number(net, "delta").doubleValue(),
        number(net, "maxPivots").intValue(),
        number(net, "pivotStep").intValue());
  }

  /** Builds the {@link MetricSpace} named by {@code distance.name} over the given orbit rows. */
  public MetricSpace makeSpace(double[][] elements) {
    switch (distanceName) {
      case "D1": return new D1(elements, eMax);
      case "D2": return new D2(elements);
      case "DD": return new D_D(elements, eMax);
      case "DH": return new D_H(elements, eMax);
      case "Rho": return new Rho(elements);
      default: throw new IllegalArgumentException("Unknown distance.name: " + distanceName);
    }
  }

  /**
   * A configured data source: which {@link DataLoader} to use, which file to read, and an
   * optional {@link DataLoader.Filter} to restrict which rows are kept.
   */
  public static final class DataSet {
    public final DataLoader loader;
    public final File file;
    public final DataLoader.Filter filter;

    private DataSet(DataLoader loader, File file, DataLoader.Filter filter) {
      this.loader = loader;
      this.file = file;
      this.filter = filter;
    }

    static DataSet of(Map<String, Object> m) {
      Map<String, Object> filter = section(m, "filter");
      return new DataSet(loader(string(m, "dataLoader")), new File(string(m, "dataFile")),
          filter == null ? null : filter(filter));
    }

    private static DataLoader loader(String name) {
      switch (name) {
        case "Neodys": return new DataLoader.Neodys();
        case "Mpcorb": return new DataLoader.Mpcorb();
        case "Cams": return new DataLoader.Cams();
        default: throw new IllegalArgumentException("Unknown dataLoader: " + name);
      }
    }

    private static DataLoader.Filter filter(Map<String, Object> m) {
      switch (string(m, "type")) {
        case "EFilter": return new DataLoader.EFilter(number(m, "eMin").doubleValue(), number(m, "eMax").doubleValue());
        case "MainBelt": return m.containsKey("times")
            ? new DataLoader.MainBelt(number(m, "times").longValue())
            : new DataLoader.MainBelt();
        case "Reduce": return new DataLoader.Reduce(number(m, "times").longValue());
        default: throw new IllegalArgumentException("Unknown filter.type: " + m.get("type"));
      }
    }
  }

  private static PrintStream openOut(Object v) throws IOException {
    String s = v == null ? "System.out" : v.toString();
    return "System.out".equals(s) ? System.out : new PrintStream(new FileOutputStream(s));
  }

  @SuppressWarnings("unchecked")
  private static Map<String, Object> section(Map<String, Object> root, String key) {
    return (Map<String, Object>) root.get(key);
  }

  private static Map<String, Object> requireSection(Map<String, Object> root, String key) {
    Map<String, Object> s = section(root, key);

    if (s == null)
      throw new IllegalArgumentException("Missing required config section: " + key);

    return s;
  }

  private static Number number(Map<String, Object> m, String key) {
    Object v = m.get(key);

    if (!(v instanceof Number))
      throw new IllegalArgumentException("Missing or non-numeric config key: " + key);

    return (Number) v;
  }

  private static String string(Map<String, Object> m, String key) {
    Object v = m.get(key);

    if (v == null)
      throw new IllegalArgumentException("Missing required config key: " + key);

    return v.toString();
  }

  private static boolean bool(Map<String, Object> m, String key, boolean def) {
    Object v = m.get(key);
    return v == null ? def : (Boolean) v;
  }
}
