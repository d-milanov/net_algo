// This algorithm calculates all distances which do not exceed a fixed value epsilon,
// for a given set of points and a distance function r(x, y).
// The algorithm works faster than the complete search in the case when the distance
// function r(x, y) complies with the relaxed triangle inequality:
//                          r(x, y) <= K(r(x, z) + r(y, z)),
// for any three data points x, y, z, and a constant K.
//
// Author: Danila Milanov (danila.milanov@gmail.com)
// Year:   2026

package net.exmachine.app.orbits;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 * Loads orbital elements from tab-separated and CSV catalog formats into {@link OrbitSpace}
 * rows, i.e. {@code double[]} arrays indexed by the column constants declared on
 * {@link OrbitSpace} ({@code a}, {@code q}, {@code p}, {@code e}, {@code i}, {@code OM},
 * {@code om} and their cached sines/cosines).
 */
public abstract class DataLoader {
  public interface Filter {
    boolean accept(String name, double a, double e, double i, double OM, double om);
    default boolean acceptQ(String name, double q, double e, double i, double OM, double om){
      return accept(name, Math.abs(e - 1.0) < Util.PRECISION ? Double.POSITIVE_INFINITY : q / (1-e), e, i, OM, om);
    }
  }

  /**
   * Loads orbit rows from the given file, keeping only the ones accepted by {@code filter}.
   *
   * @param f      the catalog file to read
   * @param filter accepts or rejects a row by its orbital elements, or {@code null} to accept
   *               every row
   * @return the accepted orbits, one row per element, in {@link OrbitSpace} column layout
   */
  public abstract double[][] load(File f, Filter filter) throws IOException;

  public static class Neodys extends DataLoader {
    @Override
    public double[][] load(File f, Filter filter) throws IOException {
      return loadDataNeodys(f, filter);
    }
  }

  public static class Mpcorb extends DataLoader {
    @Override
    public double[][] load(File f, Filter filter) throws IOException {
      return loadDataMpcorb(f, filter);
    }
  }

  public static class Cams extends DataLoader {
    @Override
    public double[][] load(File f, Filter filter) throws IOException {
      return loadDataCams(f, filter);
    }
  }

  static double[][] loadDataNeodys(File f, Filter filter) throws IOException {
    return loadData(f, Pattern.compile(
        "^'([^']+)'\\s+\\d+\\.\\d+\\s+" +
        "([\\d.E+-]+)\\s+" +
        "([\\d.E+-]+)\\s+" +
        "([\\d.E+-]+)\\s+" +
        "([\\d.E+-]+)\\s+" +
        "([\\d.E+-]+)"),
        new int[] {1, 2, 3, 4, 5, 6},
        false,
        filter);
  }

  static double[][] loadDataMpcorb(File f, Filter filter) throws IOException {
    return loadData(f, Pattern.compile(
        "^(\\S+)\\s+" +                  // Des'n
        "\\d+\\.\\d+\\s+" +              // H
        "\\d+\\.\\d+\\s+" +              // G
        "\\S+\\s+" +                     // Epoch
        "\\d+\\.\\d+\\s+" +              // M
        "([\\d\\.E\\+-]+)\\s+" +         // Peri
        "([\\d\\.E\\+-]+)\\s+" +         // Node
        "([\\d\\.E\\+-]+)\\s+" +         // Incl
        "([\\d\\.E\\+-]+)\\s+" +         // e
        "\\d+\\.\\d+\\s+" +              // n
        "([\\d\\.E\\+-]+)"),             // a
        new int[] {1, 6, 5, 4, 3, 2},
        false,
        filter);
  }

  static double[][] loadDataCams(File f, Filter filter) throws IOException {
    return loadData(f, Pattern.compile(
        "^([^,]+)," +      //#
        "(?:[^,]+,){58}" +
        "([^,]+)," +       //q
        "(?:[^,]+,){4}" +
        "([^,]+)," +       //e
        "(?:[^,]+,)" +
        "([^,]+)," +       //i
        "(?:[^,]+,)" +
        "([^,]+)," +       //om
        "(?:[^,]+,)" +
        "([^,]+),"),       //Om
        new int[] {1, 2, 3, 4, 6, 5},
        true,
        filter);
  }

  static double[][] loadData(File f, Pattern p, int[] name_a_e_i_Om_om, boolean isQ, Filter filter)
      throws IOException {
    List<double[]> res = new ArrayList<>();
    PrintStream dbgOut = System.out;

    try (BufferedReader r = new BufferedReader(new InputStreamReader(openInput(f)))) {
      String l;
      int bCnt = 0,
          vCnt = 0;

      while ((l = r.readLine()) != null) {
        Matcher m = p.matcher(l);

        if (!m.find())
          continue;

        String name = m.group(name_a_e_i_Om_om[0]),
               aqs  = m.group(name_a_e_i_Om_om[1]),
               es   = m.group(name_a_e_i_Om_om[2]),
               is   = m.group(name_a_e_i_Om_om[3]),
               OMs  = m.group(name_a_e_i_Om_om[4]),
               oms  = m.group(name_a_e_i_Om_om[5]);

        double aq, e, i, OM, om;

        try {
          aq = Double.parseDouble(aqs);
          e = Double.parseDouble(es);
          i = Double.parseDouble(is);
          OM = Double.parseDouble(OMs);
          om = Double.parseDouble(oms);
        } catch (NumberFormatException ex) {
          dbgOut.println("Cannot parse data: " + l);
          continue;
        }

        if (filter == null || 
                isQ && filter.acceptQ(name, aq, e, i, OM, om) || !isQ && filter.accept(name, aq, e, i, OM, om)) {
          res.add(isQ ? rowFromQE(aq, e, i, OM, om) : rowFromAE(aq, e, i, OM, om));
          bCnt++;
        }

        vCnt++;
      }

      dbgOut.printf("Loaded %d orbits (%.3f%% of all valid rows)%n",
          bCnt, vCnt != 0 ? 100.0 * bCnt / (double) vCnt : 0.0);
    }

    return res.toArray(new double[0][]);
  }

  /**
   * Opens {@code f} for reading. If its name ends in {@code .zip} and it contains exactly one
   * entry, the entry's decompressed content is read instead of the raw archive bytes; closing
   * the returned stream also closes the underlying {@link ZipFile}. Any other file, including a
   * {@code .zip} with zero or more than one entry, is read as plain text.
   */
  private static InputStream openInput(File f) throws IOException {
    if (!f.getName().toLowerCase(Locale.ROOT).endsWith(".zip"))
      return new FileInputStream(f);

    ZipFile zip = new ZipFile(f);
    Enumeration<? extends ZipEntry> entries = zip.entries();

    if (entries.hasMoreElements()) {
      ZipEntry entry = entries.nextElement();

      if (!entries.hasMoreElements())
        return zipEntryStream(zip, entry);
    }

    zip.close();

    return new FileInputStream(f);
  }

  /** Wraps a zip entry's stream so closing it also closes the owning {@code zip}. */
  private static InputStream zipEntryStream(ZipFile zip, ZipEntry entry) throws IOException {
    return new FilterInputStream(zip.getInputStream(entry)) {
      @Override
      public void close() throws IOException {
        try {
          super.close();
        } finally {
          zip.close();
        }
      }
    };
  }

  private static double[] rowFromAE(double a, double e, double iDeg, double OMDeg, double omDeg) {
    return row(a * (1.0 - e), e, iDeg, OMDeg, omDeg);
  }

  private static double[] rowFromQE(double q, double e, double iDeg, double OMDeg, double omDeg) {
    return row(q, e, iDeg, OMDeg, omDeg);
  }

  private static double[] row(double q, double e, double iDeg, double OMDeg, double omDeg) {
    double i = Math.toRadians(iDeg),
           OM = Math.toRadians(OMDeg),
           om = Math.toRadians(omDeg),
           p = q * (1.0 + e);

    double[] row = new double[14];

    row[OrbitSpace.a] = e == 1.0 ? Double.POSITIVE_INFINITY : q / (1.0 - e);
    row[OrbitSpace.q] = q;
    row[OrbitSpace.p] = p;
    row[OrbitSpace.e] = e;
    row[OrbitSpace.i] = i;
    row[OrbitSpace.OM] = OM;
    row[OrbitSpace.om] = om;
    row[OrbitSpace.sinI] = Math.sin(i);
    row[OrbitSpace.cosI] = Math.cos(i);
    row[OrbitSpace.sinom] = Math.sin(om);
    row[OrbitSpace.cosom] = Math.cos(om);
    row[OrbitSpace.sinOM] = Math.sin(OM);
    row[OrbitSpace.cosOM] = Math.cos(OM);
    row[OrbitSpace.sp] = Math.sqrt(p);

    return row;
  }

  public static class AFilter implements Filter {
    final double a_min;
    final double a_max;
    final Reduce reduce;

    public AFilter(double a_min, double a_max) {
      this(a_min, a_max, 1L);
    }

    public AFilter(double a_min, double a_max, long times) {
      this.a_min = a_min;
      this.a_max = a_max;
      this.reduce = new Reduce(times);
    }

    @Override
    public boolean accept(String name, double a, double e, double i, double OM, double om) {
      return a >= a_min && a < a_max &&
             reduce.accept(name, a, e, i, OM, om);
    }

    @Override
    public String toString() {
      return String.format("%f <= a < %f, times = %d", a_min, a_max, reduce.times);
    }
  }

  public static class MainBelt extends AFilter {
    public MainBelt() {
      this(1L);
    }

    public MainBelt(long times) {
      super(1.523679, 5.2044, times);
    }
  }

  public static class EFilter implements Filter {
    final double e_min;
    final double e_max;

    public EFilter(double e_min, double e_max) {
      this.e_min = e_min;
      this.e_max = e_max;
    }

    @Override
    public boolean accept(String name, double a, double e, double i, double OM, double om) {
      return e >= e_min && e < e_max;
    }

    @Override
    public String toString() {
      return String.format("%f <= a < %f", e_min, e_max);
    }
  }

  public static class Reduce implements Filter {
    final long times;
    long cnt;

    public Reduce(long times) {
      this.times = times;
    }

    @Override
    public boolean accept(String name, double a, double e, double i, double OM, double om) {
      return ((cnt++) % times) == 0;
    }

    @Override
    public String toString() {
      return String.format("times = %d", times);
    }
  }
}
