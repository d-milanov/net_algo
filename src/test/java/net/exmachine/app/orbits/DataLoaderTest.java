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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class DataLoaderTest {
  @TempDir
  Path tempDir;

  private File writeLines(String fileName, String... lines) throws IOException {
    Path file = tempDir.resolve(fileName);
    Files.write(file, List.of(lines), StandardCharsets.UTF_8);
    return file.toFile();
  }

  private File writeZip(String fileName, Map<String, String> entries) throws IOException {
    File file = tempDir.resolve(fileName).toFile();

    try (ZipOutputStream zip = new ZipOutputStream(new FileOutputStream(file))) {
      for (Map.Entry<String, String> entry : entries.entrySet()) {
        zip.putNextEntry(new ZipEntry(entry.getKey()));
        zip.write(entry.getValue().getBytes(StandardCharsets.UTF_8));
        zip.closeEntry();
      }
    }

    return file;
  }

  // --- Neodys ---

  @Test
  void neodysParsesSemiMajorAxisBasedRowIntoOrbitSpaceLayout() throws IOException {
    double a = 1.4581, e = 0.2226, iDeg = 10.828, OMDeg = 304.401, omDeg = 178.9;
    File f = writeLines("neodys.txt",
        "'433 Eros' 56000.0 " + a + " " + e + " " + iDeg + " " + OMDeg + " " + omDeg);

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(1, rows.length);
    assertOrbitRow(rows[0], a * (1 - e), e, iDeg, OMDeg, omDeg);
    assertEquals(a, rows[0][OrbitSpace.a], 1e-9);
  }

  @Test
  void neodysSkipsLinesThatDoNotMatchTheExpectedFormat() throws IOException {
    File f = writeLines("neodys.txt",
        "this line has no quoted name and is not orbital data",
        "'433 Eros' 56000.0 1.4581 0.2226 10.828 304.401 178.9");

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(1, rows.length);
  }

  @Test
  void neodysSkipsLinesWithMalformedNumericFields() throws IOException {
    // "1.2.3" matches the [\d\.E\+-]+ character class but is not a parsable double, so this
    // row reaches (and is dropped by) the NumberFormatException handler rather than the regex.
    File f = writeLines("neodys.txt",
        "'Bad Row' 56000.0 1.2.3 0.2226 10.828 304.401 178.9",
        "'433 Eros' 56000.0 1.4581 0.2226 10.828 304.401 178.9");

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(1, rows.length);
  }

  // --- zip archives ---

  @Test
  void loadsFromZipFileWhenArchiveContainsExactlyOneEntry() throws IOException {
    double a = 1.4581, e = 0.2226, iDeg = 10.828, OMDeg = 304.401, omDeg = 178.9;
    Map<String, String> entries = new LinkedHashMap<>();
    entries.put("neodys.txt", "'433 Eros' 56000.0 " + a + " " + e + " " + iDeg + " " + OMDeg + " " + omDeg);
    File f = writeZip("neodys.zip", entries);

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(1, rows.length);
    assertOrbitRow(rows[0], a * (1 - e), e, iDeg, OMDeg, omDeg);
  }

  @Test
  void zipExtensionIsRecognizedRegardlessOfCase() throws IOException {
    Map<String, String> entries = new LinkedHashMap<>();
    entries.put("neodys.txt", "'433 Eros' 56000.0 1.4581 0.2226 10.828 304.401 178.9");
    File f = writeZip("neodys.ZIP", entries);

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(1, rows.length);
  }

  @Test
  void fallsBackToReadingArchiveBytesAsTextWhenZipHasMultipleEntries() throws IOException {
    Map<String, String> entries = new LinkedHashMap<>();
    entries.put("first.txt", "'433 Eros' 56000.0 1.4581 0.2226 10.828 304.401 178.9");
    entries.put("second.txt", "'Extra' 56000.0 1.0 0.1 1.0 1.0 1.0");
    File f = writeZip("neodys.zip", entries);

    // None of the compressed archive's raw bytes happen to form a line matching the Neodys
    // format, so falling back to reading it as text yields no rows rather than an exception.
    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(0, rows.length);
  }

  @Test
  void fallsBackToReadingArchiveBytesAsTextWhenZipIsEmpty() throws IOException {
    File f = writeZip("empty.zip", Map.of());

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(0, rows.length);
  }

  // --- Mpcorb ---

  @Test
  void mpcorbParsesSemiMajorAxisBasedRowIntoOrbitSpaceLayout() throws IOException {
    double peri = 178.9, node = 304.401, incl = 10.828, e = 0.2226, a = 1.4581;
    File f = writeLines("mpcorb.txt",
        "K14X00A 19.7 0.15 K1670 123.4567 " + peri + " " + node + " " + incl + " " + e + " 0.2134 " + a);

    double[][] rows = new DataLoader.Mpcorb().load(f, null);

    assertEquals(1, rows.length);
    assertOrbitRow(rows[0], a * (1 - e), e, incl, node, peri);
    assertEquals(a, rows[0][OrbitSpace.a], 1e-9);
  }

  // --- Cams ---

  private static String camsLine(String name, double q, double e, double iDeg, double omDeg, double OMDeg) {
    String[] fields = new String[71];
    Arrays.fill(fields, "x");
    fields[0] = name;
    fields[59] = String.valueOf(q);
    fields[64] = String.valueOf(e);
    fields[66] = String.valueOf(iDeg);
    fields[68] = String.valueOf(omDeg);
    fields[70] = String.valueOf(OMDeg);
    return String.join(",", fields) + ",";
  }

  @Test
  void camsParsesPerihelionDistanceBasedRowIntoOrbitSpaceLayout() throws IOException {
    double q = 0.85, e = 0.6, iDeg = 30.0, omDeg = 120.0, OMDeg = 200.0;
    File f = writeLines("cams.csv", camsLine("20230101_00001", q, e, iDeg, omDeg, OMDeg));

    double[][] rows = new DataLoader.Cams().load(f, null);

    assertEquals(1, rows.length);
    assertOrbitRow(rows[0], q, e, iDeg, OMDeg, omDeg);
    assertEquals(q / (1 - e), rows[0][OrbitSpace.a], 1e-9);
  }

  @Test
  void camsWithParabolicEccentricityYieldsInfiniteSemiMajorAxis() throws IOException {
    File f = writeLines("cams.csv", camsLine("parabolic", 1.2, 1.0, 10.0, 0.0, 0.0));

    double[][] rows = new DataLoader.Cams().load(f, null);

    assertEquals(1, rows.length);
    assertEquals(Double.POSITIVE_INFINITY, rows[0][OrbitSpace.a]);
  }

  // --- filtering ---

  @Test
  void nullFilterKeepsEveryParsableRow() throws IOException {
    File f = writeLines("neodys.txt",
        "'A' 1.0 1.0 0.1 5.0 10.0 20.0",
        "'B' 1.0 2.0 0.5 15.0 30.0 40.0");

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(2, rows.length);
  }

  @Test
  void eccentricityFilterKeepsOnlyMatchingRows() throws IOException {
    // EFilter only looks at 'e', which is passed through identically by both accept() and
    // acceptQ(), so this assertion holds regardless of which of the two loadData() routes to.
    File f = writeLines("neodys.txt",
        "'Low' 1.0 1.0 0.1 5.0 10.0 20.0",
        "'High' 1.0 2.0 0.7 15.0 30.0 40.0");

    double[][] rows = new DataLoader.Neodys().load(f, new DataLoader.EFilter(0.0, 0.5));

    assertEquals(1, rows.length);
    assertEquals(0.1, rows[0][OrbitSpace.e], 1e-9);
  }

  @Test
  void filterRejectingEverythingReturnsNoRows() throws IOException {
    File f = writeLines("neodys.txt", "'A' 1.0 1.0 0.1 5.0 10.0 20.0");

    DataLoader.Filter rejectAll = (name, a, e, i, OM, om) -> false;
    double[][] rows = new DataLoader.Neodys().load(f, rejectAll);

    assertEquals(0, rows.length);
  }

  @Test
  void emptyFileYieldsNoRows() throws IOException {
    File f = writeLines("empty.txt");

    double[][] rows = new DataLoader.Neodys().load(f, null);

    assertEquals(0, rows.length);
  }

  // --- Filter.acceptQ ---

  @Test
  void acceptQDerivesSemiMajorAxisFromPerihelionDistance() {
    double[] capturedA = new double[1];
    DataLoader.Filter filter = new DataLoader.Filter() {
      @Override
      public boolean accept(String name, double a, double e, double i, double OM, double om) {
        capturedA[0] = a;
        return true;
      }
    };

    filter.acceptQ("x", 2.0, 0.5, 0, 0, 0);

    assertEquals(4.0, capturedA[0], 1e-9);
  }

  @Test
  void acceptQTreatsParabolicEccentricityAsInfiniteSemiMajorAxis() {
    double[] capturedA = new double[1];
    DataLoader.Filter filter = new DataLoader.Filter() {
      @Override
      public boolean accept(String name, double a, double e, double i, double OM, double om) {
        capturedA[0] = a;
        return true;
      }
    };

    filter.acceptQ("x", 2.0, 1.0, 0, 0, 0);

    assertEquals(Double.POSITIVE_INFINITY, capturedA[0]);
  }

  // --- AFilter / MainBelt / EFilter / Reduce ---

  @Test
  void aFilterAcceptsHalfOpenSemiMajorAxisRange() {
    DataLoader.AFilter filter = new DataLoader.AFilter(1.0, 2.0);

    assertTrue(filter.accept("x", 1.0, 0, 0, 0, 0));
    assertTrue(filter.accept("x", 1.999, 0, 0, 0, 0));
    assertFalse(filter.accept("x", 2.0, 0, 0, 0, 0));
    assertFalse(filter.accept("x", 0.999, 0, 0, 0, 0));
  }

  @Test
  void aFilterOnlyConsumesReduceCounterForInBoundsRows() {
    DataLoader.AFilter filter = new DataLoader.AFilter(0.0, 10.0, 2);

    assertTrue(filter.accept("x", 1.0, 0, 0, 0, 0));    // 1st in-bounds row -> kept
    assertFalse(filter.accept("x", 20.0, 0, 0, 0, 0));  // out of bounds, not counted
    assertFalse(filter.accept("x", 1.0, 0, 0, 0, 0));   // 2nd in-bounds row -> dropped
    assertTrue(filter.accept("x", 1.0, 0, 0, 0, 0));    // 3rd in-bounds row -> kept
  }

  @Test
  void mainBeltUsesStandardSemiMajorAxisBounds() {
    DataLoader.MainBelt filter = new DataLoader.MainBelt();

    assertTrue(filter.accept("x", 2.5, 0, 0, 0, 0));
    assertFalse(filter.accept("x", 1.0, 0, 0, 0, 0));
    assertFalse(filter.accept("x", 6.0, 0, 0, 0, 0));
  }

  @Test
  void eFilterAcceptsHalfOpenEccentricityRange() {
    DataLoader.EFilter filter = new DataLoader.EFilter(0.2, 0.8);

    assertTrue(filter.accept("x", 0, 0.2, 0, 0, 0));
    assertTrue(filter.accept("x", 0, 0.7999, 0, 0, 0));
    assertFalse(filter.accept("x", 0, 0.8, 0, 0, 0));
    assertFalse(filter.accept("x", 0, 0.1999, 0, 0, 0));
  }

  @Test
  void reduceAcceptsEveryNthCall() {
    DataLoader.Reduce reduce = new DataLoader.Reduce(3);

    assertTrue(reduce.accept("x", 0, 0, 0, 0, 0));
    assertFalse(reduce.accept("x", 0, 0, 0, 0, 0));
    assertFalse(reduce.accept("x", 0, 0, 0, 0, 0));
    assertTrue(reduce.accept("x", 0, 0, 0, 0, 0));
  }

  private static void assertOrbitRow(double[] row, double q, double e, double iDeg, double OMDeg, double omDeg) {
    double i = Math.toRadians(iDeg),
           OM = Math.toRadians(OMDeg),
           om = Math.toRadians(omDeg),
           p = q * (1.0 + e);

    assertEquals(q, row[OrbitSpace.q], 1e-9);
    assertEquals(e, row[OrbitSpace.e], 1e-9);
    assertEquals(i, row[OrbitSpace.i], 1e-9);
    assertEquals(OM, row[OrbitSpace.OM], 1e-9);
    assertEquals(om, row[OrbitSpace.om], 1e-9);
    assertEquals(p, row[OrbitSpace.p], 1e-9);
    assertEquals(Math.sqrt(p), row[OrbitSpace.sp], 1e-9);
    assertEquals(Math.sin(i), row[OrbitSpace.sinI], 1e-9);
    assertEquals(Math.cos(i), row[OrbitSpace.cosI], 1e-9);
    assertEquals(Math.sin(om), row[OrbitSpace.sinom], 1e-9);
    assertEquals(Math.cos(om), row[OrbitSpace.cosom], 1e-9);
    assertEquals(Math.sin(OM), row[OrbitSpace.sinOM], 1e-9);
    assertEquals(Math.cos(OM), row[OrbitSpace.cosOM], 1e-9);
  }
}
