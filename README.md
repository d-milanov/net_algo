# net_algo

A Java implementation of the algorithm presented in the article
["Relaxed triangle inequality for the orbital similarity criterion by Southworth and Hawkins and its variants"](https://link.springer.com/article/10.1007/s10569-019-9884-6).
Please cite it if you use this code in your public work.

This algorithm calculates all distances which do not exceed a fixed value `epsilon`, for a given
set of points and a distance function `r(x, y)`. The algorithm works faster than the complete
search in the case when the distance function `r(x, y)` complies with the relaxed triangle
inequality:

```
r(x, y) <= K(r(x, z) + r(y, z)),
```

for any three data points `x`, `y`, `z`, and a constant `K`.

## Build

```sh
./gradlew build
```

On Windows, use `gradlew.bat build` instead.

## Test

```sh
./gradlew test
```

## Usage

To use the library in your own code, implement the
[`MetricSpace`](src/main/java/net/exmachine/algo/distance/MetricSpace.java) interface and call
[`NetAlgo.calculate`](src/main/java/net/exmachine/algo/distance/NetAlgo.java).

The algorithm was tested on datasets of orbital elements of celestial bodies: meteoroids and
asteroids. The package `net.exmachine.app.orbits` contains a few popular distance functions
between orbits of celestial bodies, an orbital data loading utility class, and the performance
test itself.

To run the speed comparison test, build the project and execute:

```sh
java -jar build/libs/net_algo-1.0-SNAPSHOT.jar conf/sample_conf.yaml
```

Inspect the `conf/sample_conf.yaml` file for algorithm and test parameter descriptions.

## Important files

- `conf` -- various run configurations
- `data` -- orbital elements data files
  - `CAMS-v2-2013.zip` -- meteoroid orbits from the CAMS network (Jenniskens, P., et al. "The
    established meteor showers as observed by CAMS." Icarus 266 (2016): 331-354.)
  - `MPCORB.zip` -- asteroid orbits from the Minor Planet Center (https://www.minorplanetcenter.net/data)
  - `neodys.cat` -- NEAs orbits from the NEODyS-2 catalog (https://newton.spacedys.com/neodys/index.php?pc=5)
- `results` -- results of a test execution

## Brief algorithm description

Given a set of elements and a distance threshold `eps`,
[`NetAlgo`](src/main/java/net/exmachine/algo/distance/NetAlgo.java) finds every pair of elements
whose distance is less than `eps`, without comparing every possible pair.

It first selects a small set of representative "pivot" elements
([`Pivots`](src/main/java/net/exmachine/algo/distance/Pivots.java)) and, around each pivot,
builds a `delta`-neighbourhood: the set of elements within `delta` of that pivot. By the relaxed
triangle inequality, if two elements are within `eps` of each other and one of them lies in a
pivot's neighbourhood, the other element must lie within `K * (delta + eps)` of that same pivot.
This lets the algorithm restrict most pairwise comparisons to elements sharing a pivot's
neighbourhood; only the few elements not covered by any neighbourhood fall back to direct,
pairwise comparison.
