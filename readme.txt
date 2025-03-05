The net_algo package contains java implementation of the algorithm presented in the article

Relaxed triangle inequality for the orbital similarity criterion by Southworth and Hawkins and its variants https://link.springer.com/article/10.1007/s10569-019-9884-6
Please cite it if you use this code in your public work.

The algorithm calculates all distances which do not exceed a fixed number epsilon, for a given set of points and a distance function r(x, y).
The algorithm works faster than the complete search in the case when the distancefunction r(x, y) complies the relaxed triangle inequality:
                         r(x, y) <= K(r(x, z) + r(y, z)),
for any three data points x, y, z, and a constant number K.

The algorithm was applied to datasets of orbital elements of celestial bodies: meteoroids and asteroids. 
A few popular distance functions between orbits of celestial bodies can be found in the 'distance' subpackage.

To build:
gradle build

To test:
run.sh 
or
run.bat
or
java  -cp app/build/libs/app.jar  celestial.orbits.net_algo.Test <your config>


Important files
 - app/src/main/java -- the algorithm java code
 - app/src/main/python -- python test for distance calculation 
 - conf -- various run configurations
 - data -- orbital elements data files (you may need to unzip them)
  - CAMS-v2-2013.csv -- meteoroid orbits from the CAMS network (Jenniskens, P., et al. "The established meteor showers as observed by CAMS." Icarus 266 (2016): 331-354.)
  - MPCORB.DAT -- asteroid orbits from the Minor Planet Center (https://www.minorplanetcenter.net/data)
  - neodys.cat -- NEAs orbits from the NEODyS-2 catalog (https://newton.spacedys.com/neodys/index.php?pc=5)
 - results -- resuts of a test execution 

The algorithm code can be found in the Net.java. Pivots selection algorithm is in the Pivots.java.
The distance source files folder contains various orbital distance functions implementation.
The Test.java encodes a program that tests the Net algorithm performance against the complete
search approach, implemented in the BruteForce.java.
The python code includes oribital distance functions implementation. The python functions were 
used as a test for the java distance calculation.

The algorithm, briefly:

Firstly, the distance matrix is built for the reduced sample (each 20-th element
is taken) with the complete search algorithm. A collection of pivots is
selected from the most dense regions of the sample, found from the reduced
sample matrix approximation. A delta-neighbourhoods were built for all pivots,
using a predefined parameter delta. For each element x of the sample a subset of
'adjacent' pivots p_1,... p_n(x) is defined, so that rho(x, p_i) <= K(epsilon + delta). Then,
the full proximity matrix is constructed as follows. For each x of the sample,
only the distances to the neighborhood points of adjacent pivots and to the
'heap' of points, not included into any neighbourhood are calculated.
