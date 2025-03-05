// The net_algo package contains java implementation of the algorithm for a distance matrix construction.
// This algorithm calculates all distances which do not exceed a fixed number epsilon, for a given set of 
// points and a distance function r(x, y).
// The algorithm works faster than the complete search in the case when the distance
// function r(x, y) complies the relaxed triangle inequality:
//                          r(x, y) <= K(r(x, z) + r(y, z)),
// for any three data points x, y, z, and a constant number K.
//
// The algorithm was applied to datasets of orbital elements of celestial bodies: 
// meteoroids and asteroids. A few popular distance functions between orbits of celestial 
// bodies can be found in the 'distance' subpackage.
//
// Please cite our article 'Relaxed triangle inequality for the orbital similarity criterion by 
// Southworth and Hawkins and its variants' (https://link.springer.com/article/10.1007/s10569-019-9884-6)
// if you use this code for your public work.
//
// Author: Danila Milanov (danila.milanov@gmail.com)
// Year:   2018

package celestial.orbits.net_algo.distance;

import celestial.orbits.net_algo.Body;
import celestial.orbits.net_algo.util.Util;

// D-criterion by Drummond
// Drummond J.: A test of comet and meteor shower associations. Icarus 45(3), 545-553 (1981)
public class D_D extends Distance
{
    // Creates Drummond criterion instance for the given maximal value of an eccentricity E.
    public D_D( double E )
    {
        super( Util.sqrt( 1 + Math.max( 1.0, E*E )));
    }


    public double calcEps2( Body b1, Body b2, double eps2 )
    {
        count++;

        double rq = (b1.q - b2.q) / (b1.q + b2.q),
               S2 = rq*rq;

        if( S2 >= eps2 )
            return Double.POSITIVE_INFINITY;

        double re = (b1.e - b2.e) / (b1.e + b2.e);

        S2 += re*re;

        if( S2 >= eps2 )
            return Double.POSITIVE_INFINITY;

        double I = Util.acos( DistUtil.cosI( b1, b2 )) / Math.PI;

        S2 += I*I;


        if( S2 >= eps2 )
            return Double.POSITIVE_INFINITY;


        double se_2 = (b1.e + b2.e)/2.0,
               P = Util.acos( DistUtil.cosP( b1, b2 )) / Math.PI;

        S2 += se_2*se_2*P*P;

        return S2 < eps2 ? Util.sqrt( S2 ) : Double.POSITIVE_INFINITY;
    }
}