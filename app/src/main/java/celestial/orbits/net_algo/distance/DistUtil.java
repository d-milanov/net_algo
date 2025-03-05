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

class DistUtil
{
    static double cosI( Body b1, Body b2 )
    {
        return b1.cosI * b2.cosI + b1.sinI * b2.sinI * cos_d_OM( b1, b2 );
    }

    static double cosP( Body b1, Body b2 )
    {
        double s1 = b1.sinI,
               s2 = b2.sinI,
               c1 = b1.cosI,
               c2 = b2.cosI,
               so1 = b1.sinom,
               so2 = b2.sinom,
               co1 = b1.cosom,
               co2 = b2.cosom;

        return s1*s2*so1*so2 +
               (co1*co2 + c1*c2*so1*so2)*cos_d_OM( b1, b2 ) +
               (c2*co1*so2 - c1*so1*co2)*sin_d_OM( b1, b2 );
    }

    static double cos_d_OM( Body b1, Body b2 )
    {
        return b1.cosOM * b2.cosOM + b1.sinOM * b2.sinOM;
    }

    static double sin_d_OM( Body b1, Body b2 )
    {
        return b1.sinOM * b2.cosOM - b1.cosOM * b2.sinOM;
    }

    static double cos_d_om( Body b1, Body b2 )
    {
        return b1.cosom * b2.cosom + b1.sinom * b2.sinom;
    }

    static double sin_d_om( Body b1, Body b2 )
    {
        return b1.sinom * b2.cosom - b1.cosom * b2.sinom;
    }

    static double cos_sum_i( Body b1, Body b2 )
    {
        return b1.cosI * b2.cosI - b1.sinI * b2.sinI;
    }
}
