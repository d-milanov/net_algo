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

// The rho_2 orbital distance by Kholshevnikov
// Kholshevnikov K., Kokhirova G., Babadzhanov P., Khamroev U.: Metrics in the space of orbits 
// and their application to searching for celestial objects of
// common origin. Monthly Notices of the Royal Astronomical Society 462(2), 2275-2283 (2016), formula (15)
public class Rho extends Distance
{
    public Rho()
    {
        super( 1.0 );
    }

    public double calcEps2( Body b1, Body b2, double eps2 )
    {
        count++;
        double  sp1p2 = 2.0*b1.sp*b2.sp,
                S2 = b1.p + b2.p - sp1p2*DistUtil.cosI(b1, b2);
        
        if( S2 >= eps2 )
            return Double.POSITIVE_INFINITY;

        S2 += b1.e*b1.e*b1.p + b2.e*b2.e*b2.p 
                - sp1p2*b1.e*b2.e*DistUtil.cosP(b1, b2);

        return S2 < eps2 ? Util.sqrt( S2 ) : Double.POSITIVE_INFINITY;
    }
}