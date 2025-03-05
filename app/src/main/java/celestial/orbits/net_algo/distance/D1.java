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

// D-criterion by Southworth & Hawkins
// Southworth R., Hawkins G.: Statistics of meteor streams. Smithsonian Contributions
// to Astrophysics 7, 261 (1963)
public class D1 extends Distance
{
    // Creates Southworth & Hawkins criterion instance for the given maximal value of an eccentricity E.
    public D1( double E )
    {
        super( Util.sqrt( 2 + 4*E*E ));
    }

    public double calcEps2( Body b1, Body b2, double eps2 )
    {
        count++;

        double dq = b1.q - b2.q,
               S2 = dq*dq;

        if( S2 >= eps2 )
            return Double.POSITIVE_INFINITY;

        double de = b1.e - b2.e;

        S2 += de*de;

        if( S2 >= eps2 )
            return Double.POSITIVE_INFINITY;


        double cosI = DistUtil.cosI( b1, b2 ),
               sin_2_I_2 = (1.0 - cosI) /2.0;

        S2 += 4.0*sin_2_I_2;

        if( S2 >= eps2 )
            return Double.POSITIVE_INFINITY;


        double se = b1.e + b2.e,
               dOM = b1.OM - b2.OM,
               cos_delta = DistUtil.cos_d_om( b1, b2 ),
               sin_delta = DistUtil.sin_d_om( b1, b2 ),
               ksi;

        if( sin_2_I_2 < 1.0 - 1e-5 )
        {
            double t = ( 1.0 + DistUtil.cos_sum_i( b1, b2 )) * 
                        ( 1.0 - DistUtil.cos_d_OM( b1, b2 )) / (1.0 - sin_2_I_2 );

            if( t < -1e-4 )
                throw new RuntimeException( "t < 0: " + t + " bodies: " + b1 + "; " + b2 );
            
            ksi = t > 0 ? Util.sqrt( t ) / 2.0 : 0.0;
        }
        else
            ksi = Math.abs( DistUtil.sin_d_OM( b1, b2 ));

        double ksi2 = ksi*ksi;

        if( b1.i + b2.i > Math.PI )
            ksi = -ksi;

        if( dOM < 0.0 )
            ksi = -ksi;

        if( Math.abs( dOM ) > Math.PI )
            ksi = -ksi;

        double sin_2_P_2 = 0.5 * 
                (1.0 - cos_delta * (1.0 - 2.0*ksi2) + 2.0*sin_delta*ksi*Util.sqrt(1.0 - ksi2));

        S2 += se*se*sin_2_P_2;

        return S2 < eps2 ? Util.sqrt( S2 ) : Double.POSITIVE_INFINITY;
    }
}