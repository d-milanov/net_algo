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

import java.util.*;
import celestial.orbits.net_algo.Body;
import celestial.orbits.net_algo.util.Util;

// Distance root class. 
public abstract class Distance
{
    // Constant in a relaxed triangle inequality
    public final double K;

    // Calculations count for performance analysis
    public long count;

    Distance( double K )
    {
        this.K = K;
    }

    // Returns the distance between two bodies if it's square is less than eps2.
    // Otherwise returns Double.POSITIVE_INFINITY.
    public abstract double calcEps2( Body b1, Body b2, double eps2 );

    double calc( Body b1, Body b2 )
    {
        return calcEps2( b1, b2, Double.POSITIVE_INFINITY );
    }

    // Creates Distance instance by name and parameters
    public static Distance create( Map<String, String> params )
    {
        String name = params.get( "dist" ),
               E = params.get( "dist.E" );

        if( name == null )
            throw new IllegalArgumentException( "dist config parameter is not defined" );

        if( E == null && (name.equals( "D1" ) || name.equals( "D_D" ) || name.equals( "D_H" )))
            throw new IllegalArgumentException( "distance " + name + " requires E config parameter to be defined" );

        switch( name )
        {
            case "D1": return new D1( Double.parseDouble( E ));
            case "D2": return new D2();
            case "DD":
            case "D_D": return new D_D( Double.parseDouble( E ));
            case "Rho": return new Rho();
            case "DH":
            case "D_H": return new D_H( Double.parseDouble( E ));
            case "Eucl": return new Eucl();
        }

        throw new IllegalArgumentException( "Unknown distance: " + name );
    }
}