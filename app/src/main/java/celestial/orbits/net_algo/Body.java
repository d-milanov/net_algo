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

package celestial.orbits.net_algo;

import java.util.*;
import celestial.orbits.net_algo.util.Util;


public class Body
{
    public final String name;
    public final int N;
    public final double a;
    public final double q;
    public final double p;
    public final double e;
    public final double i;
    public final double OM;
    public final double om;

    public final double sinI;
    public final double cosI;

    public final double sinom;
    public final double cosom;

    public final double sinOM;
    public final double cosOM;

    public final double sp;



    private Body( int N, String name, double a, double e, double i, double OM, double om )
    {
        this( N, name, a, a*(1.0 - e), e, i, OM, om );
    }

    private Body( int N, String name, double a, double q, double e, double i, double OM, double om )
    {
        this.N = N;
        this.name = name; 
        this.a = a;
        this.q = q;
        this.e = e;
        this.i = i;
        this.OM = OM;
        this.om = om;

        this.sinI = Math.sin( i );
        this.cosI = Math.cos( i );
        this.sinom = Math.sin( om );
        this.cosom = Math.cos( om );
        this.sinOM = Math.sin( OM );
        this.cosOM = Math.cos( OM );


        this.p = q*(1.0 + e);
        this.sp = Util.sqrt( p );
    }

    public static Body fromDegrees( int N, String name, double a, double e, double i, double OM, double om )
    {
        return new Body( N, name, a, e, Math.toRadians( i ), Math.toRadians( OM ), Math.toRadians( om ));
    }

    public static Body fromDegreesQ( int N, String name, double q, double e, double i, double OM, double om )
    {
        return new Body( N, name, Double.NaN, q, e, Math.toRadians( i ), Math.toRadians( OM ), Math.toRadians( om ));
    }

    public static Body fromRadians( int N, String name, double a, double e, double i, double OM, double om )
    {
        return new Body( N, name, a, e, i, OM, om );
    }

    public static Body fromRadiansQ( int N, String name, double q, double e, double i, double OM, double om )
    {
        return new Body( N, name, Double.NaN, q, e, i, OM, om );
    }

    public String toString()
    {
        return String.format( "%d %s %f(%f) %f %f %f %f", N, name, a, q, e, i, OM, om );
    }

    public static Comparator<Body> N_COMPARATOR = new Comparator<Body>()
    {
        public int compare( Body b1, Body b2 )
        {
            return b1.N - b2.N;
        }
    };

}