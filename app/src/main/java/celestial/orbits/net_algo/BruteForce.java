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
import java.io.*;

import celestial.orbits.net_algo.distance.*;
import celestial.orbits.net_algo.util.*;



public class BruteForce implements Algo
{
    private final PrintStream dbgOut;
    private final int progress;
    private final Distance dist;

    private int bodiesSize;

    public BruteForce( Map<String, String> params )
    {
        dist = Distance.create( params );
        dbgOut = Config.getOutput( params.get( "dbgOut" ));
        progress = Config.getInt( params.get( "progress" ), 0 );
    }


    public void calculate( ArrayList<Body> bodies, double eps, DistMatrix matr )
    {
        this.bodiesSize = bodies.size();
        long t0 = System.currentTimeMillis();
        double eps2 = eps*eps;

        for( int i = 0, len = bodies.size(); i < len; i++ )
        {
            Body bi = bodies.get( i );

            for( int j = i + 1; j < len; j++ )
            {
                Body bj = bodies.get( j );

                double d = dist.calcEps2( bi, bj, eps2 );

                if( d < eps )
                    matr.add( bi, bj, d );
            }

            if( dbgOut != null && progress != 0 )
                Util.traceProgress( i, len, progress, dbgOut, t0 );
        }
    }

    public String getStatMessage()
    {
        return String.format( "%d (%f) distance calculations", dist.count, 2.0 * dist.count / (((double)bodiesSize) * (bodiesSize - 1d)));
    }
}