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

package celestial.orbits.net_algo.util;

import java.util.*;
import java.util.regex.*;
import java.io.*;
import celestial.orbits.net_algo.Body;

public class Util
{
    private static final double PRECISION = 1e-6;

    public static void traceTimeDiff( long start, PrintStream out, String name )
    {
        long time = System.currentTimeMillis() - start;

        out.printf( "%s time: %d m %d s %d ms%n", name, time / (60*1000), (time / 1000)%60, time%1000 );
    }

    public static String timeStr( long time )
    {
        return String.format( "%d m %d s %d ms%n", time / (60*1000), (time / 1000)%60, time%1000 );
    }

    public static void traceProgress( int cnt, int total, int progressPercent, PrintStream out )
    {
        int r = (total*progressPercent)/100;

        if( cnt%r == 0 )
            out.printf( "%d%% complete%n", (cnt*progressPercent)/r);
    }

    public static boolean traceProgress( int cnt, int total, int progressPercent, PrintStream out, long start )
    {
        int r = (total*progressPercent)/100;

        if( cnt%r == 0 )
        {
            long time = System.currentTimeMillis() - start;
            out.printf( "%d%% complete in %d m %d s %d ms%n", (cnt*progressPercent)/r, time / (60*1000), (time / 1000)%60, time%1000 );

            return true;
        }

        return false;
    }

    public static Body[] randomBodies( int n, double maxQ, double maxE )
    {
        Body[] res = new Body[n];
        Random r = new Random();

        for( int i = 0; i < n; i++ )
        {
            res[i] = Body.fromRadiansQ( i, String.valueOf( i ), 
                                        maxQ*r.nextDouble(),
                                        maxE*r.nextDouble(),
                                        Math.PI*r.nextDouble(),
                                        2.0*Math.PI*r.nextDouble(),
                                        2.0*Math.PI*r.nextDouble());
        }

        return res;
    }

    public static double acos( double x )
    {
        if( x >= 1 )
        {
            x -= 1.0;

            if( x > PRECISION )
                throw new RuntimeException( "acos: " + x );

            return 0.0;
        }
        else if( x <= -1.0 )
        {
            x += 1.0;

            if( x < -PRECISION )
                throw new RuntimeException( "acos: " + x );

            return Math.PI;
        }

        return Math.acos( x );
    }

    public static double asin( double x )
    {
        if( x >= 1 )
        {
            x -= 1.0;

            if( x > PRECISION )
                throw new RuntimeException( "asin: " + x );

            return Math.PI/2.0;
        }
        else if( x <= -1.0 )
        {
            x += 1.0;

            if( x < -PRECISION )
                throw new RuntimeException( "asin: " + x );

            return -Math.PI/2.0;
        }

        return Math.asin( x );
    }

    public static double sqrt( double x )
    {
        if( x <= 0.0 )
        {
            if( x < -PRECISION )
                throw new RuntimeException( "sqrt: " + x );

            return 0.0;
        }

        return Math.sqrt( x );
    }
}
