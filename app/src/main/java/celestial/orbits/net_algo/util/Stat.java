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
import java.io.*;
import celestial.orbits.net_algo.*;

public class Stat
{
    private Map<Integer, Long> results = new TreeMap<>();
    private int count;
    private long start;
    private String algoMessage;

    public void print( PrintStream out )
    {
        if( results.isEmpty())
        {
            out.println( "No results" );
            return;
        }

        String[] times = cacTimes();

        out.printf( "%nAverage: %s\tBest: %s\tWorst: %s%n", 
                     times[0], times[1], times[2] );

        if( algoMessage != null )
            out.printf( "%nAlgo message:%n%s%n", algoMessage );

        if( results.size() <= 3 )
            return;

        out.printf( "%nAll:%n" );

        for( Map.Entry<Integer, Long> nt : results.entrySet())
            out.printf( "%d: %s", nt.getKey(), Util.timeStr( nt.getValue()));
    }

    public void algoStarted()
    {
        if( start > 0 )
            throw new IllegalStateException();

        start = System.currentTimeMillis();
    }

    public void algoEnded( Algo algo )
    {
        long end = System.currentTimeMillis();

        if( start <= 0 )
            throw new IllegalStateException();

        results.put( ++count, end - start );
        start = 0;

        algoMessage = algo.getStatMessage();
    }

    private String[] cacTimes()
    {
        long max = Long.MIN_VALUE,
             min = Long.MAX_VALUE,
             avg = 0l;

        for( Map.Entry<Integer, Long> nt : results.entrySet())
        {
            long t = nt.getValue();

            if( t > max )
                max = t;
            if( t < min )
                min = t;

            avg += t;
        }

        return new String[]{ Util.timeStr( avg/results.size()), 
                             Util.timeStr( min ), 
                             Util.timeStr( max )};
    }

    public double avgTime()
    {
        long avg = 0l;

        for( Long t : results.values())
            avg += t;

        return avg/(double)results.size();
    }
}