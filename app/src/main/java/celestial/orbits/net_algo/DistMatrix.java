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

// Distance matrix: the algorithm execution result storage class. 
// See Stub and Distrib inner classes for implementation variants.
public class DistMatrix
{
    final Map<Body, Map<Body, Double>> matrix = new TreeMap<>( Body.N_COMPARATOR );

    public void add( Body b1, Body b2, double d )
    {
        if( b1.N > b2.N )
        {
            Body b = b2;
            b2 = b1;
            b1 = b;
        }

        Map<Body, Double> b1bds = matrix.get( b1 );

        if( b1bds == null )
        {
            b1bds = new TreeMap<>( Body.N_COMPARATOR );

            matrix.put( b1, b1bds );
        }

        b1bds.put( b2, d );
    }

    public String toString()
    {
        StringBuilder buff = new StringBuilder();

        for( Map.Entry<Body, Map<Body, Double>> bbd : matrix.entrySet())
        {
            Body b1 = bbd.getKey();

            buff.append( b1.name ).append( String.format( "%n" ));

            for( Map.Entry<Body, Double> bd : bbd.getValue().entrySet())
            {
                Body b2 = bd.getKey();

                buff.append( '\t')
                    .append( b2.name )
                    .append( "," )
                    .append( String.format( "%.15f%n", bd.getValue()));
            }
        }

        return buff.toString();
    }

    public void print( PrintStream out )
    {
        for( Map.Entry<Body, Map<Body, Double>> bbd : matrix.entrySet())
        {
            Body b1 = bbd.getKey();

            out.println( b1.name );

            for( Map.Entry<Body, Double> bd : bbd.getValue().entrySet())
            {
                Body b2 = bd.getKey();
                out.printf( "\t%s,%.15f%n", b2.name, bd.getValue());
            }
        }
    }

    public long size()
    {
        long res = 0;

        for( Map<Body, Double> bd : matrix.values())
            res += bd.size();

        return res;
    }


    public static class Stub extends DistMatrix
    {
        long count;

        public void add( Body b1, Body b2, double d )
        {
            count++;
        }

        public String toString()
        {
            return "Size: " + String.valueOf( count );
        }

        public void print( PrintStream out )
        {
            out.println( this );
        }
    }

    public static class Distrib extends DistMatrix
    {
        final double precision;
        Map<Long, Long> distr = new TreeMap<>();

        public Distrib()
        {
            precision = 1e-6;
        }

        public void add( Body b1, Body b2, double d )
        {
            if( d < precision )
                d = precision;

            long logD = Math.round( Math.log10( d ) / precision);

            Long dCnt = distr.get( logD );

            if( dCnt == null )
                dCnt = 0l;

            distr.put( logD, ++dCnt );
        }

        public String toString()
        {
            return "Size: " + String.valueOf( distr.size());
        }

        public void print( PrintStream out )
        {
            for( Map.Entry<Long, Long> dc : distr.entrySet())
                out.printf( "%f,%d%n", (double)dc.getKey() * precision, dc.getValue());
        }
    }

    public static DeltaStat deltaStat( DistMatrix m1, DistMatrix m2 )
    {
        long m1_m2 = 0;
        double ddMax = 0d;

        for( Map.Entry<Body, Map<Body, Double>> bbds : m1.matrix.entrySet())
        {
            Body b1 = bbds.getKey();
            Map<Body, Double> bds1 = bbds.getValue(),
                              bds2 = m2.matrix.get( b1 );

            if( bds2 == null )
            {
                m1_m2 += bds1.size();
                continue;
            }

            for( Map.Entry<Body, Double> bd1 : bds1.entrySet())
            {
                Double d2 = bds2.get( bd1.getKey());

                if( d2 == null )
                {
                    m1_m2++;
                    continue;
                }

                double d1 = bd1.getValue(),
                       dd = Math.abs( d1 - d2 );

                if( dd > ddMax )
                    ddMax = dd;
            }
        }

        long m1Sz = m1.size(),
             m2Sz = m2.size(),
             m2_m1 = m2Sz - (m1Sz - m1_m2);

        return new DeltaStat( m1Sz, m2Sz, m1_m2, m2_m1, ddMax );
    }

    public static class DeltaStat
    {
        final long m1_m2;
        final long m2_m1;
        final double dMax;
        final long m1;
        final long m2;

        public DeltaStat( long m1, long m2, long m1_m2, long m2_m1, double dMax )
        {
            this.m1_m2 = m1_m2;
            this.m2_m1 = m2_m1;
            this.dMax = dMax;
            this.m1 = m1;
            this.m2 = m2;
        }

        public String toString()
        {
            return String.format( "#M1 = %d  #M2 = %d  #(M1 \\ M2) = %d  #(M2 \\ M1) = %d  max|d1 - d2| = %f", m1, m2, m1_m2, m2_m1, dMax );
        }
    }
}