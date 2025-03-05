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
import java.util.*;
import java.io.PrintStream;

// Distancecalculation correctness test program.
// Input: number of bodies (positive integer), max pericentral distance (positive double), 
//        max eccentricity (positive double), test scenario (1 or 2)
public class DistTest
{
    private static final double PRECISION = 1e-4;    

    public static void main( String args[] )
    {
        Locale.setDefault( new Locale( "us" ));

        int N = Integer.parseInt( args[0] );
        double maxQ = Double.parseDouble( args[1] ),
               maxE = Double.parseDouble( args[2] );
        int testN = Integer.parseInt( args[3] );

        if( testN == 1 )
            distTest1( N, maxQ, maxE );
        else if( testN == 2 )
            distTest2( N, maxQ, maxE );
        else
            throw new IllegalArgumentException( "testN is " + testN );
    }


    private static void distTest1( int N, double maxQ, double maxE )
    {
        Body[] bb = Util.randomBodies( N, maxQ, maxE );
        Distance D1 = new D1( maxE ),
                 D2 = new D2(),
                 DD = new D_D( maxE ),
                 DRho = new Rho(),
                 DH = new D_H( maxE );

        PrintStream out = System.out;

        for( int i = 0; i < N; i++ )
        {
            Body bi = bb[i];
            for( int j = i; j < N; j++ )
            {
                Body bj = bb[j];

                double d1 = testDist( bi, bj, D1, out ),
                       d2 = testDist( bi, bj, D2, out ),
                       dd = testDist( bi, bj, DD, out ),
                       dr = testDist( bi, bj, DRho, out ),
                       dh = testDist( bi, bj, DH, out );

                printBody( bi, out );
                out.print(',');
                printBody( bj, out );
                out.printf(",%.10f,%.10f,%.10f,%.10f,%.10f%n", d1, d2, dd, dr, dh);

            }
        }
    }

    private static void distTest2( int N, double maxQ, double maxE )
    {
        Body[] bb = Util.randomBodies( N, maxQ, maxE );
        Distance D1 = new D1( maxE ),
                 D2 = new D2(),
                 DD = new D_D( maxE ),
                 DRho = new Rho(),
                 DH = new D_H( maxE );

        PrintStream out = System.out;

        for( int i = 1; i < N; i++ )
        {
            Body bi = bb[i-1],
                 bj = bb[i];

            double d1 = testDist( bi, bj, D1, out ),
                   d2 = testDist( bi, bj, D2, out ),
                   dd = testDist( bi, bj, DD, out ),
                   dr = testDist( bi, bj, DRho, out ),
                   dh = testDist( bi, bj, DH, out );

            printBody( bi, out );
            out.print(',');
            printBody( bj, out );
            out.printf(",%.10f,%.10f,%.10f,%.10f,%.10f%n", d1, d2, dd, dr, dh);
        }
    }


    private static double testDist( Body b1, Body b2, Distance dist, PrintStream out )
    {
        double d12 = dist.calc( b1, b2 ),
               d21 = dist.calc( b2, b1 ),
               dd = Math.abs( d12 - d21 );

        if( Double.isNaN( dd ) || d12 < 0.0 || d21 < 0.0 || dd > PRECISION || b1 == b2 && d12 > PRECISION )
        {
            String err = String.format( "Error %s:%n%s%n%s%nd12=%.10f d21=%.10f dd=%.10f", dist.getClass().getSimpleName(), b1, b2, d12, d21, dd );
            out.println( err );
            throw new RuntimeException( err );
        }

        double eps = 0.0,
               eps2 = eps*eps,
               d12e = dist.calcEps2( b1, b2, eps2 ),
               d21e = dist.calcEps2( b2, b1, eps2 ),
               dde = Math.abs( d12e - d21e );
/*
        if( !Double.isInfinite( d12e ) || !Double.isInfinite( d21e ))
        {
            String err = String.format( "Error %s:%n%s%n%s%nd12e=%f d21e=%f dde=%f eps=%f", dist.getClass().getSimpleName(), b1, b2, d12e, d21e, dde, eps );
            out.println( err );
            throw new RuntimeException( err );
        }
*/
        for( eps = 0.1; eps < 2.1; eps += 0.2 )
        {
            eps2 = eps*eps;
            d12e = dist.calcEps2( b1, b2, eps2 );
            d21e = dist.calcEps2( b2, b1, eps2 );
            dde = Math.abs( d12e - d21e );

            if( d12e < 0.0 || d21e < 0.0 || 
                b1 == b2 && d12e > PRECISION || 
                d12 < eps -PRECISION && (Math.abs(d12 - d12e) > PRECISION || dde > PRECISION || Double.isNaN( dde )) || 
                d12 > eps + PRECISION && (!Double.isInfinite( d12e ) || !Double.isInfinite( d21e )))
            {
                String err = String.format( "Error %s:%n%s%n%s%nd12e=%.10f d21e=%.10f dde=%.10f eps=%.10f", dist.getClass().getSimpleName(), b1, b2, d12e, d21e, dde, eps );
                out.println( err );
                throw new RuntimeException( err );
            }
        }

        return d12;
    }

    private static void printBody( Body b, PrintStream out )
    {
        out.printf( "%d,%s,%.10f,%.10f,%.10f,%.10f,%.10f", b.N, b.name, b.q, b.e, b.i, b.OM, b.om );
    }
}