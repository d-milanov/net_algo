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

// Net algorithm performance test program.
// Input: configuration file name.
public class Test
{
    private static final boolean CALC = false;

    public static void main( String args[] )
    {
        Locale.setDefault( new Locale( "us" ));

        try
        {
            new Test().perform(new File(args[0]));
        }
        catch( IOException e )
        {
            throw new RuntimeException( e );
        }
    }

    public void perform( File confFile ) throws IOException
    {
        Config conf = Config.load( confFile );

        conf.out.println( "Testing algo" );
        conf.printAlgoParams();
        conf.out.printf( "%nagainst benchmark algo%n" );
        conf.printBenchmarkParams();

        doCheck( conf );
        doTest( conf );
    }

    private void doCheck( Config conf ) throws IOException
    {
        conf.printSection();

        if( !conf.doCheck())
        {
            conf.out.println( "Correctness check is skipped" );
            return;
        }
        else
            conf.out.println( "Correctness check" );

        Algo algo = conf.makeAlgo(),
             benchmark = conf.makeBenchmark();

        DistMatrix aMatr = new DistMatrix(),
                   bMatr = new DistMatrix();

        ArrayList<Body> data = conf.getCheckData();

        doTest( data, algo, aMatr, conf.eps, null );
        doTest( data, benchmark, bMatr, conf.eps, null );

        DistMatrix.DeltaStat delta = DistMatrix.deltaStat( aMatr, bMatr );

        conf.out.print( delta );

//        conf.printAlgoMatrix( aMatr );
//        conf.printBenchmarkMatrix( bMatr );
    }

    private void doTest( Config conf ) throws IOException
    {
        conf.printSection();

        if( !conf.doTest())
        {
            conf.out.println( "Time test is skipped" );
            return;
        }
        else
            conf.out.println( "Test" );

        ArrayList<Body> data = conf.getTestData();
//printData( data, conf.out );

        Stat aStat = new Stat(),
             bStat = new Stat();

        DistMatrix aMatr = null,
                   bMatr = null;


        for( int i = 0; i < conf.nTests; i++ )
        {
            Algo algo = conf.makeAlgo(),
                 benchmark = conf.makeBenchmark();

            aMatr = conf.makeAlgoMatrix();
            bMatr = conf.makeBenchmarkMatrix();

            if( CALC )
            {
                doTest( data, benchmark, bMatr, conf.eps, bStat );
                conf.printBenchmarkMatrix( bMatr );
                return;
            }

            doTest( data, algo, aMatr, conf.eps, aStat );
            doTest( data, benchmark, bMatr, conf.eps, bStat );
        }

        conf.out.println( "Algo stat:" );
        aStat.print( conf.out );

        conf.printSubsection();
        conf.out.println( "Benchmark stat:" );
        bStat.print( conf.out );

        conf.printSubsection();
        conf.out.println( "Avg Benchmark / Avg Algo: " + (bStat.avgTime() / aStat.avgTime()));

        conf.printAlgoMatrix( aMatr );
        conf.printBenchmarkMatrix( bMatr );
    }

    private void doTest( ArrayList<Body> data, Algo algo, DistMatrix matr, double eps, Stat stat )
    {
        if( stat != null )
            stat.algoStarted();

        algo.calculate( data, eps, matr );

        if( stat != null )
            stat.algoEnded( algo );
    }

    private void printData( ArrayList<Body> data, PrintStream out )
    {
        for( Body b : data )
            out.printf( "%f,%f,%f,%f,%f%n", b.q, b.e, b.i, b.OM, b.om );
    }
}