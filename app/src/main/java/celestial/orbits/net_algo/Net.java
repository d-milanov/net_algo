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

// The algorithm for the fast construction of a distance matrix.
public class Net implements Algo
{
    private final double delta;
    private int[][] adjacent;
    private final ArrayList<TreeSet<Body>> neighbourhoods = new ArrayList<>();
    private final TreeSet<Body> heap = new TreeSet<>( Body.N_COMPARATOR );
    private final Distance dist;
    private final int pivotStep;
    private final int maxNodes;

    private final PrintStream dbgOut;
    private int bodiesSize;
    private final int progress;

    public Net( Map<String, String> params )
    {
        dist = Distance.create( params );
        dbgOut = Config.getOutput( params.get( "dbgOut" ));
        progress = Config.getInt( params.get( "progress" ), 0 );
        pivotStep = Config.getInt( params.get( "pivotStep" ));
        delta = Config.getDouble( params.get( "delta" ));
        maxNodes = Config.getInt( params.get( "maxNodes" ));
    }


    public void calculate( ArrayList<Body> bodies, double eps, DistMatrix matr )
    {
        this.bodiesSize = bodies.size();
        long t0 = System.currentTimeMillis();

        init( bodies, eps );

        int len = bodies.size();

        for( Body b : bodies )
        {
            iterateNeighbours( b, eps, matr );

            if( dbgOut != null && progress != 0 )
                Util.traceProgress( b.N, len, progress, dbgOut, t0 );
        }
    }

    private void init( ArrayList<Body> bodies, double eps )
    {
        long t0 = System.currentTimeMillis();

        Body[] nodes = Pivots.selectPivotsMatrDel( maxNodes, bodies, pivotStep, dist, delta );

        this.adjacent = new int[bodies.size()][];

        final int len1 = nodes.length + 1;
        final int len2 = nodes.length + 2;
        final double th = dist.K * (delta + eps);
        final double th2 = th*th;
        
        heap.addAll( bodies );

        int nCnt = 0;

        for( Body n : nodes )
        {
            TreeSet<Body> neighbourhood = new TreeSet<>( Body.N_COMPARATOR );


            for( Body bn : bodies )
            {
                int adj[] = adjacent[bn.N];
                int ind = 0;

                if( adj == null )
                {
                    adj = adjacent[bn.N] = new int[len2];
                    adj[0] = -1;
                }
                else
                    while( adj[ind] >= 0 ) ind++;
                
                
                double d = dist.calcEps2( n, bn, th2 );

                if( d < th )
                {
                    adj[ind] = nCnt;
                    adj[ind + 1] = -1;

                    if( d < delta && adj[len1] == 0 )
                    {
                        neighbourhood.add( bn );
                        adj[len1] = 1;
                    }
                }
            }

            nCnt++;
            neighbourhoods.add( neighbourhood );
            heap.removeAll( neighbourhood );

            if( dbgOut != null )
                dbgOut.println( "Node " + n.name + ": " + neighbourhood.size() );
        }

        if( dbgOut != null )
        {
            dbgOut.println( "Heap: " + heap.size());
/*
            dbgOut.println( "Nodes adjacent:" );
            for( Body n : nodes )
                dbgOut.println( Arrays.toString( adjacent[n.N] ));
*/
            Util.traceTimeDiff( t0, dbgOut, "Net.ctor" );
        }
    }


    void iterateNeighbours( Body b, double eps, DistMatrix matr )
    {
        final double eps2 = eps*eps;

        int[] adj = adjacent[b.N];

        for( int i : adj )
        {
            if( i < 0 )
                break;

            Set<Body> bns = neighbourhoods.get( i ).tailSet( b, false );

            for( Body bn : bns )
            {
                double d = dist.calcEps2( b, bn, eps2 );

                if( d < eps )
                    matr.add( b, bn, d );
            }
        }

        for( Body bn : heap.tailSet( b, false ))
        {
            double d = dist.calcEps2( b, bn, eps2 );

            if( d < eps )
                matr.add( b, bn, d );
        }
            
    }

    public String getStatMessage()
    {
        return String.format( "%d (%f) distance calculations", dist.count, 2.0 * dist.count / (((double)bodiesSize) * (bodiesSize - 1d)));
    }
}