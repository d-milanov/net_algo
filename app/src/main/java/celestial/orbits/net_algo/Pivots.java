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
import java.util.regex.*;
import java.io.*;
import celestial.orbits.net_algo.distance.*;
import celestial.orbits.net_algo.util.*;


// Various pivots selection algorithms.
public class Pivots
{
    static Body[] selectPivotsRnd( int n, List<Body> bodies )
    {
        Random r = new Random();
        Body[] res = new Body[n];

        for( int i = 0; i < n; i++ )
        {
            Body b = bodies.get( r.nextInt( bodies.size()));
            res[i] = b;

            for( int j = 0; j < i; j++ )
                if( res[j] == b )
                {
                    i--;
                    break;
                }

        }

        return res;
    }

    static Body[] selectPivots1( List<Body> bodies )
    {
        return new Body[]{ Body.fromRadians( -1, "pivot1", 0.01, 0.0, Math.PI/4.0, 0.0, 0.0 )};
    }

    static Body[] selectPivotsFirst( int n, List<Body> bodies )
    {
        Body[] res = new Body[n];

        for( int i = 0; i < n; i++ )
        {
            Body b = bodies.get( i );
            res[i] = b;
        }

        return res;
    }

    static Body[] selectPivotsMean( List<Body> bodies )
    {
        double ma = 0.0,
               me = 0.0,
               mi = 0.0,
               mO = 0.0,
               mo = 0.0,
               nb = (double)bodies.size();

        for( Body b : bodies )
        {
            ma += b.a;
            me += b.e;
            mi += b.i;
            mO += b.OM;
            mo += b.om;
        }

        return new Body[]{ Body.fromRadians( -1, "pivotM", ma/nb, me/nb, mi/nb, mO/nb, mo/nb )};
    }

    static Body[] selectPivotsMatrDel( int n, ArrayList<Body> bodies, int step, Distance dist, double delta )
    {
        return selectPivotsMatrDel( n, bodies, step, dist, delta, 0.0 );
    }

    static Body[] selectPivotsMatrDel( int n, ArrayList<Body> bodies, int step, Distance dist, double delta, double minD )
    {
        Map<Body, List<Body>> neighbours = new HashMap<>();
        Map<Body, Integer> bodySz = new HashMap<>();

        double delta2 = delta*delta;

        for( int i = 0, len = bodies.size(); i < len; i += step )
        {
            Body b1 = bodies.get( i );

            for( int j = i + step; j < len; j += step )
            {
                Body b2 = bodies.get( j );

                double d = dist.calcEps2( b1, b2, delta2 );

                if( d < delta )
                    addBodies( b1, b2, neighbours, bodySz );
            }
        }

        NavigableMap<Integer, Set<Body>> neighbSz = makeNeighbSz( bodySz );
        Body[] res = new Body[n];
        Set<Body> removed = new HashSet<>();
        int cnt = 0;

        while( !neighbSz.isEmpty())
        {
            Map.Entry<Integer, Set<Body>> nSzbb = neighbSz.lastEntry();
            Set<Body> bb = nSzbb.getValue();
            int bbSz = nSzbb.getKey();

            Body b = null;
            for( Iterator<Body> i = bb.iterator(); i.hasNext(); )
            {
                Body b1 = i.next();

                i.remove();

                if( !removed.contains( b1 ))
                {
                    b = b1;
                    break;
                }
            }

            if( bb.isEmpty())
            {
                neighbSz.remove( bbSz );

                if( b == null )
                    continue;
            }

            List<Body> bNs = neighbours.get( b );

            bNs.removeAll( removed );

            int bNsz = bNs.size();

            if( bNsz < bbSz )
            {
                if( bNsz > 0 )
                    putSzB( bNsz, b, neighbSz, bodySz );
            }
            else if( checkMinDist( b, res, cnt, minD, dist ))
            {
                res[cnt++] = b;
//System.out.println( "Pivot found: " + bNsz + " - " + b.name );
//System.out.println(  b );
                if( cnt == n )
                    break;

                removed.addAll( bNs );
            }
        }

        if( cnt < n )
        {
            Body[] r = new Body[cnt];
            System.arraycopy( res, 0, r, 0, cnt );
            res = r;
        }

        return res;
    }

    private static boolean checkMinDist( Body b, Body[] bb, int len, double minD, Distance dist )
    {
        if( minD < 1E-10 )
            return true;

        double minD2 = minD*minD;

        for( int i = 0; i < len; i++ )
        {
            Body b1 = bb[i];

            if( dist.calcEps2( b, b1, minD2 ) < minD )
                return false;
        }

        return true;
    }

    private static void putSzB( int sz, Body b, NavigableMap<Integer, Set<Body>> neighbSz, Map<Body, Integer> bodySz )
    {
        Set<Body> bb = neighbSz.get( sz );

        if( bb == null )
        {
            bb = new HashSet<>();
            neighbSz.put( sz, bb );
        }

        bb.add( b );
        bodySz.put( b, sz );
    }

    private static void addBodies( Body b1, Body b2, Map<Body, List<Body>> neighbours, Map<Body, Integer> bodySz )
    {
        addBodies0( b1, b2, neighbours, bodySz );
        addBodies0( b2, b1, neighbours, bodySz );
    }

    private static void addBodies0( Body b, Body b1, Map<Body, List<Body>> neighbours, Map<Body, Integer> bodySz )
    {
        List<Body> bns = neighbours.get( b );

        if( bns == null )
        {
            bns = new LinkedList<>();
            neighbours.put( b, bns );
        }

        bns.add( b1 );
        bodySz.put( b, neighbours.size());
    }

    private static NavigableMap<Integer, Set<Body>>  makeNeighbSz( Map<Body, Integer> bodySz )
    {
        NavigableMap<Integer, Set<Body>> res = new TreeMap<>();

        for( Map.Entry<Body, Integer> bsz : bodySz.entrySet())
        {
            Body b = bsz.getKey();
            int sz = bsz.getValue();

            Set<Body> bb = res.get( sz );

            if( bb == null )
            {
                bb = new HashSet<>();
                res.put( sz, bb );
            }

            bb.add( b );
        }

        return res;
    }
}