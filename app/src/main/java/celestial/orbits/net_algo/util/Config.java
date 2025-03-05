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

public class Config
{
    public static Config load( File f ) throws IOException
    {
        BufferedReader r = null;
        Map<String, String> params = new TreeMap<>();

        try
        {
            r = new BufferedReader( new InputStreamReader( new FileInputStream( f )));

            String l;

            while( (l = r.readLine()) != null )
            {
                int cmnt = l.indexOf( '#' );

                if( cmnt >= 0 )
                    l = l.substring( 0, cmnt );

                l = l.trim();

                if( l.isEmpty())
                    continue;

                int eq = l.indexOf( "=" );

                if( eq < 0 )
                {
                    System.out.println( "Ignoring config line: " + l );
                    continue;
                }

                String key = l.substring( 0, eq ).trim(),
                       value = eq < l.length() - 1 ? l.substring( eq + 1 ).trim() : null;

                params.put( key, value );
            }
        }
        finally
        {
            if( r != null )
                r.close();
        }

        return new Config( params );
    }

    private final Map<String, String> algo = new TreeMap<>();
    private final Map<String, String> benchmark = new TreeMap<>();
    private final Map<String, String> check = new TreeMap<>();
    private final Map<String, String> test = new TreeMap<>();

    public double eps = 0.2;
    public int nTests = 1;
    public PrintStream out = System.out;
    public PrintStream dbgOut = System.out;

    private Config( Map<String, String> params )
    {
        for( Map.Entry<String, String> kv : params.entrySet())
        {
            String k = kv.getKey(),
                   v = kv.getValue();

            try
            {
                if( k.startsWith( "check." ))
                    check.put( k.substring( "check.".length()), v );
                else if( k.startsWith( "test." ))
                    test.put( k.substring( "test.".length()), v );
                else if( k.startsWith( "algo." ))
                    algo.put( k.substring( "algo.".length()), v );
                else if( k.startsWith( "benchmark." ))
                    benchmark.put( k.substring( "benchmark.".length()), v );
                else
                {
                    switch( k )
                    {
                        case "algo": algo.put( "class", v ); break;
                        case "benchmark": benchmark.put( "class", v ); break;
                        case "eps": eps = Double.parseDouble( v ); break;
                        case "nTests": nTests = Integer.parseInt( v ); break;
                        case "dbgOut": dbgOut = getOutput( v );
                        case "out": out = getOutput( v );
                    }
                }
            }
            catch( Exception e )
            {
                throw new RuntimeException( "Config creation error at " + k + " = " + v, e );
            }
        }
    }

    public Algo makeAlgo()
    {
        return makeAlgo( algo.get( "class" ), algo );
    }

    public Algo makeBenchmark()
    {
        return makeAlgo( benchmark.get( "class" ), benchmark );
    }

    private Algo makeAlgo( String cls, Map<String, String> params )
    {
        switch( cls )
        {
            case "BruteForce" : return new BruteForce( params );
            case "Net" : return new Net( params );
        }

        throw new IllegalArgumentException( "Failed to instantiate algo " + cls );
    }

    public ArrayList<Body> getCheckData() throws IOException
    {
        return getData( check );
    }

    public ArrayList<Body> getTestData() throws IOException
    {
        return getData( test );
    }

    public ArrayList<Body> getData( Map<String, String> params ) throws IOException
    {
        String lName = params.get( "dataLoader" );

        if( lName.equals( "Random" ))
        {
            Body[] bds = Util.randomBodies( getInt( params.get( "dataLoader.n" )), 
                                            getDouble( params.get( "dataLoader.maxQ" )),
                                            getDouble( params.get( "dataLoader.maxE" )));

            ArrayList<Body> res = new ArrayList<>();
            for( Body b : bds )
                res.add( b );
            return res;
        }
        else
        {
            File f = new File( params.get( "dataFile" ));
            DataLoader l = newInst( lName, "Failed to load data");
            DataLoader.Filter filter = createFilter( params );
            return l.load( f, filter );
        }
    }

    public DistMatrix makeAlgoMatrix()
    {
        return newInst( algo.get( "matrix" ), "Failed to make algo matrix" );
    }

    public DistMatrix makeBenchmarkMatrix()
    {
        return newInst( benchmark.get( "matrix" ), "Failed to make benchmark matrix" );
    }

    private <T> T newInst( String cls, String msg )
    {
        if( cls == null )
            return null;

        try
        {
            return (T)Class.forName( cls ).newInstance();
        }
        catch( Exception e )
        {
            throw new RuntimeException( msg + " " + cls, e );
        }
    }

    private DataLoader.Filter createFilter( Map<String, String> params )
    {
        String name = params.get( "dataFilter" );

        if( name == null )
            return null;

        switch( name )
        {
            case "AFilter": 
                return new DataLoader.AFilter( 
                        Double.parseDouble( params.get("dataFilter.a_min" )),
                        Double.parseDouble( params.get("dataFilter.a_max" )));
            case "EFilter": 
                return new DataLoader.EFilter( 
                        Double.parseDouble( params.get("dataFilter.e_min" )),
                        Double.parseDouble( params.get("dataFilter.e_max" )));
            case "Reduce": 
                return new DataLoader.Reduce( 
                        Long.parseLong( params.get( "dataFilter.times" )));
            case "MainBelt": 
                return new DataLoader.MainBelt( params.containsKey( "dataFilter.times" ) ?
                     Long.parseLong( params.get( "dataFilter.times" )) : 1l);
        }

        return newInst( name,  "Failed to instantiate data filter " + name );
    }


    public boolean doCheck()
    {
        return check.containsKey( "dataFile" ) && check.containsKey( "dataLoader" );
    }

    public boolean doTest()
    {
        return test.containsKey( "dataFile" ) && test.containsKey( "dataLoader" );
    }

    public void printAlgoParams()
    {
        printParams( algo );
    }

    public void printBenchmarkParams()
    {
        printParams( benchmark );
    }

    private void printParams( Map<String, String> pp )
    {
        out.print( pp );
    }

    public void printSection()
    {
        out.printf( "%n==========================================================%n" );
    }
    public void printSubsection()
    {
        out.printf( "%n=============================%n" );
    }

    public static double getDouble( String d )
    {
        return Double.parseDouble( d );
    }

    public static int getInt( String i )
    {
        return Integer.parseInt( i );
    }

    public static int getInt( String i, int dflt )
    {
        return i != null ? getInt( i ) : dflt;
    }

    private static final Map<String, PrintStream> outputs = new HashMap<>();
    {
        outputs.put( "System.out", System.out );
        outputs.put( "System.err", System.err );
    }

    public static PrintStream getOutput( String name )
    {
        if( name == null )
            return null;

        PrintStream res = outputs.get( name );

        if( res == null )
        {
            try
            {
                res = new PrintStream( name );
            }
            catch( IOException e )
            {
                throw new RuntimeException( e );
            }

            outputs.put( name, res );
        }

        return res;
    }

    public void printAlgoMatrix( DistMatrix m )
    {
        printMatrix( m, getOutput( algo.get( "matrixOut" )));
    }

    public void printBenchmarkMatrix( DistMatrix m )
    {
        printMatrix( m, getOutput( benchmark.get( "matrixOut" )));
    }

    private void printMatrix( DistMatrix m, PrintStream out )
    {
        if( out != null && m != null )
            m.print( out );
    }
}