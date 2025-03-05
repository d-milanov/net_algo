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

public abstract class DataLoader
{
    public interface Filter 
    {
        boolean accept( String name, double a, double e, double i, double OM, double om );
    }

    public abstract ArrayList<Body> load( File f, Filter filter ) throws IOException;


    public static class Neodys extends DataLoader
    {
        public ArrayList<Body> load( File f, Filter filter ) throws IOException
        {
            return loadDataNeodys( f, filter );
        }
    }

    public static class Mpcorb extends DataLoader
    {
        public ArrayList<Body> load( File f, Filter filter ) throws IOException
        {
            return loadDataMpcorb( f, filter );
        }
    }

    public static class Cams extends DataLoader
    {
        public ArrayList<Body> load( File f, Filter filter ) throws IOException
        {
            return loadDataCams( f, filter );
        }
    }

    static ArrayList<Body> loadDataNeodys( File f, Filter filter ) throws IOException
    {
        return loadData( f, Pattern.compile( 
            "^'([^']+)'\\s+\\d+\\.\\d+\\s+"+
            "([\\d\\.E\\+-]+)\\s+" +
            "([\\d\\.E\\+-]+)\\s+" +
            "([\\d\\.E\\+-]+)\\s+" +
            "([\\d\\.E\\+-]+)\\s+" +
            "([\\d\\.E\\+-]+)" ),
            new int[]{1, 2, 3, 4, 5, 6},
            false,
            filter);
    }

    static ArrayList<Body> loadDataMpcorb( File f, Filter filter ) throws IOException
    {
        return loadData( f, Pattern.compile( 
            "^(\\S+)\\s+" +                  // Des'n
            "\\d+\\.\\d+\\s+"+               // H
            "\\d+\\.\\d+\\s+" +              // G
            "\\S+\\s+" +                     // Epoch
            "\\d+\\.\\d+\\s+"+               // M
            "([\\d\\.E\\+-]+)\\s+" +         // Peri
            "([\\d\\.E\\+-]+)\\s+" +         // Node
            "([\\d\\.E\\+-]+)\\s+" +         // Incl
            "([\\d\\.E\\+-]+)\\s+" +         // e
            "\\d+\\.\\d+\\s+"+               // n
            "([\\d\\.E\\+-]+)" ),            // a
            new int[]{1, 6, 5, 4, 3, 2},
            false,
            filter );
    }

    static ArrayList<Body> loadDataCams( File f, Filter filter ) throws IOException
    {
        return loadData( f, Pattern.compile( 
            "^([^,]+)," +      //#
            "(?:[^,]+,){58}" +
            "([^,]+)," +       //q
            "(?:[^,]+,){4}" +
            "([^,]+)," +       //e
            "(?:[^,]+,)" +
            "([^,]+)," +       //i
            "(?:[^,]+,)" +
            "([^,]+)," +       //om
            "(?:[^,]+,)" +
            "([^,]+)," ),      //Om
            new int[]{1, 2, 3, 4, 6, 5},
            true,
            filter);
    }


    static ArrayList<Body> loadData( File f, Pattern p, int[] name_a_e_i_Om_om, boolean isQ, Filter filter ) throws IOException
    {
        ArrayList<Body> res = new ArrayList<>();
        PrintStream dbgOut = System.out;

        BufferedReader r = null;

        try
        {
            r = new BufferedReader( new InputStreamReader( new FileInputStream( f )));

            String l;
            int bCnt = 0,
                vCnt = 0;

            while( (l = r.readLine()) != null )
            {
                Matcher m = p.matcher( l );

                if( !m.find())
                    continue;

                String name = m.group( name_a_e_i_Om_om[0] ),
                       aqs  = m.group( name_a_e_i_Om_om[1] ),
                       es   = m.group( name_a_e_i_Om_om[2] ),
                       is   = m.group( name_a_e_i_Om_om[3] ),
                       OMs  = m.group( name_a_e_i_Om_om[4] ),
                       oms  = m.group( name_a_e_i_Om_om[5] );

                double aq, e, i, OM, om;

                try
                {
                    aq = Double.parseDouble( aqs );
                    e = Double.parseDouble( es );
                    i = Double.parseDouble( is );
                    OM = Double.parseDouble( OMs );
                    om = Double.parseDouble( oms );
                }
                catch( NumberFormatException ex )
                {
                    dbgOut.println( "Cannot parse data: " + l );
                    continue;
                }

                if( filter == null || filter.accept( name, aq, e, i, OM, om ))
                    res.add( isQ ? Body.fromDegreesQ( bCnt++, name, aq, e, i, OM, om ) :
                             Body.fromDegrees( bCnt++, name, aq, e, i, OM, om ));
                vCnt++;
            }

            dbgOut.printf( "Loaded %d orbits (%.3f%% of all valid rows)%n", bCnt, vCnt != 0 ? 100.0*bCnt / (double)vCnt : 0.0 );
        }
        finally
        {
            if( r != null )
                r.close();
        }

        return res;
    }

    public static class  AFilter implements Filter 
    {
        final double a_min;
        final double a_max;
        final Reduce reduce;

        AFilter( double a_min, double a_max )
        {
            this( a_min, a_max, 1l );
        }

        AFilter( double a_min, double a_max, long times )
        {
            this.a_min = a_min;
            this.a_max = a_max;
            this.reduce = new Reduce( times );
        }

        public boolean accept( String name, double a, double e, double i, double OM, double om )
        {
            return a >= a_min && a < a_max &&
                   reduce.accept( name, a, e, i, OM, om );
        }

        public String toString()
        {
            return String.format( "%f <= a < %f, times = %d", a_min, a_max, reduce.times );
        }
    }

    public static class MainBelt extends AFilter
    {
        MainBelt()
        {
            this( 1l );
        }

        MainBelt( long times )
        {
            super( 1.523679, 5.2044, times );
        }
    }

    public static class  EFilter implements Filter 
    {
        final double e_min;
        final double e_max;

        EFilter( double e_min, double e_max )
        {
            this.e_min = e_min;
            this.e_max = e_max;
        }

        public boolean accept( String name, double a, double e, double i, double OM, double om )
        {
            return e >= e_min && e < e_max;
        }

        public String toString()
        {
            return String.format( "%f <= a < %f", e_min, e_max );
        }
    }

    public static class Reduce implements Filter 
    {
        final long times;
        long cnt;

        Reduce( long times )
        {
            this.times = times;
        }

        public boolean accept( String name, double a, double e, double i, double OM, double om )
        {
            return ((cnt++)%times) == 0;
        }

        public String toString()
        {
            return String.format( "times = %d", times );
        }
    }
}
