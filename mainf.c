#include <stdio.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include "timer.c"
#include "timer.h"

extern double g (double *t, size_t dim, void *params);

double dipole_approx (double r);

double gaussian (double *x, int dim);

int main (void)
{
    double res, err;

    size_t dim = 6;
    double x1[] = { 0., 0., 0., 0., 0., 0., };
    double xu[] = { 1., 1., 1., 1., 1., 1., };
    double distmin = 1.001;
    double distmax = 4.;
    double dist;
    int np = 20;
    double nt = (distmax - distmin) / (np - 1);
    double vegas[20], dipole[20], distance[20];

    gsl_rng *r = gsl_rng_alloc (gsl_rng_taus2);
    unsigned long seed = 1UL;

    gsl_rng_set (r, seed);

    size_t calls = 1000000;
   
    dist = distmin;

    gsl_monte_function G = { &g, dim, &dist };

    gsl_monte_vegas_state *sv = gsl_monte_vegas_alloc (dim);

    gsl_monte_vegas_init (sv);

    // Vegas Integration

    /*commented so does not output in order to create a proper res
    printf ("#  Stat      Dist       E(r)        ErrEst      Dipolapprox\n");
   */

    timer_start ();

    for (int i = 0; i < np; i++)
    {
        gsl_monte_vegas_integrate (&G, x1, xu, dim, calls / 5, r, sv, &res,
            &err);
        do
        {
            gsl_monte_vegas_integrate (&G, x1, xu, dim, calls, r, sv, &res,
                &err);
        }
        while (fabs (gsl_monte_vegas_chisq (sv) - 1.0) > 0.2);
        /*commented so does not output in order to create a proper res
        printf ("   %.2f    %.6f   %.6f    %.6f      %.6f\n",
            gsl_monte_vegas_chisq (sv), dist, res, err, -2. / pow (dist, 3.));
        */

        fflush (stdout);
        dist += nt;
        vegas[i] = res;
    distance[i] = dist;
    dipole[i] = -2. / pow (dist, 3.);
    }

    timer_stop();
    /*commented so does not output in order to create a proper res
    double t1 = timer_stop ();
    */

    gsl_monte_vegas_free (sv);


    // Monte Carlo Integration

    double sum;
    double x[6];

    long i, j, nn;

    /*commented so does not output in order to create a proper res
    printf ("# Iter        Npoints      AbsErr\n");
    */

    nn = 1000000;

    timer_start ();
    double home[20];
  
    dist = distmin;
    for (j = 0; j < np; j++)
    {
        sum = 0.;
        for (i = 0; i < nn; i++)
        {
            // 6-dimensional random point
            for (int k = 0; k < (int) dim; k++)
            {
                x[k] = gsl_rng_uniform (r);
            }
            sum += g (x, dim, &dist);
        }
        res = sum/nn;
        /*commented so does not output in order to create a proper res
        printf ("   %.6f   %.6f   %.6f\n",
             dist, res, -2. / pow (dist, 3.));
        */
        dist += nt;
        home[j] = res;
    }   

    timer_stop ();
    /*commented so does not output in order to create a proper res
    double t2 = timer_stop();
    */

    /*commented so does not output in order to create a proper res
    printf ("Time for Vegas:   %6f\n", t1);
    printf ("Time for HomeMC:  %6f\n", t2);
    */

    gsl_rng_free (r);

    //calculating error
    double homeerr = 0.0;

    for (int e = 0; e < np; e++)
    {
    homeerr += fabs(home[e] - vegas[e]);
    }

    /*commented so does not output in order to create a proper res
    double homeerr_avg = homeerr/np;
    double differ = homeerr_avg/err;
    */

    printf("#    Dist               Vegas           Home       Dipolapprox\n");
    for( int l = 0; l < np; l++)
    {
    double dd = distance[l];
    double vv = fabs(vegas[l]);
    double hh = fabs(home[l]);
    double di = fabs(dipole[l]);
    printf("   %.6f         %.6f       %.6f      %.6f\n", dd, vv, hh, di);
    }

    /*commented so does not output in order to create a proper res
    printf("\nAverage error for Home integration: %.6f\n", homeerr_avg);
    printf("Average error for Vegas integration: %.6f\n", err);
    printf("Ratio of the two errors: %.6f\n", differ);
    */


    return 0;
}
