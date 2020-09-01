#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/sysinfo.h>

#include "init.h"
#include "solve.h"

static double max = 0;

struct arg
{
    double *AA, *A;
    int *IA, *I;
    double *X, *B;
    double *r, *u, *v;
    int *row;
    double *resid;
    int nx, ny;
    double hx, hy;
    int max_nz, M;
    int q, p;
    double eps;
    double *buf;
    double (*f) (double, double);
    pthread_barrier_t *barrier;
};

void *thread_function (void *args)
{
    arg *d = (arg*) args;
    double *AA = d -> AA;
    int *IA = d -> IA;
    double *A = d -> A;
    int *I = d -> I;
    double *X = d -> X, *B = d -> B;
    double *r = d -> r, *u = d -> u, *v = d -> v;
    int *row = d -> row;
    double *resid = d -> resid;
    int nx = d -> nx, ny = d -> ny;
    double hx = d -> hx, hy = d -> hy;
    int max_nz = d -> max_nz;
    int M = d -> M;
    int q = d -> q;
    int p = d -> p;
    double eps = d -> eps;
    double *buf = d -> buf;
    double (*f) (double, double) = d -> f;
    int res;
    double residual;

    init_matrix (AA, IA, A, I, row, nx, ny, hx, hy, M, max_nz, q, p);
    
    pthread_barrier_wait (d -> barrier);
    
    init_vector (f, B, nx, ny, hx, hy, M, q, p);
    
    pthread_barrier_wait (d -> barrier);

    if (q == 0) 
    {
        /*printf ("A:\n");
        print_matrix (AA, IA, M);
        print_msr_matrix (AA, IA, M + max_nz + 1, M);
        printf ("\nB: ");
        print_vector (B, M);*/
    }

    res = solve (AA, IA, X, B, r, u, v, M, 10000, q, p, buf, d -> barrier, eps);

    if (q == 0)
    {
        /*printf ("\nX: ");
        print_vector (X, M);*/
        printf ("\nIterations: %d\n", res);
    }

    pthread_barrier_wait (d -> barrier);

    calc_residual (f, X, resid, nx, ny, hx, hy, q, p);

    pthread_barrier_wait (d -> barrier);
    
    if (q == 0)
    {
        residual = choose_max (resid, nx * ny);
        printf ("Residual: %.5e\n", residual);
        max = residual;
    }

    return 0;
}

double approximation (double (*f)(double, double), double hx, double hy,
                        int nx, int ny, int M, int p, double *X, double eps)
{
    int N, max_nz;
    double *AA, *A = 0, *B;
    double *r, *u, *v, *resid, *buf;
    int *IA, *I = 0, *row;
    pthread_t * threads;
    pthread_barrier_t barrier;
    int i;
    
    N = nx * ny;

    row = new int[M];
    resid = new double[N];

    max_nz = get_non_zeros (A, I, row, nx, ny, hx, hy, M);
    
    AA = new double[M + max_nz + 1];
    IA = new int[M + max_nz + 1];
    B = new double[M];
    A = new double[7 * p];
    I = new int[6 * p];
    buf = new double[p];
    r = new double[M];
    u = new double[M];
    v = new double[M];

    arg *d = new arg[p];
    threads = new pthread_t [p * sizeof(pthread_t)];

    if (!threads)
    {
        fprintf (stderr, "Not enough memory!\n");
        return 0;
    }

    pthread_barrier_init (&barrier, 0, p);

    for (i = 0; i < p; i++)
    {
        d[i].AA = AA;
        d[i].IA = IA;
        d[i].X = X;
        d[i].B = B;
        d[i].A = A;
        d[i].I = I;
        d[i].r = r;
        d[i].u = u;
        d[i].v = v;
        d[i].row = row;
        d[i].resid = resid;
        d[i].nx = nx;
        d[i].ny = ny;
        d[i].hx = hx;
        d[i].hy = hy;
        d[i].max_nz = max_nz;
        d[i].M = M;
        d[i].q = i;
        d[i].p = p;
        d[i].buf = buf;
        d[i].f = f;
        d[i].eps = eps;
        d[i].barrier = &barrier;
    }
 
    for (i = 1; i < p; i++)
    {
        if (pthread_create (threads + i, 0, thread_function, d + i))
        {
            fprintf (stderr, "Cannot create thread %d!\n", i);
            return 0;
        }
    }

    thread_function (d + 0);

    pthread_barrier_destroy (&barrier);

    delete [] threads;
    delete [] d;
    delete [] AA;
    delete [] IA;
    delete [] B;
    delete [] A;
    delete [] I;
    delete [] row;
    delete [] resid;
    delete [] buf;
    delete [] r;
    delete [] u;
    delete [] v;

    return max;
}
