#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/sysinfo.h>

#include "init.h"
#include "solve.h"

#define MAX_IT 50

void init_vector_7 (double *x, double x1, double x2, double x3, 
                                 double x4, double x5, double x6, double x7)
{
    x[0] = x1;
    x[1] = x2;
    x[2] = x3;
    x[3] = x4;
    x[4] = x5;
    x[5] = x6;
    x[6] = x7;
}

double go_through_triangle (double *x, double *y, double (*f) (double, double), 
                                double *c, int nx, int ny, int i, int j, bool h)
{
    int l, k1, k2, k3;
    double curr, max = 0.;

    get_k (nx, ny, i, j, &k1);
    
    if (!h)
    {
        get_k (nx, ny, i, j + 1, &k2);
        get_k (nx, ny, i + 1, j + 1, &k3);
    }
    else
    {
        get_k (nx, ny, i + 1, j + 1, &k2);
        get_k (nx, ny, i + 1, j, &k3);
    }

    for (l = 0; l < 7; l++)
    {
        switch (l)
        {
            case 0: 
                curr = fabs (f (x[0], y[0]) - c[k1]);
                break;
            case 1:
                curr = fabs (f (x[1], y[1]) - c[k2]);
                break;
            case 2:
                curr = fabs (f (x[2], y[2]) - c[k3]);
                break; 
            case 3: 
                curr = fabs (f (x[3], y[3]) - (c[k1] + c[k2]) / 2.);
                break;
            case 4:
                curr = fabs (f (x[4], y[4]) - (c[k2] + c[k3]) / 2.);
                break;
            case 5:
                curr = fabs (f (x[5], y[5]) - (c[k1] + c[k3]) / 2.);
                break;
            case 6:
                curr = fabs (f (x[6], y[6]) - (c[k1] + c[k2] + c[k3]) / 3.);
                break;   
        }

        if (curr > max)
          max = curr;
    }
    
    return max;
}

void calc_residual (double (*f) (double, double), double *c, double *r, int nx, int ny, double hx, double hy, int q, int p)
{
    int i, j, k, N;
    double res1, res2;

    N = nx * ny;
    
    double *x = new double[7];
    double *y = new double[7];

    int i1 = q * N; i1 /= p;
    int i2 = (q + 1) * N; i2 /= p;

    for (k = i1; k < i2; k++)
    {
        i = k / ny;
        j = k - i * ny;

        init_vector_7 (x, i * hx, i * hx, i * hx + hx, i * hx, i * hx + hx / 2., i * hx + hx / 2., i * hx + hx / 3.);
        init_vector_7 (y, j * hy, j * hy + hy, j * hy + hy, j * hy + hy / 2., j * hy + hy, j * hy + hy / 2., j * hy + 2 * hy / 3.);

        res1 = go_through_triangle (x, y, f, c, nx, ny, i, j, 0);

        init_vector_7 (x, i * hx, i * hx + hx, i * hx + hx, i * hx + hx / 2., i * hx + hx, i * hx + hx / 2., i * hx + 2 * hx / 3.);
        init_vector_7 (y, j * hy, j * hy + hy, j * hy, j * hy + hy / 2., j * hy + hy / 2., j * hy, j * hy + hy / 3.);

        res2 = go_through_triangle (x, y, f, c, nx, ny, i, j, 1);

        r[k] = (res1 > res2 ? res1 : res2);
    }

    delete [] x;
    delete [] y;
}

double choose_max (double *r, int N)
{
    int i;
    double max = r[0];

    for (i = 1; i < N; i++)
      if (r[i] > max)
        max = r[i];
    
    return max;
}

void msr_matrix_mult_vector (double *AA, int *IA, int M, double *x, double *b, int q, int p)
{
    int i, j;
    int l, m;
    double s;

    int i1 = q * M; i1 /= p;
    int i2 = (q + 1) * M; i2 /= p;

    for (i = i1; i < i2; i++)
    {
        s = AA[i] * x[i];

        l = IA[i + 1] - IA[i];
        m = IA[i];

        for (j = 0; j < l; j++)
          s += AA[m + j] * x[IA[m + j]];

        b[i] = s;     
    }
}

double scalar_product (double *x, double *y, int M, int q, int p, double *buf, pthread_barrier_t *barrier)
{
    int i, j;
    double s = 0;

    int i1 = q * M; i1 /= p;
    int i2 = (q + 1) * M; i2 /= p;

    for (j = i1; j < i2; j++)
      s += x[j] * y[j];

    buf[q] = s;
    s = 0;

    pthread_barrier_wait (barrier);

    for (i = 0; i < p; i++)
    {
        s += buf[i];
    }

    return s;
}

void linear_comb (double *x, double *y, double tau, int M, int q, int p)
{
    int i;

    int i1 = q * M; i1 /= p;
    int i2 = (q + 1) * M; i2 /= p;

    for (i = i1; i < i2; i++)
      x[i] -= tau * y[i];
}

void apply_preconditioner (double *AA, double *v, double *r, int M, int q, int p)
{
    /* M is diag of AA */
    int i;

    int i1 = q * M; i1 /= p;
    int i2 = (q + 1) * M; i2 /= p;

    for (i = i1; i < i2; i++)
      v[i] = r[i] / AA[i];
}

double norm_matrix_M (double *A, int n)
{
    /* M is diag of A */
    int j;
    double max = A[0];

    for (j = 1; j < n; j++)
    {
        if (A[j] > max) max = A[j];
    }

    return max;
}

void init_0 (double *x, int M, int q, int p)
{
    int i;

    int i1 = q * M; i1 /= p;
    int i2 = (q + 1) * M; i2 /= p;

    for (i = i1; i < i2; i++)
      x[i] = 0.;
}

int solve_stage (double *AA, int *IA, double *X, double *B, double *r, double *u, double *v,
                     int M, int max_it, int q, int p, double *buf, pthread_barrier_t *barrier, double norm, double eps)
{
    int it;
    double c1, c2, tau;

    /* r = Ax; */
    msr_matrix_mult_vector (AA, IA, M, X, r, q, p);
    pthread_barrier_wait (barrier);
    /* r -= b; */
    linear_comb (r, B, 1., M, q, p);
    pthread_barrier_wait (barrier);

    for (it = 0; it < max_it; it++)
    {
        /* Mv = r; */
        apply_preconditioner (AA, v, r, M, q, p);
        pthread_barrier_wait (barrier);

        /* u = Av; */
        msr_matrix_mult_vector (AA, IA, M, v, u, q, p);
        pthread_barrier_wait (barrier);

        /* c1 = (u, r); */
        c1 = scalar_product (u, r, M, q, p, buf, barrier);
        pthread_barrier_wait (barrier);

        /* c2 = (u, u); */
        c2 = scalar_product (u, u, M, q, p, buf, barrier);
        pthread_barrier_wait (barrier);

        if (c1 <= eps * eps * norm || c2 <= eps * eps * norm * norm)
          break;

        tau = c1 / c2;

        /* x = x - tau * v; */
        linear_comb (X, v, tau, M, q, p);

        /* r = r - tau * u; */
        linear_comb (r, u, tau, M, q, p);
    }

    if (it >= max_it)
      return -1;
    else
      return it;
}

double norm_matrix_msr (double *A, int M)
{
    int i;
    double max = A[0];

    for (i = 1; i < M; i++)
    {
        if (A[i] > max)
          max = A[i];
    }

    return max;
}

int solve (double *AA, int *IA, double *X, double *B, double *r, double *u, double *v, int M, int max_it,
               int q, int p, double *buf, pthread_barrier_t *barrier, double eps)
{
    int it, res, IT = 0;
    double norm;

    /* X = 0.; */
    init_0 (X, M, q, p);
    norm = norm_matrix_msr (AA, M);
    pthread_barrier_wait (barrier);

    for (it = 0; it <= max_it; it += MAX_IT)
    {
        res = solve_stage (AA, IA, X, B, r, u, v, M, MAX_IT, q, p, buf, barrier, norm, eps);
        pthread_barrier_wait (barrier);

        IT += res;

        if (res >= 0)
          break;
    }

    if (it >= max_it)
      return -1;
    else
      return IT;
}
