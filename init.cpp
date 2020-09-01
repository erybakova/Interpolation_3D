#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/sysinfo.h>

#include "init.h"

void get_k (int nx, int ny, int i, int j, int *k)
{
    (void) nx;
    *k = i * (ny + 1) + j;
}

void get_ij (int nx, int ny, int *i, int *j, int k)
{
    (void) nx;
    *j = k % (ny + 1);
    *i = (k - (*j)) / (ny + 1);
}

int get_akl_row (double *A, int *I, int nx, int ny, double hx, double hy, int k, bool h)
{
    int i, j, nz;

    get_ij (nx, ny, &i, &j, k);

    if ((i >= 1 && i <= nx - 1) && (j >= 1 && j <= ny - 1))
        nz = 6;
    else if ((i >= 1 && i <= nx - 1) || (j >= 1 && j <= ny - 1))
        nz = 4;
    else if ((i == 0 && j == 0) || (i == nx && j == ny))
        nz = 3;
    else nz = 2;

    if (!h) return nz;

    if ((i >= 1 && i <= nx - 1) && (j >= 1 && j <= ny - 1))
    {
        A[0] = hx * hy / 2.; 
        A[1] = hx * hy / 12.;
        A[2] = hx * hy / 12.;
        A[3] = hx * hy / 12.;
        A[4] = hx * hy / 12.;
        A[5] = hx * hy / 12.;
        A[6] = hx * hy / 12.;

        get_k (nx, ny, i - 1, j - 1, I + 0);
        get_k (nx, ny, i - 1, j,     I + 1);
        get_k (nx, ny, i,     j - 1, I + 2);
        get_k (nx, ny, i,     j + 1, I + 3);
        get_k (nx, ny, i + 1, j,     I + 4);
        get_k (nx, ny, i + 1, j + 1, I + 5);
    }
    else if (i >= 1 && i <= nx - 1 && j == 0)
    {
        A[0] = hx * hy / 4.; 
        A[1] = hx * hy / 24.;
        A[2] = hx * hy / 12.;
        A[3] = hx * hy / 24.;
        A[4] = hx * hy / 12.;

        get_k (nx, ny, i - 1, j    , I + 0);
        get_k (nx, ny, i,     j + 1, I + 1);
        get_k (nx, ny, i + 1,     j, I + 2);
        get_k (nx, ny, i + 1, j + 1, I + 3);
    }
    else if (i >= 1 && i <= nx - 1 && j == ny)
    {
        A[0] = hx * hy / 4.; 
        A[1] = hx * hy / 12.;
        A[2] = hx * hy / 24.;
        A[3] = hx * hy / 12.;
        A[4] = hx * hy / 24.;

        get_k (nx, ny, i - 1, j - 1, I + 0);
        get_k (nx, ny, i - 1,     j, I + 1);
        get_k (nx, ny, i,     j - 1, I + 2);
        get_k (nx, ny, i + 1,     j, I + 3);
    }
    else if (i == 0 && j >= 1 && j <= ny - 1)
    {
        A[0] = hx * hy / 4.; 
        A[1] = hx * hy / 24.;
        A[2] = hx * hy / 24.;
        A[3] = hx * hy / 12.;
        A[4] = hx * hy / 12.;

        get_k (nx, ny, i,     j - 1, I + 0);
        get_k (nx, ny, i,     j + 1, I + 1);
        get_k (nx, ny, i + 1,     j, I + 2);
        get_k (nx, ny, i + 1, j + 1, I + 3);
    }
    else if (i == nx && j >= 1 && j <= ny - 1)
    {
        A[0] = hx * hy / 4.; 
        A[1] = hx * hy / 12.;
        A[2] = hx * hy / 12.;
        A[3] = hx * hy / 24.;
        A[4] = hx * hy / 24.;

        get_k (nx, ny, i - 1, j - 1, I + 0);
        get_k (nx, ny, i - 1,     j, I + 1);
        get_k (nx, ny, i,     j - 1, I + 2);
        get_k (nx, ny, i,     j + 1, I + 3);
    }
    else if (i == 0 && j == 0)
    {
        A[0] = hx * hy / 6.; 
        A[1] = hx * hy / 24.;
        A[2] = hx * hy / 24.;
        A[3] = hx * hy / 12.;

        get_k (nx, ny, i,     j + 1, I + 0);
        get_k (nx, ny, i + 1,     j, I + 1);
        get_k (nx, ny, i + 1, j + 1, I + 2);
    }
    else if (i == nx && j == ny)
    {
        A[0] = hx * hy / 6.; 
        A[1] = hx * hy / 12.;
        A[2] = hx * hy / 24.;
        A[3] = hx * hy / 24.;

        get_k (nx, ny, i - 1, j - 1, I + 0);
        get_k (nx, ny, i - 1,     j, I + 1);
        get_k (nx, ny, i,     j - 1, I + 2);
    }
    else if (i == 0 && j == ny)
    {
        A[0] = hx * hy / 12.; 
        A[1] = hx * hy / 24.;
        A[2] = hx * hy / 24.;

        get_k (nx, ny, i, j - 1, I + 0);
        get_k (nx, ny, i + 1, j, I + 1);
    }
    else /* if (i == nx && j == 0) */
    {
        A[0] = hx * hy / 12.; 
        A[1] = hx * hy / 24.;
        A[2] = hx * hy / 24.;

        get_k (nx, ny, i - 1, j, I + 0);
        get_k (nx, ny, i, j + 1, I + 1);
    }

    return nz;
}

int get_non_zeros (double *A, int *I, int *row, int nx, int ny, double hx, double hy, int M)
{
    int k, nz, curr = 0;
    
    for (k = 0; k < M; k++)
    {
        row[k] = curr + M + 1;

        nz = get_akl_row (A, I, nx, ny, hx, hy, k, 0);

        curr += nz;
    }

    return curr;
}

void init_matrix (double *AA, int *IA, double *A, int *I, int *row,
                      int nx, int ny, double hx, double hy, int M, int max_nz, int q, int p)
{
    int k, nz = 0;
    int i1 = q * M; i1 /= p;
    int i2 = (q + 1) * M; i2 /= p;

    for (k = i1; k < i2; k++)
    {
        nz = get_akl_row (A + q * 7, I + q * 6, nx, ny, hx, hy, k, 1);

        AA[k] = (A + q * 7)[0];
        IA[k] = row[k];
        
        memcpy (AA + row[k], A + q * 7 + 1, nz * sizeof(double));
        memcpy (IA + row[k], I + q * 6, nz * sizeof(int));
    }
    if (q == 0)
      IA[M] = M + 1 + max_nz;
}

void init_vector (double (*f) (double, double), double *B, int nx, int ny, double hx, double hy, int M, int q, int p)
{
    int i, j, k, l;
    double total;
    
    double *a = new double[6 * 6];
    double *c = new double[6];
    double *x = new double[6];
    double *y = new double[6];
    int *s = new int[6];

    init_s (s, 6);

    int i1 = q * M; i1 /= p;
    int i2 = (q + 1) * M; i2 /= p;

    for (k = i1; k < i2; k++)
    {
        total = 0.;

        get_ij (nx, ny, &i, &j, k);

        if ((i >= 1 && i <= nx - 1) && (j >= 1 && j <= ny - 1))
        {
            /* 6 triangles: 1, 2, 3, 4, 5, 6 */
            for (l = 1; l <= 6; l++)
              total += int_over_triangle_l (l, a, c, x, y, s, f, hx, hy, i, j);
        }
        else if (i >= 1 && i <= nx - 1 && j == 0)
        {
            /* 3 triangles: 1, 5, 6 */
            for (l = 1; l <= 6; l++)
              if (l == 1 || l == 5 || l == 6) 
                total += int_over_triangle_l (l, a, c, x, y, s, f, hx, hy, i, j);
        }
        else if (i >= 1 && i <= nx - 1 && j == ny)
        {
            /* 3 triangles: 2, 3, 4 */
            for (l = 1; l <= 6; l++)
              if (l == 2 || l == 3 || l == 4) 
                total += int_over_triangle_l (l, a, c, x, y, s, f, hx, hy, i, j);
        }
        else if (i == 0 && j >= 1 && j <= ny - 1)
        {
            /* 3 triangles: 1, 2, 6 */
            for (l = 1; l <= 6; l++)
              if (l == 1 || l == 2 || l == 6) 
                total += int_over_triangle_l (l, a, c, x, y, s, f, hx, hy, i, j);
        }
        else if (i == nx && j >= 1 && j <= ny - 1)
        {
            /* 3 triangles: 3, 4, 5 */
            for (l = 1; l <= 6; l++)
              if (l == 3 || l == 4 || l == 5) 
                total += int_over_triangle_l (l, a, c, x, y, s, f, hx, hy, i, j);
        }
        else if (i == 0 && j == 0)
        {
            /* 2 triangles: 1, 6 */
            for (l = 1; l <= 6; l++)
              if (l == 1 || l == 6) 
                total += int_over_triangle_l (l, a, c, x, y, s, f, hx, hy, i, j);
        }
        else if (i == nx && j == ny)
        {
            /* 2 triangles: 3, 4 */
            for (l = 1; l <= 6; l++)
              if (l == 3 || l == 4) 
                total += int_over_triangle_l (l, a, c, x, y, s, f, hx, hy, i, j);
        }
        else if (i == 0 && j == ny)
        {
            /* 1 triangle: 2 */
            total += int_over_triangle_l (2, a, c, x, y, s, f, hx, hy, i, j);
        }
        else /* if (i == nx && j == 0) */
        {
            /* 1 triangle: 5 */
            total += int_over_triangle_l (5, a, c, x, y, s, f, hx, hy, i, j);
        }

        B[k] = total;
    }

    delete [] a;
    delete [] c;
    delete [] x;
    delete [] y;
    delete [] s;
}

void init_vector_6 (double *x, double s, double x1, double x2, double x3, 
                                 double x4, double x5, double x6)
{
    x[0] = s + x1;
    x[1] = s + x2;
    x[2] = s + x3;
    x[3] = s + x4;
    x[4] = s + x5;
    x[5] = s + x6;
}

void init_vector_c (double *c, double (*f) (double, double), double *x, double *y)
{
    int i;

    for (i = 0; i < 6; i++)
      c[i] = f (x[i], y[i]);
}

void init_matrix_6 (double *a, double *x, double *y)
{
    int i;
  
    for (i = 0; i < 6; i++)
    {
        (a + i * 6)[0] = x[i] * x[i];
        (a + i * 6)[1] = x[i] * y[i];
        (a + i * 6)[2] = y[i] * y[i];
        (a + i * 6)[3] = x[i];
        (a + i * 6)[4] = y[i];
        (a + i * 6)[5] = 1.;
    }
}

double int_over_triangle_l (int l, double *a, double *c, double *x, double *y, int *s,
                               double (*f) (double, double), double hx, double hy, int i, int j)
{
   /* double m1, m2, m3, m4, m5, m6;

    m1 = hx * hx * hx * hy;
    m2 = hx * hx * hy * hy;
    m3 = hx * hy * hy * hy;
    m4 = hx * hx * hy;
    m5 = hx * hy * hy;
    m6 = hx * hy;*/

    switch (l)
    {
        case 1:
            init_vector_6 (x, i * hx, 0, hx / 2., hx, hx, hx, hx / 2.);
            init_vector_6 (y, j * hy, 0, hy / 2., hy, hy / 2., 0, 0);
            init_matrix_6 (a, x, y);
            init_vector_c (c, f, x, y);

            solve (a, c, s, 6);

            /*return c[0] * m1 / 20. + c[1] * m2 / 40. + c[2] * m3 / 60. 
                   + c[3] * m4 / 12. + c[4] * m5 / 24. + c[5] * m6 / 6.;*/
            return (1. / 120) * hx * hy * (20 * c[5] + 10 * c[3] * (hx + 2 * hx * i) 
                       + 2 * c[0] * hx * hx * (3 + 10 * i * (1 + i)) 
                           + hy * (5 * c[4] * (1 + 4 * j) + 2 * c[2] * hy * (1 + 5 * j * (1 + 2 * j)) 
                               + c[1] * hx * (3 + 10 * j + 5 * i * (1 + 4 * j))));
        case 2:
            init_vector_6 (x, i * hx, 0, hx / 2., hx, hx / 2., 0, 0);
            init_vector_6 (y, j * hy, 0, 0, 0, -hy / 2., -hy, -hy / 2.);
            init_matrix_6 (a, x, y);
            init_vector_c (c, f, x, y);

            solve (a, c, s, 6);

            /*return c[0] * m1 / 60. - c[1] * m2 / 120. + c[2] * m3 / 60. 
                   + c[3] * m4 / 24. - c[4] * m5 / 24. + c[5] * m6 / 6.;*/
            return (1. / 120) * hx * hy * (20 * c[5] - 5 * c[4] * hy + 5 * c[3] * (hx + 4 * hx * i) 
                       + 2 * c[0] * hx * hx * (1 + 5 * i * (1 + 2 * i)) 
                           + hy * (2 * c[2] * hy + 20 * c[4] * j + 10 * c[2] * hy * j * (-1 + 2 * j) 
                               + c[1] * hx * (-1 + 5 * j + 5 * i * (-1 + 4 * j))));
        case 3:
            init_vector_6 (x, i * hx, 0, 0, 0, -hx / 2., -hx, -hx / 2.);
            init_vector_6 (y, j * hy, 0, -hy / 2., -hy, -hy, -hy, -hy / 2.);
            init_matrix_6 (a, x, y);
            init_vector_c (c, f, x, y);

            solve (a, c, s, 6);

            /*return c[0] * m1 / 60. + c[1] * m2 / 40. + c[2] * m3 / 20. 
                   - c[3] * m4 / 24. - c[4] * m5 / 12. + c[5] * m6 / 6.;*/
            return (1. / 120) * hx * hy * (20 * c[5] - 10 * c[4] * hy + 5 * c[3] * hx * (-1 + 4 * i) 
                       + 2 * c[0] * hx * hx * (1 + 5 * i * (-1 + 2 * i)) 
                           + hy * (20 * c[4] * j + 2 * c[2] * hy * (3 + 10 * (-1 + j) * j) 
                               + c[1] * hx * (3 - 5 * j + 10 * i * (-1 + 2 * j))));
        case 4:
            init_vector_6 (x, i * hx, 0, -hx / 2., -hx, -hx, -hx, -hx / 2.);
            init_vector_6 (y, j * hy, 0, -hy / 2., -hy, -hy / 2., 0, 0);
            init_matrix_6 (a, x, y);
            init_vector_c (c, f, x, y);

            solve (a, c, s, 6);

            /*return c[0] * m1 / 20. + c[1] * m2 / 40. + c[2] * m3 / 60. 
                   - c[3] * m4 / 12. - c[4] * m5 / 24. + c[5] * m6 / 6.;*/
            return (1. / 120) * hx * hy * (20 * c[5] - 5 * c[4] * hy + 10 * c[3] * hx * (-1 + 2 * i) 
                       + 2 * c[0] * hx * hx * (3 + 10 * (-1 + i) * i) 
                           + hy * (2 * c[2] * hy + 20 * c[4] * j + 10 * c[2] * hy * j * (-1 + 2 * j) 
                               + c[1] * hx * (3 - 10 * j + 5 * i * (-1 + 4 * j))));
        case 5:
            init_vector_6 (x, i * hx, 0, -hx / 2., -hx, -hx / 2., 0, 0);
            init_vector_6 (y, j * hy, 0, 0, 0, hy / 2., hy, hy / 2.);
            init_matrix_6 (a, x, y);
            init_vector_c (c, f, x, y);

            solve (a, c, s, 6);

            /*return c[0] * m1 / 60. - c[1] * m2 / 120. + c[2] * m3 / 60. 
                   - c[3] * m4 / 24. + c[4] * m5 / 24. + c[5] * m6 / 6.;*/
            return (1. / 120) * hx * hy * (20 * c[5] + 5 * c[3] * hx * (-1 + 4 * i) 
                       + 2 * c[0] * hx * hx * (1 + 5 * i * (-1 + 2 * i))
                           + hy * (5 * c[4] * (1 + 4 * j) + c[1] * hx * (-1 + 5 * i - 5 * j + 20 * i * j) 
                               + 2 * c[2] * hy * (1 + 5 * j * (1 + 2 * j))));
        case 6:
            init_vector_6 (x, i * hx, 0, 0, 0, hx / 2., hx, hx / 2.);
            init_vector_6 (y, j * hy, 0, hy / 2., hy, hy, hy, hy / 2.);
            init_matrix_6 (a, x, y);
            init_vector_c (c, f, x, y);

            solve (a, c, s, 6);

            /*return c[0] * m1 / 60. + c[1] * m2 / 40. + c[2] * m3 / 20. 
                   + c[3] * m4 / 24. + c[4] * m5 / 12. + c[5] * m6 / 6.;*/
            return (1. / 120) * hx * hy * (20 * c[5] + 10 * c[4] * hy + 5 * c[3] * (hx + 4 * hx * i) 
                       + 2 * c[0] * hx * hx * (1 + 5 * i * (1 + 2 * i)) 
                           + hy * (20 * c[4] * j + 2 * c[2] * hy * (3 + 10 * j * (1 + j)) 
                               + c[1] * hx * (3 + 5 * j + 10 * i * (1 + 2 * j))));
        default:
            printf ("Unknown error!\n");

            return -1;
    }
}

double division (double *a, double *b, int i, int n, double z)
{
    int j = 0;
    double f = a[i * n + i];

    if (fabs(f) < 1e-14 * z)
      return -1;
    else
    {
        a[i * n + i] = 1.;

        for (j = i + 1; j < n; j++)
          a[i * n + j] /= f;

        b[i] /= f;
    }

    return 0;
}

void subtraction (double *a, double *e, int i, int n)
{
    int str, j;
    double f;

    for (str = i + 1; str < n; str++)
    {
        f = a[str * n + i];
        a[str * n + i] = 0.;

        for (j = i + 1; j < n; j++)
            a[str * n + j] -= f * a[i * n + j];

        e[str] -= f * e[i];
    }
}

void answer (double *a, double *b, int n)
{
    int i, j;

    for (i = n - 2; i >= 0; i--)
        for (j = n - 1; j > i; j--)
            b[i] -= a[i * n + j] * b[j];
}

void swap (double *a, double *b, int *s, int n, int i)
{
    int j;

    for (j = i + 1; j < n; j++)
    {
        /* swap i row with j row */
        if (fabs (a[j * n + i]) > 1e-14)
        {
            swap_rows (a, b, n, i, j);
            s[i] = j;
            break;
        }
    }
}

void swap_rows (double *a, double *b, int n, int i, int j)
{
    int k;
    double curr;

    for (k = 0; k < n; k++)
    {
        curr = a[i * n + k];
        a[i * n + k] = a[j * n + k];
        a[j * n + k] = curr;
    }

    curr = b[i];
    b[i] = b[j];
    b[j] = curr;
}

void swap_back (double *a, double *b, int *s, int n)
{
    int k;

    for (k = n - 1; k >= 0; k--)
    {
        if (s[k] >= 0)
          swap_rows (a, b, n, k, s[k]);
    }
}

int solve (double *a, double *b, int *s, int n)
{
    int i, res;
    double z = norm_matrix (a, n);

    /* forward */
    for (i = 0; i < n; i++)
    {
        res = division (a, b, i, n, z);
        if (res < 0)
        {
            swap (a, b, s, n, i);
            division (a, b, i, n, z);
        }
        if (i != n - 1) subtraction (a, b, i, n);
    }
    
    /* reverse */    
    answer (a, b, n);

    //swap_back (a, b, s, n);

    return 0;
}

void init_s (int *s, int n)
{
    int i;
    
    for (i = 0; i < n; i++)
      s[i] = -1;
}

double norm_matrix (double *a, int n)
{
    int i, j;
    double s, max = 0.;

    for (j = 0; j < n; j++)
    {
        s = 0.;

        for (i = 0; i < n; i++)
            s += fabs(a[i * n + j]);

        if (s > max) max = s;
    }
    return max;
}

double norm_msr_matrix (double *AA, int *IA, int M)
{
    int i, j;
    int l, m;
    double max = AA[0];

    for (i = 0; i < M; i++)
    {
        if (AA[i] > max)
            max = AA[i]; 

        l = IA[i + 1] - IA[i];
        m = IA[i];

        for (j = 0; j < l; j++)
          if (AA[m + j] > max)
            max = AA[m + j];     
    }

    return max;
}

void print_msr_matrix (double *AA, int *IA, int n, int M)
{
    int i;

    printf ("AA: ");
    for (i = 0; i < n; i++)
    {
        if (i == M) printf ("| * | ");
        else printf ("%.3g ", AA[i]);
    }
    printf ("\n\n");

    printf ("IA: ");
    for (i = 0; i < n; i++)
      printf ("%d ", IA[i]);

    printf ("\n");
}

void print_matrix (double *AA, int *IA, int M)
{
    int i, j, m, curr;

    for (i = 0; i < M; i++)
    {
        m = IA[i];

        for (j = 0, curr = 0; curr < M; curr++)
        {
            if (curr == i)
            {
                printf ("%6.3f ", AA[i]);
            }
            else if (curr == IA[j + m])
            { 
                printf ("%6.3f ", AA[m + j]);
                j++;
            }
            else
            {
                printf ("%6d ", 0);
            }
        }
        printf ("\n");
    }
    printf ("\n");
}

void print_vector (double *X, int M)
{
    int i;

    for (i = 0; i < M; i++)
      printf ("%.3f ", X[i]);

    printf ("\n");
}

void print (double *a, int n)
{
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf ("%4.3g ", a[i * n + j]);
        }
        printf ("\n");
    }
    printf ("\n");
}
