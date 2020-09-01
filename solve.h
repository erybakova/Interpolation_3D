void init_vector_7 (double *x, double x1, double x2, double x3, 
                                 double x4, double x5, double x6, double x7);
double go_through_triangle (double *x, double *y, double (*f) (double, double), 
                                double *c, int nx, int ny, int i, int j, bool h);
void calc_residual (double (*f) (double, double), double *c, double *r, int nx, int ny, double hx, double hy, int q, int p);
double choose_max (double *r, int M);
void msr_matrix_mult_vector (double *AA, int *IA, int M, double *x, double *b, int q, int p);
double scalar_product (double *x, double *y, int M, int q, int p, double *buf, pthread_barrier_t *barrier);
void linear_comb (double *x, double *y, double tau, int M, int q, int p);
void apply_preconditioner (double *AA, double *v, double *r, int M, int q, int p);
double norm_matrix_M (double *A, int n);
void init_0 (double *x, int M, int q, int p);
int solve_stage (double *AA, int *IA, double *X, double *B, double *r, double *u, double *v,
                 int M, int max_it, int q, int p, double *buf, pthread_barrier_t *barrier, double norm, double eps);
int solve (double *AA, int *IA, double *X, double *B, double *r, double *u, double *v, int M, int max_it,
               int q, int p, double *buf, pthread_barrier_t *barrier, double eps);
double norm_matrix_msr (double *A, int M);
