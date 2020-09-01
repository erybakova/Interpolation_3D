void get_k (int nx, int ny, int i, int j, int *k);
void get_ij (int nx, int ny, int *i, int *j, int k);
int get_akl_row (double *A, int *I, int nx, int ny, double hx, double hy, int k, bool h);
int get_non_zeros (double *A, int *I, int *row, int nx, int ny, double hx, double hy, int M);
void init_matrix (double *AA, int *IA, double *A, int *I, int *row, int nx, int ny,
                      double hx, double hy, int M, int max_nz, int q, int p);
void init_vector (double (*f) (double, double), double *X, int nx, int ny, double hx, double hy, int M, int q, int p);
void init_vector_6 (double *x, double s, double x1, double x2, double x3, 
                                 double x4, double x5, double x6);
void init_vector_c (double *c, double (*f) (double, double), double *x, double *y);
void init_matrix_6 (double *a, double *x, double *y);
double int_over_triangle_l (int l, double *a, double *c, double *x, double *y, int *s, 
                               double (*f) (double, double), double hx, double hy, int i, int j);
double division (double *a, double *e, int i, int n, double z);
void subtraction (double *a, double *e, int i, int n);
void answer (double *a, double *b, int n);
void swap (double *a, double *b, int *s, int n, int i);
void swap_rows (double *a, double *b, int n, int i, int j);
void swap_back (double *a, double *b, int *s, int n);
void init_s (int *s, int n);
int solve (double *a, double *b, int *s, int n);
void init_s (int *s, int n);
double norm_matrix (double *a, int n);
double norm_msr_matrix (double *AA, int *IA, int M);
void print_msr_matrix (double *AA, int *IA, int n, int M);
void print_matrix (double *AA, int *IA, int M);
void print_vector (double *X, int M);
void print (double *a, int n);
