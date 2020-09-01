void *thread_function (void *args);
double approximation (double (*f)(double, double), double hx, double hy,
                        int nx, int ny, int M, int p, double *X, double eps);
