#include "glwidget.h"
#include "approximation.h"

#include <stdio.h>
#include <algorithm>
#include <math.h>

#define LEN 55

double f_0 (double, double)
{
    return 1;
}
double f_1 (double x, double)
{
    return x;
}
double f_2 (double, double y)
{
    return y;
}
double f_3 (double x, double y)
{
    return x + y;
}
double f_4 (double x, double y)
{
    return sqrt(x * x + y * y);
}
double f_5 (double x, double y)
{
    return x * x + y * y;
}
double f_6 (double x, double y)
{
    return exp (x * x - y * y);
}
double f_7 (double x, double y)
{
    return 1. / (25 * (x * x + y * y) + 1);
}

myGLWidget:: myGLWidget (QWidget *parent/* = 0 */) : QGLWidget(parent)
{
    xRot = -60.;
    yRot = 0.;
    zRot = 25.;
    nSca = 1.;
}

void myGLWidget:: init (double a_, double b_, double alpha_,
                            int nx_, int ny_, int k_, double eps_, int p_)
{
    a = a_;
    b = b_;
    alpha = alpha_;
    p = p_;
    nx = nx_;
    ny = ny_;
    nx1 = nx_;
    ny1 = ny_;
    hx = a / nx;
    hy = b / ny;
    M = (nx + 1) * (ny + 1);
    k = k_;
    eps = eps_;
    s = 0;
    h = 0;
    max_abs = 0;
    f_name = "f(x,y) = 1";
    f = f_0;
    content = 0;
    create_vector (M);
    recalc = 1;
    step_x = 1;
    step_y = 1;
}

void myGLWidget:: initializeGL ()
{
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    // отключает режим сглаживания цветов
    //glShadeModel(GL_FLAT);
    //glEnable(GL_CULL_FACE);
}

void myGLWidget:: resizeGL (int nWidth, int nHeight)
{
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();

    GLfloat ratio = (GLfloat)nHeight / (GLfloat)nWidth;

    if (nWidth >= nHeight)
        glOrtho (-1. / ratio, 1. / ratio, -1., 1., -10., 20. /* * sqrt (a * a + b * b + 2 * a * b * cos (alpha))*/);
    else
        glOrtho (-1., 1., -1. * ratio, 1. * ratio, -10., 20. /* * sqrt (a * a + b * b + 2 * a * b * cos (alpha))*/);

    // параметры видимости ортогональной проекции
    // плоскости отсечения (левая, правая, верхняя, нижняя, передняя, задняя)

    glViewport (0, 0, (GLint)nWidth, (GLint)nHeight);
}

void myGLWidget:: create_vector (int M)
{
    X = new double[M];
}
void myGLWidget:: free_vector ()
{
    delete [] X;
}

void myGLWidget:: set_X ()
{
    residual = approximation (f, hx, hy, nx, ny, M, p, X, eps);
}

void myGLWidget:: paintGL ()
{
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode (GL_MODELVIEW);

    glLoadIdentity ();

    glScaled (nSca, nSca, nSca);
    glTranslated (0.d, 0.d, 0.d);
    glRotated (xRot, 1.d, 0.d, 0.d);
    glRotated (yRot, 0.d, 1.d, 0.d);
    glRotated (zRot, 0.d, 0.d, 1.d);

    draw_axes ();
    draw_shadow (0.5, 0.5);
    draw_figure ();
    print_info ();
}

void myGLWidget:: draw_figure ()
{
    switch (content)
    {
        case func:
            paint_function (0.5, 0.5);
            break;
        case approx:
            paint_approximation (0.5, 0.5);
            break;
        case error_graph:
            paint_error (0.5, 0.5);
            break;
    }
}

void myGLWidget:: draw_shadow (double a, double b)
{
    double x, y;

    glBegin (GL_QUADS);
    glColor3d (0.8, 0.8, 0.8);
    x = 0.;
    y = 0.;
    x += y * cos (alpha * M_PI / 180.);
    y *= sin (alpha * M_PI / 180.);
    glVertex3d (x, y, 0.);
    x = 0.;
    y = b;
    x += y * cos (alpha * M_PI / 180.);
    y *= sin (alpha * M_PI / 180.);
    glVertex3d (x, y, 0.);
    x = a;
    y = b;
    x += y * cos (alpha * M_PI / 180.);
    y *= sin (alpha * M_PI / 180.);
    glVertex3d (x, y, 0.);
    x = a;
    y = 0;
    x += y * cos (alpha * M_PI / 180.);
    y *= sin (alpha * M_PI / 180.);
    glVertex3d (x, y, 0.);
    glEnd ();
}

void myGLWidget:: draw_axes ()
{
    //glLineWidth(2.d);
    glColor3d (0., 0., 0.);

    // OX
    glBegin (GL_LINES);
    glVertex3d (1., 0., 0.);
    glVertex3d (-1., 0., 0.);

    // стрелочки
    glVertex3d (0.97, 0.03, 0.);
    glVertex3d (1., 0., 0.);

    glVertex3d (0.97, -0.03, 0.);
    glVertex3d (1., 0., 0.);

    // буковка
    glVertex3d (1.03, 0.03, 0.);
    glVertex3d (1.06, -0.03, 0.);

    glVertex3d (1.03, -0.03, 0.);
    glVertex3d (1.06, 0.03, 0.);
    // OY
    glVertex3d (0., 1., 0.);
    glVertex3d (0., -1., 0.);

    // стрелочки
    glVertex3d (-0.03, 0.97, 0.);
    glVertex3d (0., 1., 0.);

    glVertex3d (0.03, 0.97, 0.);
    glVertex3d (0., 1., 0.);

    // буковка
    glVertex3d (-0.02, 1.01, 0.);
    glVertex3d (0.02, 1.05, 0.);

    glVertex3d (-0.02, 1.05, 0.);
    glVertex3d (0.04, 1.04, 0.);
    // OZ
    glVertex3d (0., 0., 0.95);
    glVertex3d (0., 0., -0.95);

    // стрелочки
    glVertex3d (0.02, 0., 0.93);
    glVertex3d (0., 0., 0.95);

    glVertex3d (-0.02, 0., 0.93);
    glVertex3d (0., 0., 0.95);

    // буковка
    glVertex3d (-0.02, 0., 0.96);
    glVertex3d (0.02, 0., 0.96);

    glVertex3d (-0.02, 0., 1.);
    glVertex3d (0.02, 0., 1.);

    glVertex3d (-0.02, 0., 0.96);
    glVertex3d (0.02, 0., 1.);

    glVertex3d (-0.015, 0., 0.98);
    glVertex3d (0.015, 0., 0.98);
    glEnd ();
}

void myGLWidget:: print_info ()
{
    char buf[60];

    switch (content)
    {
        case func:
            glColor3d (0., 0., 1.);
            renderText (width () / 10., height () / 10., "Function");
            break;
        case approx:
            glColor3d (0., 0.9, 0.);
            renderText (width () / 10., height () / 10., "Approximation");
            sprintf (buf, "Residual = %.3e", residual);
            renderText (width () / 10., 6 * height () / 20., buf);
            break;
        case error_graph:
            glColor3f (1.f, 0.6f, 0.0f);
            renderText (width () / 10., height () / 10., "Error");
            sprintf (buf, "Residual = %.3e", residual);
            renderText (width () / 10., 6 * height () / 20., buf);
            break;
    }

    glColor3d (0., 0., 0.);

    renderText (width () / 10., 3 * height () / 20., f_name);

    sprintf (buf, "nx = %d, ny = %d", nx, ny);
    renderText (width () / 10., 4 * height () / 20., buf);

    sprintf (buf, "max |F(x,y)| = %.3g", max_abs);
    renderText (width () / 10., 5 * height () / 20., buf);
}

void myGLWidget:: paint_function (double a, double b)
{
    int i, j;
    double	x, y, z;
    double x_min = 0, y_min = 0;

    max_abs = find_max_f();

    double coef = 4 * max_abs;

    for (j = 0; j < ny; j += step_y)
    {
        for (i = 0; i < nx; i += step_x)
        {
            glBegin (GL_TRIANGLES);
            //glColor3d (1. * (nx - i) / nx, 1. * j / ny, 0.);
            glColor3d (0., 0.5, 1.);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            //glColor3d (1. * (nx - i) / nx, 1. * j / ny, 0.5);
            glColor3d (0., 0., 1.);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            glEnd ();

            glBegin (GL_LINES);
            glColor3d (0.0, 0.0, 0.0);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = f(x, y) / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            glEnd ();
        }
    }
}

void myGLWidget::paint_approximation (double a, double b)
{
    int i, j;
    double	x, y, z;
    double x_min = 0, y_min = 0;

    if (recalc)
    {
        set_X ();
        recalc = 0;
        max_abs = find_max_X_i();
    }

    double coef = 4 * max_abs;

    for (j = 0; j < ny; j += step_y)
    {
        for (i = 0; i < nx; i += step_x)
        {
            glBegin (GL_TRIANGLES);
            //glColor3d (1. * (nx - i) / nx, 1. * j / ny, 0.);
            glColor3f(0.0f, 1.0f, 0.0f);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[i * (ny + 1) + j] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[i * (ny + 1) + j + step_y] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[(i + step_x) * (ny + 1) + j + step_y] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            //glColor3d (1. * (nx - i) / nx, 1. * j / ny, 0.5);
            glColor3f (0.f, 0.9f, 0.0f);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[i * (ny + 1) + j] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[(i + step_x) * (ny + 1) + j] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[(i + step_x) * (ny + 1) + j + step_y] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            glEnd ();
            glBegin (GL_LINES);
            glColor3d (0.0, 0.0, 0.0);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[i * (ny + 1) + j] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[i * (ny + 1) + j + step_y] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[(i + step_x) * (ny + 1) + j + step_y] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[i * (ny + 1) + j] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[(i + step_x) * (ny + 1) + j] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = X[(i + step_x) * (ny + 1) + j + step_y] / coef;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            glEnd ();
        }
    }
}

void myGLWidget::paint_error (double a, double b)
{
    int i, j;
    double	x, y, z;
    double x_min = 0, y_min = 0;
    double min = 0.01;

    if (recalc)
    {
        set_X ();
        recalc = 0;
    }

    double coef = 20.;

    for (j = 0; j < ny; j += step_y)
    {
        for (i = 0; i < nx; i += step_x)
        {
            glBegin (GL_TRIANGLES);
            //glColor3d (1. * (nx - i) / nx, 1. * j / ny, 0.);
            glColor3f (1.f, 0.8f, 0.0f);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[i * (ny + 1) + j]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[i * (ny + 1) + j + step_y]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[(i + step_x) * (ny + 1) + j + step_y]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            //glColor3d (1. * (nx - i) / nx, 1. * j / ny, 0.5);
            glColor4f (1.0f, 1.0f, 0.0f, 0.0f);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[i * (ny + 1) + j]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[(i + step_x) * (ny + 1) + j]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[(i + step_x) * (ny + 1) + j + step_y]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            glEnd ();
            glBegin (GL_LINES);
            glColor3d (0.0, 0.0, 0.0);
            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[i * (ny + 1) + j]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[i * (ny + 1) + j + step_y]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[(i + step_x) * (ny + 1) + j + step_y]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * i / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[i * (ny + 1) + j]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * j / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[(i + step_x) * (ny + 1) + j]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            x = a * (i + step_x) / nx;
            y = b * (j + step_y) / ny;
            x += y * cos (alpha * M_PI / 180.);
            y *= sin (alpha * M_PI / 180.);
            z = fabs(f(x,y) - X[(i + step_x) * (ny + 1) + j + step_y]) / coef;
            if (z <= 0 && z >= 0)
                z = min;
            x += x_min;
            y += y_min;
            glVertex3d (x, y, z);

            glEnd ();
        }
    }
}

double myGLWidget::find_max_f()
{
    int i, j;
    double x, y, z;
    double max = 0;

    max_abs = 0;

    for (j = 0; j <= ny; j++)
    {
        for (i = 0; i <= nx; i++)
        {
            x = a * i / nx;
            y = b * j / ny;
            z = f(x, y);

            if (fabs (z) > max)
                max = fabs (z);
        }
    }

    return max;
}

double myGLWidget::find_max_X_i()
{
    int j;
    double max = 0;

    max_abs = 0;

    for (j = 0; j < M; j++)
    {
        if (fabs (X[j]) > max)
            max = fabs (X[j]);
    }

    return max;
}

double myGLWidget::find_max_error()
{
    int i, j;
    double x, y, z;
    double min = 1;

    max_abs = 0;

    for (j = 0; j <= ny; j++)
    {
        for (i = 0; i <= nx; i++)
        {
            x = a * i / nx;
            y = b * j / ny;
            z = fabs(f(x,y) - X[i * (ny + 1) + j]);

            if (fabs (z) < min)
                min = fabs (z);
        }
    }

    return min;
}

void myGLWidget:: change_content ()
{
    s = 0;
    h = 0;
    content = (content + 1) % 3;

    switch (content)
    {
        case func:
            method_name = "Function";
            break;
        case approx:
            method_name = "Approximation";
            break;
        case error_graph:
            method_name = "Error graph";
            break;
    }

    //recalc = 1;
}

void myGLWidget:: change_func ()
{
    s = 0;
    h = 0;
    k = k % 8;

    switch (k)
    {
        case 0:
          f_name = "f(x,y) = 1";
          f = f_0;
          break;
        case 1:
          f_name = "f(x,y) = x";
          f = f_1;
          break;
        case 2:
          f_name = "f(x,y) = y";
          f = f_2;
          break;
        case 3:
          f_name = "f(x,y) = x + y";
          f = f_3;
          break;
        case 4:
          f_name = "f(x,y) = (x^2 + y^2)^0.5";
          f = f_4;
          break;
        case 5:
          f_name = "f(x,y) = x^2 + y^2";
          f = f_5;
          break;
        case 6:
          f_name = "f(x,y) = exp(x^2 - y^2)";
          f = f_6;
          break;
        case 7:
          f_name = "f(x,y) = 1 / (25 * (x^2 + y^2) + 1)";
          f = f_7;
          break;
    }

    k++;

    recalc = 1;
}

void myGLWidget:: scale_plus ()
{
    nSca = nSca * 1.1;

    s++;
}

void myGLWidget:: scale_minus ()
{
    nSca = nSca / 1.1;

    s--;
}

void myGLWidget:: rotate_left ()
{
    zRot += 15.;
}

void myGLWidget:: rotate_right()
{
    zRot -= 15.;
}

void myGLWidget:: p_plus_1()
{
    p++;

    recalc = 1;
}

void myGLWidget:: p_minus_1()
{
    p--;

    recalc = 1;
}

void myGLWidget:: nx_ny_div_2()
{
    free_vector ();

    nx /= 2;
    ny /= 2;

    step_x /= 2;
    step_y /= 2;

    if (nx < 5)
      nx = 5;

    if (ny < 5)
      ny = 5;

    if (step_x < 1)
      step_x = 1;

    if (step_y < 1)
      step_y = 1;

    hx = a / nx;
    hy = b / ny;

    M = (nx + 1) * (ny + 1);

    create_vector (M);

    recalc = 1;
}

void myGLWidget:: nx_ny_mult_2()
{
    free_vector ();

    nx *= 2;
    ny *= 2;

    step_x *= 2;
    step_y *= 2;

    if (nx == nx1) step_x = 1;
    if (ny == ny1) step_y = 1;

    hx = a / nx;
    hy = b / ny;

    M = (nx + 1) * (ny + 1);

    create_vector (M);

    recalc = 1;
}

void myGLWidget:: keyPressEvent (QKeyEvent* p)
{
    switch (p -> key())
    {
      case Qt::Key_0:
         change_func();
         break;
      case Qt::Key_1:
         change_content();
         break;
      case Qt::Key_2:
         scale_plus();
         break;
      case Qt::Key_3:
         scale_minus();
         break;
      case Qt::Key_4:
         nx_ny_mult_2();
         break;
      case Qt::Key_5:
         nx_ny_div_2();
         break;
      case Qt::Key_8:
         rotate_right();
         break;
      case Qt::Key_9:
         rotate_left();
         break;
      case Qt::Key_Escape:
         this -> close();
         break;
    }

    updateGL();
}

int parse (int argc, char**argv, double *a, double *b, double *alpha,
               int *nx, int *ny, int *k, double *eps, int *p)
{
    if (argc != 7)
        return -1;

    char buf[LEN];
    FILE *fp;
    char sep[] = " \n", *sym;
    bool h = 0;

    fp = fopen (argv[1], "r");
    if (!fp)
        return -2;

    while (fgets (buf, LEN, fp))
    {
        if (buf[0] < '0') continue;

        if (!h)
        {
            sym = strtok (buf, sep);
            *a = atof (sym);
            sym = strtok (NULL, sep);
            *b = atof (sym);

            h = 1;
        }
        else
        {
            sym = strtok (buf, sep);
            *alpha = atof (sym);
        }
    }

    *nx = atoi (argv[2]);
    *ny = atoi (argv[3]);
    *k = atoi (argv[4]);
    *eps = atof (argv[5]);
    *p = atoi (argv[6]);

    return 0;
}
