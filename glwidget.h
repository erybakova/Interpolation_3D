#ifndef _my_widget
#define _my_widget

#include <qgl.h>
#include <QKeyEvent>

#include <math.h>
#include <qnamespace.h>

enum content
{
    func,
    approx,
    error_graph
};

class myGLWidget : public QGLWidget
{
  Q_OBJECT

  private:
    double a, b;
    int p;
    int nx, ny;
    double hx, hy;
    double alpha;
    int M;          // число точек сетки
    double *X;      // массив с коэффициентами
    double eps;     // точность решения системы
    int k;          // номер функции
    int s;          // текущий масштаб
    int h;          // текущее возмущение
    int content;    // текущее содержимое виджета
    double max_abs; // максимум модуля
    const char *f_name;
    const char *method_name;
    double (*f) (double, double);
    bool recalc;
    int step_x, step_y;
    int nx1, ny1;
    double residual;

    double xRot; // угол поворота вокруг оси X
    double yRot; // угол поворота вокруг оси Y
    double zRot; // угол поворота вокруг оси Z
    double nSca; // масштабирование объекта

  protected:
    void initializeGL ();
    void resizeGL (int nWidth, int nHeight);
    void paintGL ();
    void keyPressEvent (QKeyEvent* p);

  public:
    myGLWidget (QWidget* parent = 0);
    void init (double a_, double b_, double alpha_,
                 int nx_, int ny_, int k_, double eps_, int p_);
    void set_X ();
    double find_max_f ();
    double find_max_X_i ();
    double find_max_error ();
    void print_info ();
    void create_vector (int M);
    void free_vector ();
    void paint_function (double a, double b);
    void paint_approximation (double a, double b);
    void paint_error (double a, double b);

    void scale_plus ();
    void scale_minus ();
    void rotate_left ();
    void rotate_right ();
    void draw_axes ();
    void draw_figure ();
    void draw_shadow (double a, double b);

  public slots:
    void change_func ();
    void change_content ();
    void nx_ny_div_2 ();
    void nx_ny_mult_2 ();
    void p_plus_1 ();
    void p_minus_1 ();
};

int parse (int argc, char**argv, double *a, double *b, double *alpha,
               int *nx, int *ny, int *k, double *eps, int *p);

#endif
