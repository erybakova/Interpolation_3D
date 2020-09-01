#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <qgl.h>

#include "glwidget.h"
#include "approximation.h"
#include "init.h"
#include "solve.h"

int main(int argc, char** argv)
{
    double a, b, alpha, eps;
    int nx, ny, k, p;

    int res = parse (argc, argv, &a, &b, &alpha, &nx, &ny, &k, &eps, &p);

    if (res < 0)
    {
        switch (res)
        {
            case -1:
              printf("Usage: ./a.out input.txt nx ny k eps p\n");
              return 0;
            case -2:
              printf("Cannot open input.txt!\n");
              return 0;
        }
    }

    //printf ("a %.3f, b %.3f, k %d, eps %e\n", a, b, k, eps);

    QApplication app(argc, argv);

    myGLWidget my_widget;

    my_widget.init(a, b, alpha, nx, ny, k, eps, p);
    app.setActiveWindow(&my_widget);

    my_widget.show();
    app.exec();

    my_widget.free_vector ();

    return 0;
}
