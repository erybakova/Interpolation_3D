# Interpolation_3D
Represents an approximation of a function of 2 variables x and y, defined into the parallelogram, by building the grid and defining Courant functions in the vertices.

How to set the approximation area (parallelogram):

For simplicity let`s assume that the parallelogram is located in the first coordinate quarter, the lower base lies on the axis of the abscissus 
and the lower left vertex coincides with the origin.

Let 
    a and b be the lengths of parallelogram sides; 
    nx and ny - the number of split segments on the sides, respectively (so we get the grid of size nx*ny);
    alpha - degree measure of the angle between the sides a and b (from 0 to pi/2).
    
Input: ./a.out input.txt nx ny k eps p

Here 

     input.txt - file that defines a, b and alpha;
     k - the number of the initial approximate function;
     eps - accuracy;
     p - the number of threads;
     
An example of running: ./a.out input.txt 5 5 0 1e-14 4

Available values of k:

k = 0: f(x, y) = 1
k = 1: f(x, y) = x
k = 2: f(x, y) = y
k = 3: f(x, y) = x + y
k = 4: f(x) = √(x^2 + y^2)
k = 5: f(x) = x^2 + y^2
k = 6: f(x) = e^(x^2 - y^2)
k = 7: f(x) = 1/(25(x^2 + y^2) + )
