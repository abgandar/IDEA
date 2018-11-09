#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;



template <class T> 

T rhs (double t, T x)
{

    return x; // cos(t*(2*M_PI));

}


template <class T>

T Euler (double t0, double tf, double h, T X0 )
{
    T X = X0;
    int steps = (tf-t0)/h;

    for(int i=0; i<steps; i+=1)
    {
        X = X + h*rhs(t0+h*i, X);
    }

    return X;

}




int main( void )
{
    
    DA::init( 20, 1 ); // initialize DACE for 20th-order computations in 1 variable
    
    double eps = 1e-3;
    double h   = 0.00001;
    double t0  = 0;
    double tf  = 1;
    
    DA x0 = 1 + eps*DA(1); // initialize x as DA
    DA x;

    double y0 = cons(x0);
    double y;
   
    x = Euler (t0, tf, h, x0);
    y = Euler (t0, tf, h, y0);

    cout << tf  << endl << x0 << endl << x << endl;
    cout << x0 * exp(tf) << endl;
    cout << norm( x0*exp(tf) - x, 2 ) << endl;
    cout << y << endl;
    
    
}
