#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;



template <class T> 

AlgebraicVector<T> rhs (double t, AlgebraicVector<T> x)
{

    return AlgebraicVector<T> ( {x[1], -x[0]} ); // cos(t*(2*M_PI));


    /*AlgebraicVector<T> v(2);
    v[0] =  x[1];
    v[1] = -x[0]; // cos(t*(2*M_PI));
    return v;*/


}


template <class T>

AlgebraicVector<T> Euler (double t0, double tf, double h, AlgebraicVector<T> X0 )
{
    AlgebraicVector<T> X = X0;
    int steps = (tf-t0)/h;

    for(int i=0; i<steps; i+=1)
    {
        X = X + h*rhs(t0+h*i, X);
    }

    return X;

}


template <class T>
AlgebraicVector<T> ExactSolution (double t0, double tf, AlgebraicVector<T> X0 )
{
    return AlgebraicVector<T>( { cos(tf-t0)*X0[0]+sin(tf-t0)*X0[1], -sin(tf-t0)*X0[0]+cos(tf-t0)*X0[1] } );
}



int main( void )
{
    
    DA::init( 20, 2 ); // initialize DACE for 20th-order computations in 1 variable
    
    double eps = 1e-3;
    double h   = 0.0001;
    double t0  = 0;
    double tf  = 2*M_PI;
    
    AlgebraicVector<DA> x0 ( {1 + eps*DA(1), 1 + eps*DA(2)} ); // initialize x as DA
    AlgebraicVector<DA> x;

    AlgebraicVector<double> y0 = x0.eval<double>({ 1, 0.5 });
    AlgebraicVector<double> y;
   
    x = Euler (t0, tf, h, x0);
    y = Euler (t0, tf, h, y0);

    cout << tf  << endl << x0 << endl << x << endl;
    cout << vnorm( ExactSolution(t0, tf, x0) - x ) << endl;

    cout << vnorm(y-x.eval<double>({1,0.5})) << endl;
}
