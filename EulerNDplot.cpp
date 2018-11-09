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

    return AlgebraicVector<T> ( {x[1], -x[0]} ) * vnorm(x); // cos(t*(2*M_PI));


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



void Draw (AlgebraicVector<DA> v, int N=100)
{
    AlgebraicVector<double> a(N), bp(N,1.), bm(N,-1.), x, y;
    for (int i=0; i<N; i++){
        a[i] = -1. + i*2./(N-1);
    }

    x.insert( x.end(),  a.begin(),  a.end() );
    y.insert( y.end(), bp.begin(), bp.end() );

    x.insert( x.end(), bp.begin(), bp.end() );
    a = -a;
    y.insert( y.end(), a.begin(),  a.end() );
    
    x.insert( x.end(),  a.begin(),  a.end() );
    y.insert( y.end(), bm.begin(), bm.end() );
    a = -a;

    x.insert( x.end(), bm.begin(), bm.end() );
    y.insert( y.end(),  a.begin(),  a.end() );

    auto r = v.eval<AlgebraicVector<double>> ({x, y});

    for(int n = 0; n < r[0].size(); n++)
        cout << r[0][n] << "   " << r[1][n] << endl;
    cout << endl;



}



int main( void )
{
    
    DA::init( 3, 2 ); // initialize DACE for 20th-order computations in 1 variable
    
    double eps = 3e-1;
    double h   = 0.0001;
    double t0  = 0;
    AlgebraicVector<double> tf ( {0.2, 0.4, 0.6, 0.8, 1.} ); tf *= 4; //2 * M_PI;

    AlgebraicVector<DA> x0 ( {1 + eps*DA(1), 1 + eps*DA(2)} ); // initialize x as DA
    AlgebraicVector<DA> x = x0;

    for(int i=0; i<tf.size(); i++)
    {
        x  = Euler (t0, tf[i], h, x);
        t0 = tf[i];

        Draw(x);
    }
    
    //cout << tf  << endl << x0 << endl << x << endl;
    
}