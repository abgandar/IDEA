#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;

int main( void )
{
    
    DA::init( 20, 2 ); // initialize DACE for 20th-order computations in 1 variable
    
    DA x = DA(1)+ DA(2); // initialize x as DA
   
    DA y = sin(x); // compute y = sin(x)
    
    // print x and y to screen
    cout << "x" << endl << x << endl;
    cout << "x_1" << endl << DA(1) << endl;
    cout << "x_2" << endl << DA(2) << endl;
    cout << "y = sin(x)" << endl << y;
    cout << "y_{2,3}" << endl << y.getCoefficient({2, 3}) << endl;
    cout << "y(1,2)" << endl << y.eval<double>({2., 1.}) << endl;   // note: the <double> is required only when using the curly braces syntax
    //cout << "y(1,2)" << endl << y.eval<DA>({2., 1.}) << endl;
    //AlgebraicVector<double> p = {2., 1.};  cout << "y(1,2)" << endl << y.eval(p) << endl;
    
}
