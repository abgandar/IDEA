#include "propagator.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace DACE;




int main(void)
{

    // initialize DACE for 5th-order computations in 6 variables
    DA::init(5, 6);

    
    // Initial conditions
    AlgebraicVector<double> x0(6);
    
    //x0[0] = 1.;
    //x0[1] = 0.;
    //x0[2] = 0.;
    //x0[3] = 0.;
    //x0[4] = 1.;
    //x0[5] = 0.;

    x0[0] =     0.000000;
    x0[1] = -5888.972700;
    x0[2] = -3400.00000;
    x0[3] =    10.691338;
    x0[4] =     0.000000;
    x0[5] =     0.000000;

    AlgebraicVector<DA> y0(6);
    double eps = 1e-3;

    y0[0] = x0[0] + eps * DA(1);
    y0[1] = x0[1] + eps * DA(2);
    y0[2] = x0[2] + eps * DA(3);
    y0[3] = x0[3] + eps * DA(4);
    y0[4] = x0[4] + eps * DA(5);
    y0[5] = x0[5] + eps * DA(6);


    // Set up propagators and propagation parameters
    Cowell<double> PropDouble;
    Cowell<DA> PropDA;
    
    double h = 3;
    double t0 = 0.;
    //double tf = 288.12768941 * 24. * 3600.;
    double tf = 288.12768941 * 24. * 3600. / 100;


    // Perform propagations
    AlgebraicVector<double> x(6);
    double t_double;

    cout << endl << "Propagating with 'double'..." << endl;
    PropDouble.Propagate (x0.extract(0,2), x0.extract(3,5), t0, tf, h, x, t_double);
    cout << "  Done!" << endl;

    AlgebraicVector<DA> y(6);
    DA t_DA;
/*
    cout << endl << "Propagating with 'DA'..." << endl;
    PropDA.Propagate(y0.extract(0, 2), y0.extract(3, 5), t0, tf, h, y, t_DA);
    cout << "  Done!" << endl;
*/

    // Reference Position
    //AlgebraicVector<double> Ref_Sol ({-24219.05011593605201960070788, 227962.10637302200887202306088, 129753.44240008247047344318589});
    AlgebraicVector<double> Ref_Sol ({-0.006472054691587e5, 2.290183472899976e5, 1.322860577190283e5});

    // Print Cartesian state vector at end of integration
    cout << endl;
    cout << " t:  " << t_double << endl;
    cout << " x:  " << x[0] << endl;
    cout << " y:  " << x[1] << endl;
    cout << " z:  " << x[2] << endl;
    cout << " vx: " << x[3] << endl;
    cout << " vx: " << x[4] << endl;
    cout << " vx: " << x[5] << endl;
    cout << endl;
    cout << " t-Error: " << tf-t_double << endl;
    cout << " x-Error: " << 1e3 * vnorm(x.extract(0,2) - Ref_Sol) << " m" << endl;
    cout << endl;
  /*  
    Interv I = Interv(y, 0);
    I = Interv(y, 0);   cout << " x:  " << I.Mean << "\t[" << I.Lower << ", " << I.Upper << "]\t w:" << I.Width << endl;
    I = Interv(y, 1);   cout << " y:  " << I.Mean << "\t[" << I.Lower << ", " << I.Upper << "]\t w:" << I.Width << endl;
    I = Interv(y, 2);   cout << " z:  " << I.Mean << "\t[" << I.Lower << ", " << I.Upper << "]\t w:" << I.Width << endl;
    I = Interv(y, 3);   cout << " vx: " << I.Mean << "\t[" << I.Lower << ", " << I.Upper << "]\t w:" << I.Width << endl;
    I = Interv(y, 4);   cout << " vy: " << I.Mean << "\t[" << I.Lower << ", " << I.Upper << "]\t w:" << I.Width << endl;
    I = Interv(y, 5);   cout << " vz: " << I.Mean << "\t[" << I.Lower << ", " << I.Upper << "]\t w:" << I.Width << endl;
    cout << endl;
*/
}