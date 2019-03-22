/*
    ToDo:
    - Check RK78 step-size control
    - Integrators:
      * Taylor (Alex)
      * Adams family
      * Everhart
    - Problems:
      * elliptic, circular
      * 3 body (decide on test case)
    - Measures:
      * Automate code for order plotting
      * ...
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <dace/dace.h>

#include "StiefelScheiffele.h"
#include "VanDerPol.h"

#include "RK4.h"
#include "RK8.h"
#include "RK78.h"
#include "Leapfrog.h"
#include "Midpoint.h"
#include "BulirschStoer.h"
#include "EverhartRadau.h"

using namespace std;
using namespace DACE;

int main(void)
{
    // initialize DACE for 5th-order computations in 6 variables
    DA::init(3, 6);

    // Initial conditions (km, km/s)
    AlgebraicVector<double> x0 = { 0.000000, -5888.972700, -3400.00000, 10.691338, 0.000000, 0.000000 };
    const double eps = 1e-1;
    AlgebraicVector<DA> y0 = x0 + eps * AlgebraicVector<DA>::identity(6);

    double h = 100;
    double tol = 1e-8;
    double t0 = 0.;
    double tf = 288.12768941 * 24. * 3600.;

    // Perform propagations
    AlgebraicVector<double> x(6);
    vector<AlgebraicVector<double>> resx;

    cout << endl << "Propagating with 'double'..." << endl;
    x = EverhartRadauVar(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    //x = BulirschStoer<double, 2>(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    //x = RK8(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    cout << "  Done!" << endl;

    AlgebraicVector<DA> y(6);
    vector<AlgebraicVector<DA>> resy;

    cout << endl << "Propagating with 'DA'..." << endl;
//    y = Leapfrog(t0, tf, h, tol, y0, StiefelScheiffele<DA>, resy, false);
    cout << "  Done!" << endl;

    // Reference Position
    AlgebraicVector<double> Ref_Sol ({-24219.05011593605201960070788, 227962.10637302200887202306088, 129753.44240008247047344318589});
    //AlgebraicVector<double> Ref_Sol ({-0.006472054691587e5, 2.290183472899976e5, 1.322860577190283e5});

    // Print Cartesian state vector at end of integration
    cout << endl;
    cout << " x:  " << x[0] << endl;
    cout << " y:  " << x[1] << endl;
    cout << " z:  " << x[2] << endl;
    cout << " vx: " << x[3] << endl;
    cout << " vx: " << x[4] << endl;
    cout << " vx: " << x[5] << endl;
    cout << endl;
    cout << " x-Error: " << 1e3 * vnorm(x.extract(0,2) - Ref_Sol) << " m" << endl;
    cout << endl;
    cout << "Number of fixed steps expected: " << int((tf-t0)/h)+1 << endl;
    cout << "Number of function evaluations: " << funcs << endl;
    /*cout << endl;
    cout << " x:  " << y[0] << endl;
    cout << " y:  " << y[1] << endl;
    cout << " z:  " << y[2] << endl;
    cout << " vx: " << y[3] << endl;
    cout << " vy: " << y[4] << endl;
    cout << " vz: " << y[5] << endl;
    cout << endl;*/

    return 0;
}


int main_not_so_much(void)
{
    // Initial conditions (km, km/s)
    AlgebraicVector<double> x0 = { 2.0, 0.0 };

    double h = 0.005;
    double tol = 0.0001;
    double t0 = 0.;
    double tf = 0.01;

    // Perform propagations
    AlgebraicVector<double> x(2);
    vector<AlgebraicVector<double>> resx;

    cout << endl << "Propagating with 'double'..." << endl;
    x = BulirschStoer<double, 4>(t0, tf, h, tol, x0, VanDerPol<double>, resx, false);
//    x = RK78(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    cout << "  Done!" << endl;

    // Print Cartesian state vector at end of integration
    cout << endl;
    cout << " x: " << std::setprecision(16) << x[0] << endl;
    cout << " y: " << std::setprecision(16) << x[1] << endl;
    cout << endl;
    cout << "Number of fixed steps expected: " << int((tf-t0)/h)+1 << endl;

    return 0;
}