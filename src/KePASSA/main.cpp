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
#include "EllipticOrbit.h"

#include "RK4.h"
#include "RK8.h"
#include "RK78.h"
#include "Adams.h"
#include "Leapfrog.h"
#include "Midpoint.h"
#include "BulirschStoer.h"
#include "EverhartRadau.h"

#include "test_expansion.h"

using namespace std;
using namespace DACE;

int main(void)
{
    // initialize DACE for 5th-order computations in 6 variables
    DA::init(6, 6);
    DA::setEps(0.0);  // required since we don't use scaled units

    cout.precision(16);
    cout.setf( ios::fixed, ios::floatfield );

    // Initial conditions (km, km/s), copied from the problem (StiefelScheiffele)
    AlgebraicVector<double> x0 = SS::x0;
    const double epsR = SS::epsR, epsV = SS::epsV;
    double t0 = SS::t0;
    double tf = SS::tf;
    AlgebraicVector<double> Ref_Sol = SS::Ref_Sol;

    // local settings, derived quantities
    //#define INTEGRATOR Leapfrog
    //#define INTEGRATOR RK4
    //#define INTEGRATOR Midpoint
    //#define INTEGRATOR Adams
    //#define INTEGRATOR BulirschStoer
    //#define INTEGRATOR RK8
    #define INTEGRATOR EverhartRadau
    double h = 100.;
    double tol = 1e-13;
    AlgebraicVector<DA> y0(6);
    y0[0] = x0[0] + epsR*DA(1);
    y0[1] = x0[1] + epsR*DA(2);
    y0[2] = x0[2] + epsR*DA(3);
    y0[3] = x0[3] + epsV*DA(4);
    y0[4] = x0[4] + epsV*DA(5);
    y0[5] = x0[5] + epsV*DA(6);

    // Perform propagations
    AlgebraicVector<double> x(6);
    vector<AlgebraicVector<double>> resx;

    cout << endl << "Propagating with 'double'..." << endl;
    funcs = 0;
    //x = EverhartRadauVar(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    //x = BulirschStoer<double, 2>(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    //x = Adams<double, 8, 1>(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    x = INTEGRATOR(t0, tf, h, tol, x0, StiefelScheiffele<double>, resx, false);
    cout << "  Done!" << endl;

    AlgebraicVector<DA> y(6);
    vector<AlgebraicVector<DA>> resy;

//    cout << endl << "Propagating with 'DA'..." << endl;
//    y = Leapfrog(t0, tf, h, tol, y0, StiefelScheiffele<DA>, resy, false);
//    cout << "  Done!" << endl;

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

    // test DA expansion
    test_expansion(INTEGRATOR, y0, t0, tf, h, tol, 10065);

    return 0;
}



int main_elliptic(void)
{
    // initialize DACE for 5th-order computations in 6 variables
    DA::init(3, 6);

    double rp  = 1;
    double ecc = 0.5;
    double t0  = 0.;
    double tf;
    double h   = 1e-2;
    double tol = 1e-6;

    AlgebraicVector<double> x0(4);
    AlgebraicVector<double> x(4);
    vector<AlgebraicVector<double>> resx;

    Elliptic_IC (rp, ecc, x0, tf);

    // Print Cartesian state vector at start of integration
    cout << endl;
    cout << " Initial Values:" << endl;
    cout << " x:  " << std::setprecision(16) << x0[0] << endl;
    cout << " y:  " << std::setprecision(16) << x0[1] << endl;
    cout << " vx: " << std::setprecision(16) << x0[2] << endl;
    cout << " vy: " << std::setprecision(16) << x0[3] << endl;
    cout << endl;

    // Perform integration
    x = Adams<double, 8, 1>(t0, tf, h, tol, x0, Elliptic_ODE<double>, resx, false);

    // Print Cartesian state vector at end of integration
    cout << " Final Values:" << endl;
    cout << " x:  " << std::setprecision(16) << x[0] << endl;
    cout << " y:  " << std::setprecision(16) << x[1] << endl;
    cout << " vx: " << std::setprecision(16) << x[2] << endl;
    cout << " vy: " << std::setprecision(16) << x[3] << endl;
    cout << endl;

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