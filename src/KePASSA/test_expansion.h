#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <dace/dace.h>

#include "StiefelScheiffele.h"
//#include "VanDerPol.h"
//#include "EllipticOrbit.h"

#include "RK4.h"
#include "RK8.h"
#include "RK78.h"
#include "Adams.h"
#include "Leapfrog.h"
#include "Midpoint.h"
#include "BulirschStoer.h"
#include "EverhartRadau.h"

using namespace std;
using namespace DACE;

// borrowed and improved from the intertubes (and checked to be correct, but not very elegant, stable, or efficient)
template <int n>
double determinant( const double matrix[n][n] )
{
    double det = 0.0;
    double submatrix[n-1][n-1];

    for( int x = 0; x < n; x++ )
    {
        int subi = 0; 
        for( int i = 1; i < n; i++ )
        {
            int subj = 0;
            for( int j = 0; j < n; j++ )
            {
                if (j == x) continue;
                submatrix[subi][subj] = matrix[i][j];
                subj++;
            }
            subi++;
        }
        det += ((x%2) ? -1.0 : 1.0) * matrix[0][x] * determinant( submatrix );
    }
    return det;
}

template <>
double determinant<2>( const double matrix[2][2] )
{
    return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
}

// compute determinant at given point
double det(const compiledDA &cda, vector<double> pt)
{
    auto res = cda.eval(pt);
    double mat[6][6];
    for( unsigned int i = 0; i < 36; i++ )
        mat[i/6][i%6] = res[i];
    return determinant(mat);
}

void test_expansion(AlgebraicVector<DA> (*integrate)(double, double, double, double, AlgebraicVector<DA>, AlgebraicVector<DA> (*)(AlgebraicVector<DA>, double), vector<AlgebraicVector<DA>>&, bool), AlgebraicVector<DA> x0, const double t0, const double tf, const double h, const double tol, const unsigned int N, const unsigned int Nskip = 0)
{
    // test accuracy of expansion by integrating DA for set time, then evaluate a cloud of N points, and compare to "correct" results
    vector<AlgebraicVector<DA>> resx;
    auto x = integrate(t0, tf, h, tol, x0, StiefelScheiffele<DA>, resx, false);
    auto xc = x.compile();
    AlgebraicVector<DA> Dx(36);
    for( unsigned int i = 0; i < 6; i++ )
        for( unsigned int j = 0; j < 6; j++ )
            Dx[i*6+j] = x[i].deriv(j+1)/(j < 3 ? SS::epsR : SS::epsV);    // variable numbers are 1 based, don't forget to undo scaling
    auto Dxc = Dx.compile();
    // read N data points from file and evaluate at each
    vector<double> pt(6), res(6), diff(6);
    double maxdet = 0.0;
    ifstream is("points.dat"), iref("reference.dat");
    ofstream os("results.data"), odiff("diff.data"), odet("det.data");
    os.precision(16);
    os.setf( std::ios::fixed, std::ios::floatfield );
    odiff.precision(16);
    odiff.setf( std::ios::fixed, std::ios::floatfield );
    odet.precision(16);
    odet.setf( std::ios::fixed, std::ios::floatfield );
    // skip first Nskip points
    for( unsigned int i = 0; i < 6*Nskip; i++ )
    {
        double temp;
        is >> temp;
    }

    // evaluate all others and write results to file
    for( unsigned int i = 0; i < N; i++ )
    {
        for( unsigned int j = 0 ; j < 6; j++ )
            is >> pt[j];
        res = xc.eval(pt);
        for( unsigned int j = 0 ; j < 6; j++ )
        {
            double ref;
            iref >> ref;
            ref *= 1e-3;    // reference in m => km
            os << res[j] << endl;
            odiff << res[j]-ref << endl;
            diff[j] = max(fabs(res[j]-ref), diff[j]);
        }
        const double d = det(Dxc, pt);
        odet << d << endl;
        maxdet = max(fabs(1.0-d), maxdet);
    }
    // print results
    cout << "Maximum error in each component [m, m/s]:" << endl << diff[0]*1e3 << ", \n" << diff[1]*1e3 << ", \n" << diff[2]*1e3 << ", \n" << diff[3]*1e3 << ", \n" << diff[4]*1e3 << ", \n" << diff[5]*1e3 << endl;
    cout << "Maximum error in determinant:" << endl << maxdet << endl;
}
