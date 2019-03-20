#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace DACE;


template <class T>
class Propagator
{
  public:
    double Time0;

    const double GM = 3.98601e5;
    const double Lc = 1.360004184565669e+05;
    const double Wc = 1.258806060043139e-05;

    virtual AlgebraicVector<T> RHS(double IndepVar, const AlgebraicVector<T> &DepVars) = 0;
    virtual void Setup(const AlgebraicVector<T> &Pos, const AlgebraicVector<T> &Vel, double &Time, AlgebraicVector<T> &DepVars, T IndepVar) = 0;
    virtual void Cartesian(const AlgebraicVector<T> &DepVars, double &IndepVar, AlgebraicVector<T> &State, T &Time) = 0;
    virtual void Propagate(const AlgebraicVector<T> &R0, const AlgebraicVector<T> &V0, double &t0, double tf, double h, AlgebraicVector<T> &State, T &Time) = 0;
};


// Auxiliary:
#include "interv.h"


// Formulations:
#include "Cowell.h"


// Integrators:
#include "euler.h"
#include "rkf45.h"
#include "rkf78.h"
