#include <dace/dace.h>
using namespace std;
using namespace DACE;

template <class T>
AlgebraicVector<T> VanDerPol(AlgebraicVector<T> x, double t)
{
    AlgebraicVector<T> res(2);
    res[0] = x[1];
    res[1] = ((1-x[0]*x[0])*x[1]-x[0])/0.1;
    return res;
}

