#include <dace/dace.h>

using namespace DACE;
using namespace std;



class Interv
{
  public:
    double Lower;
    double Upper;
    double Mean;
    double Width;

    // Constructor
    Interv(AlgebraicVector<DA> y, int idx);
};



Interv::Interv(AlgebraicVector<DA> y, int idx)
{
    double val0, valm, valp;

    AlgebraicVector<double> zeros(y.size(), 0.);
    AlgebraicVector<double> ones(y.size(), 1.);

    val0 = y[idx].eval<double>(zeros);
    valm = y[idx].eval<double>(-ones);
    valp = y[idx].eval<double>(ones);

    Lower = min(valm, valp);
    Upper = max(valm, valp);
    Mean = val0;
    Width = Upper - Lower;
}