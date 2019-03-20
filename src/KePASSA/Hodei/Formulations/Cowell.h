#include <dace/dace.h>

using namespace DACE;




template <class T>
class Cowell : public Propagator<T>
{
  public:
    AlgebraicVector<T> RHS(double IndepVar, const AlgebraicVector<T> &DepVars);
    void Setup(const AlgebraicVector<T> &Pos, const AlgebraicVector<T> &Vel, double &Time, AlgebraicVector<T> &DepVars, T IndepVar);
    void Cartesian(const AlgebraicVector<T> &DepVars, double &IndepVar, AlgebraicVector<T> &State, T &Time);
    void Propagate(const AlgebraicVector<T> &R0, const AlgebraicVector<T> &V0, double &t0, double tf, double h, AlgebraicVector<T> &State, T &Time);
};




template <class T>
AlgebraicVector<T> Cowell<T>::RHS(double IndepVar, const AlgebraicVector<T> &DepVars)
{
    AlgebraicVector<T> deriv(6);
    AlgebraicVector<T> State(6);
    AlgebraicVector<T> ap(3);
    T Time;

    // Recover dimensional Cartesian State at current integration step
    Cowell<T>::Cartesian (DepVars, IndepVar, State, Time);

    // Evaluate dimensional perturbing accelerations
    //ap = SS2B_perturbation (State.extract(0,2), State.extract(3,5), Time);

//----------
    AlgebraicVector<T> R = State.extract(0,2);
    AlgebraicVector<T> V = State.extract(3,5);

    // J2 Perturbation
    T Rho = vnorm(R);
    
    const double R_p = 6371.2200;
    const double J2 = 0.00108265;
    const double GM = 3.98601e5;

    T cost_J2 = -1.5 * GM * J2 * R_p*R_p / pow(Rho, 7);
    T Rho2 = Rho * Rho;
    T z2 = R[2] * R[2];

    AlgebraicVector<T> Contribution(3);

    Contribution[0] = cost_J2 * R[0] * (    Rho2 - 5 * z2);
    Contribution[1] = cost_J2 * R[1] * (    Rho2 - 5 * z2);
    Contribution[2] = cost_J2 * R[2] * (3 * Rho2 - 5 * z2);

    ap = ap + Contribution;

    // Lunar Perturbation
    const double mu_m = 4.90266e3;
    const double d_m = 384400;
    const double w_m = 2.665315780887e-6;

    AlgebraicVector<T> RM(3);
    AlgebraicVector<T> RMrel(3);

    RM[0] =  d_m * sin(w_m * Time);
    RM[1] = -d_m * 0.5 * sqrt(3.) * cos(w_m * Time);
    RM[2] = -d_m * 0.5 * cos(w_m * Time);

    RMrel = R - RM;
    
    T Rho_M = vnorm(RMrel);
    T Rho_M3 = Rho_M * Rho_M * Rho_M;
    double d_m3 = d_m * d_m * d_m;

    Contribution[0] = -mu_m * (RMrel[0] / Rho_M3 + RM[0] / d_m3);
    Contribution[1] = -mu_m * (RMrel[1] / Rho_M3 + RM[1] / d_m3);
    Contribution[2] = -mu_m * (RMrel[2] / Rho_M3 + RM[2] / d_m3);

    ap = ap + Contribution;

//----------

    // Make perturbing accelerations non-dimensional
    ap = ap / ( Cowell<T>::Lc * Cowell<T>::Wc * Cowell<T>::Wc );

    const double mu = 1.;

    T r = vnorm(DepVars.extract(0,2));
    T aux = mu / (r*r*r);

    deriv[0] =  DepVars[3];
    deriv[1] =  DepVars[4];
    deriv[2] =  DepVars[5];
    deriv[3] = -aux * DepVars[0] + ap[0];
    deriv[4] = -aux * DepVars[1] + ap[1];
    deriv[5] = -aux * DepVars[2] + ap[2];

    return deriv;
}

template <class T>
void Cowell<T>::Setup(const AlgebraicVector<T> &Pos, const AlgebraicVector<T> &Vel, double &Time, AlgebraicVector<T> &DepVars, T IndepVar)
{
    
    // Get Semi-Major Axis
    //T Energy =  0.5 * vnorm(Vel) * vnorm(Vel) - Cowell<T>::GM / vnorm(Pos);
    //T SMA    = -0.5 * Cowell<T>::GM / Energy;

    // Characteristic Values
    //Cowell<T>::Lc = abs(SMA);
    //Cowell<T>::Wc = sqrt(Cowell<T>::GM / (Cowell<T>::Lc * Cowell<T>::Lc * Cowell<T>::Lc) );
    //Cowell<T>::Time0 = Time;

    // Characteristic Values
    //AlgebraicVector<double> zeros(2*Pos.size(), 0.);
    
    //Cowell<T>::Lc = 1.; //vnorm(Pos.eval<double>(zeros));
    //Cowell<T>::Wc = 1.; //sqrt(Cowell<T>::GM / (Cowell<T>::Lc * Cowell<T>::Lc * Cowell<T>::Lc));
    
    Cowell<T>::Time0 = Time;

    IndepVar = 0.;
    
    DepVars[0] = Pos[0] /  Cowell<T>::Lc;
    DepVars[1] = Pos[1] /  Cowell<T>::Lc;
    DepVars[2] = Pos[2] /  Cowell<T>::Lc;
    DepVars[3] = Vel[0] / (Cowell<T>::Lc * Cowell<T>::Wc);
    DepVars[4] = Vel[1] / (Cowell<T>::Lc * Cowell<T>::Wc);
    DepVars[5] = Vel[2] / (Cowell<T>::Lc * Cowell<T>::Wc);

}

template <class T>
void Cowell<T>::Cartesian(const AlgebraicVector<T> &DepVars, double &IndepVar, AlgebraicVector<T> &State, T &Time)
{

    State[0] = DepVars[0] * Cowell<T>::Lc;
    State[1] = DepVars[1] * Cowell<T>::Lc;
    State[2] = DepVars[2] * Cowell<T>::Lc;
    State[3] = DepVars[3] * Cowell<T>::Lc * Cowell<T>::Wc;
    State[4] = DepVars[4] * Cowell<T>::Lc * Cowell<T>::Wc;
    State[5] = DepVars[5] * Cowell<T>::Lc * Cowell<T>::Wc;

    Time = Cowell<T>::Time0 + IndepVar / Cowell<T>::Wc;
}

template <class T>
void Cowell<T>::Propagate(const AlgebraicVector<T> &R0, const AlgebraicVector<T> &V0, double &t0, double tf, double h, AlgebraicVector<T> &State, T &Time)
{
    double IndepVar0;
    double IndepVarf = IndepVar0 + (tf-t0) * Cowell<T>::Wc;
    double Stepsize  = h * Cowell<T>::Wc;

    AlgebraicVector<T> DepVars0(6);
    AlgebraicVector<T> DepVars(6);

    // Get initial Cartesian state vector in non-dimensional form
    Cowell<T>::Setup(R0, V0, t0, DepVars0, IndepVar0);
    
    //cout << endl;
    //cout << " x:  " << DepVars0[0] << endl;
    //cout << " y:  " << DepVars0[1] << endl;
    //cout << " z:  " << DepVars0[2] << endl;
    //cout << " vx: " << DepVars0[3] << endl;
    //cout << " vx: " << DepVars0[4] << endl;
    //cout << " vx: " << DepVars0[5] << endl;
    //cout << endl;
    
    // Perform numerical integration of the eqs. of motion
    DepVars = RKF45(IndepVar0, IndepVarf, Stepsize, DepVars0, *this);

    // Recover Cartesian State at the end of integration
    Cowell<T>::Cartesian(DepVars, IndepVarf, State, Time);

}
