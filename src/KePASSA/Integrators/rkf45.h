template <class T>
AlgebraicVector<T> RKF45(double t0, double tf, double h, const AlgebraicVector<T> &X0, Propagator<T> &P)
{

    const double C[6]    = { 0., 1./5., 3./10., 3./5., 1., 7./8. };
    const double B[6]    = { 37./378., 0., 250./621., 125./594., 0., 512./1771. };
    const double Bm[6]   = { 2825./27648., 0., 18575./48384., 13525./55296., 277./14336., 1./4. };
    const double A[6][5] = {{0., 0., 0., 0., 0.}, 
                            {1./5.,0., 0., 0., 0.},
                            {3./40., 9./40., 0., 0., 0.},
                            {3./10., -9./10., 6./5., 0., 0.},
                            {-11./54., 5./2., -70./27., 35./27., 0.},
                            {1631./55296., 175./512., 575./13824., 44275./110592., 253./4096.}};
    

    unsigned int n = X0.size();
    unsigned int maxsteps  = (tf - t0) / h;
    
    AlgebraicVector<T> X(n);
    AlgebraicVector<T> Xk(n);
    AlgebraicVector<T> K0(n);
    AlgebraicVector<T> K1(n);
    AlgebraicVector<T> K2(n);
    AlgebraicVector<T> K3(n);
    AlgebraicVector<T> K4(n);
    AlgebraicVector<T> K5(n);
    
    double t;
    double progress;

    Xk = X0;

    for (int steps = 0; steps < maxsteps; steps++)
    {

        // Update value of independent variables
        t = t0 + h * steps;


        // Evaluate Ki for every intermediate stage
        X = Xk;

        //>> Stage 1
        K0 = h * P.RHS(t + C[0]*h, X);

        //>> Stage 2
        X = Xk + A[1][0] * K0;
        K1 = h * P.RHS(t + C[1]*h, X);
        
        //>> Stage 3
        X = Xk + A[2][0] * K0 + A[2][1] * K1;
        K2 = h * P.RHS(t + C[2]*h, X);
        
        //>> Stage 4
        X = Xk + A[3][0] * K0 + A[3][1] * K1 + A[3][2] * K2;
        K3 = h * P.RHS(t + C[3]*h, X);
        
        //>> Stage 5
        X = Xk + A[4][0] * K0 + A[4][1] * K1 + A[4][2] * K2 + A[4][3] * K3;
        K4 = h * P.RHS(t + C[4]*h, X);
        
        //>> Stage 6
        X = Xk + A[5][0] * K0 + A[5][1] * K1 + A[5][2] * K2 + A[5][3] * K3 + A[5][4] * K4;
        K5 = h * P.RHS(t + C[5]*h, X);


        // Estimate y(x+h) of the highest order available
        Xk = Xk + B[0] * K0
                + B[1] * K1
                + B[2] * K2
                + B[3] * K3
                + B[4] * K4
                + B[5] * K5;
        

        progress = 100 * (t0 + h * steps) / (tf-t0); 
        if( steps % 10000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return Xk;
};
