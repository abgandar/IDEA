template <class T>
AlgebraicVector<T> RKF78(double t0, double tf, double h, const AlgebraicVector<T> &X0, Propagator<T> &P)
{

    const double C[13]     = { 0., 2./27., 1./9., 1./6., 5./12., 1./2., 5./6., 1./6., 2./3., 1./3., 1., 0., 1. };
    const double B[13]     = { 0., 0., 0., 0., 0., 34./105., 9./35., 9./35., 9./280., 9./280., 0., 41./840., 41./840. };
    const double Bm[13]    = { 41./840., 0., 0., 0., 0., 34./105., 9./35., 9./35., 9./280., 9./280., 41./840., 0., 0. };
    const double A[13][12] = {{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
                             {2./27., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
                             {1./36., 1./12., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
                             {1./24., 0., 1./8., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
                             {5./12., 0., -25./16., 25./16., 0., 0., 0., 0., 0., 0., 0., 0.}, 
                             {1./20., 0., 0., 1./4., 1./5., 0., 0., 0., 0., 0., 0., 0.}, 
                             {-25./108., 0., 0., 125./108., -65./27., 125./54., 0., 0., 0., 0., 0., 0.}, 
                             {31./300., 0., 0., 0., 61./225., -2./9., 13./900., 0., 0., 0., 0., 0.}, 
                             {2., 0., 0., -53/6., 704./45., -107/9., 67/90., 3., 0., 0., 0., 0.}, 
                             {-91./108., 0., 0., 23./108., -976./135., 311./54., -19./60., 17./6., -1./12., 0., 0., 0.}, 
                             {2383./4100., 0., 0., -341./164., 4496./1025., -301./82., 2133./4100., 45./82., 45./164., 18./41., 0., 0.}, 
                             {3./205., 0., 0., 0., 0., -6./41., -3./205., -3./41., 3./41., 6./41., 0., 0.}, 
                             {-1777./4100., 0., 0., -341./164., 4496./1025., -289./82., 2193./4100., 51./82., 33./164., 12./41., 0., 1.},  };    


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
    AlgebraicVector<T> K6(n);
    AlgebraicVector<T> K7(n);
    AlgebraicVector<T> K8(n);
    AlgebraicVector<T> K9(n);
    AlgebraicVector<T> K10(n);
    AlgebraicVector<T> K11(n);
    AlgebraicVector<T> K12(n);
    
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

        //>> Stage 7
        X = Xk + A[6][0] * K0 + A[6][1] * K1 + A[6][2] * K2 + A[6][3] * K3 + A[6][4] * K4
               + A[6][5] * K5;
        K6 = h * P.RHS(t + C[6]*h, X);

        //>> Stage 8
        X = Xk + A[7][0] * K0 + A[7][1] * K1 + A[7][2] * K2 + A[7][3] * K3 + A[7][4] * K4
               + A[7][5] * K5 + A[7][6] * K6;
        K7 = h * P.RHS(t + C[7]*h, X);

        //>> Stage 9
        X = Xk + A[8][0] * K0 + A[8][1] * K1 + A[8][2] * K2 + A[8][3] * K3 + A[8][4] * K4
               + A[8][5] * K5 + A[8][6] * K6 + A[8][7] * K7;
        K8 = h * P.RHS(t + C[8]*h, X);

        //>> Stage 10
        X = Xk + A[9][0] * K0 + A[9][1] * K1 + A[9][2] * K2 + A[9][3] * K3 + A[9][4] * K4
               + A[9][5] * K5 + A[9][6] * K6 + A[9][7] * K7 + A[9][8] * K8;
        K9 = h * P.RHS(t + C[9]*h, X);

        //>> Stage 11
        X = Xk + A[10][0] * K0 + A[10][1] * K1 + A[10][2] * K2 + A[10][3] * K3 + A[10][4] * K4
               + A[10][5] * K5 + A[10][6] * K6 + A[10][7] * K7 + A[10][8] * K8 + A[10][9] * K9;
        K10 = h * P.RHS(t + C[10]*h, X);

        //>> Stage 12
        X = Xk + A[11][0]  * K0 + A[11][1] * K1 + A[11][2] * K2 + A[11][3] * K3 + A[11][4] * K4
               + A[11][5]  * K5 + A[11][6] * K6 + A[11][7] * K7 + A[11][8] * K8 + A[11][9] * K9
               + A[11][10] * K10;
        K11 = h * P.RHS(t + C[11]*h, X);

        //>> Stage 13
        X = Xk + A[12][0]  * K0  + A[12][1]  * K1 + A[12][2] * K2 + A[12][3] * K3 + A[12][4] * K4
               + A[12][5]  * K5  + A[12][6]  * K6 + A[12][7] * K7 + A[12][8] * K8 + A[12][9] * K9
               + A[12][10] * K10 + A[12][11] * K11;
        K12 = h * P.RHS(t + C[12]*h, X);



        // Estimate y(x+h) of the highest order available
        Xk = Xk + B[0]  * K0
                + B[1]  * K1
                + B[2]  * K2
                + B[3]  * K3
                + B[4]  * K4
                + B[5]  * K5
                + B[6]  * K6
                + B[7]  * K7
                + B[8]  * K8
                + B[9]  * K9
                + B[10] * K10
                + B[11] * K11
                + B[12] * K12;
        

        progress = 100 * (t0 + h * steps) / (tf-t0); 
        if( steps % 10000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return Xk;
};
