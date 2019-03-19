template <class T>
AlgebraicVector<T> Euler(double t0, double tf, double h, const AlgebraicVector<T> &X0, Propagator<T> &P)
{
    AlgebraicVector<T> X = X0;
    int steps = (tf - t0) / h;
    double progress;

    for (int i = 0; i < steps; i += 1)
    {
        X = X + h * P.RHS(t0 + h * i, X);

        progress = 100 * (t0 + h * i) / (tf-t0); 
        if( i % 10000 == 0) std::cout << "  Steps: " << i << "\t" << progress << " \%" << std::endl;
    }

    return X;
};
