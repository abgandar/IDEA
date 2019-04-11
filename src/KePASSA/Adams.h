#include <dace/dace.h>
using namespace std;
using namespace DACE;

template <class T, unsigned int ORDER, unsigned int ITERS>
/*  ORDER: Order of the Adams method. Coincides with # of steps for Adams-Bashforth, or # of steps + 1 for Adams-Moulton.
    ITERS: Number of 'corrector' iterations to be used, i.e.
            - If 0, the method reduces to an explicit Adams-Bashforth method.
            - If 1, the method becomes a PECE     Adams-Bashforth-Moulton method
            - If 2, the method becomes a PECECE   Adams-Bashforth-Moulton method
            - If n, the method becomes a PE(CE)^n Adams-Bashforth-Moulton method */
AlgebraicVector<T> Adams(double t0, double t1, double h, double tol, AlgebraicVector<T> x0, AlgebraicVector<T> (*RHS)(AlgebraicVector<T>, double), vector<AlgebraicVector<T>> &res, bool save)
{

    // '\beta_{i,j}' Coefficients of Adams-Bashforth methods up to order 12
    const double B_Bashforth[12][12] = {        // Note: C++ initializes all omitted values to 0.0
        { 1. },
        { -0.5, 3./2. },
        { 5./12., -16./12., 23./12. },
        { -9./24., 37./24., -59./24., 55./24. },
        { 251./720., -1274./720., 2616./720., -2774./720., 1901./720. },
        { -475./1440., 2877./1440., -7298./1440., 9982./1440., -7923./1440., 4277./1440. },
        { 19087./60480., -134472./60480., 407139./60480., -688256./60480., 705549./60480., -447288./60480., 198721./60480. },
        { -36799./120960., 295767./120960., -1041723./120960., 2102243./120960., -2664477./120960., 2183877./120960., -1152169./120960., 434241./120960. },
        { 1070017./3628800., -4832053./1814400., 19416743./1814400., -45586321./1814400., 862303./22680., -69927631./1814400., 47738393./1814400., -21562603./1814400., 14097247./3628800. },
        { -25713./89600., 20884811./7257600., -2357683./181440., 15788639./453600., -222386081./3628800., 269181919./3628800., -28416361./453600., 6648317./181440., -104995189./7257600., 4325321./1036800. },
        { 26842253./95800320., -52841941./17107200., 2472634817./159667200., -186080291./3991680., 2492064913./26611200., -82260679./623700., 3539798831./26611200., -1921376209./19958400., 1572737587./31933440., -2067948781./119750400., 2132509567./479001600. },
        { -4777223./17418240., 30082309./9123840., -17410248271./958003200., 923636629./15206400., -625551749./4561920., 35183928883./159667200., -41290273229./159667200., 35689892561./159667200., -15064372973./106444800., 12326645437./191600640., -6477936721./319334400., 4527766399./958003200. },
    };

    // '\beta_{i,j}' Coefficients of Adams-Moulton methods up to order 12
    const double B_Moulton[12][12] = {        // Note: C++ initializes all omitted values to 0.0
        { },
        { 0.5, 0.5 },
        { -1./12., 8./12., 5./12. },
        { 1./24., -5./24., 19./24., 9./24. },
        { -19./720., 106./720., -264./720., 646./720., 251./720. },
        { 27./1440., -173./1440., 482./1440., -798./1440., 1427./1440., 475./1440. },
        { -863./60480., 6312./60480., -20211./60480., 37504./60480., -46461./60480., 65112./60480., 19087./60480. },
        { 1375./120960., -11351./120960., 41499./120960., -88547./120960., 123133./120960., -121797./120960., 139849./120960., 36799./120960. },
        { -33953./3628800., 156437./1814400., -645607./1814400., 1573169./1814400., -31457./22680., 2797679./1814400., -2302297./1814400., 2233547./1814400., 1070017./3628800. },
        { 8183./1036800., -116687./1451520., 335983./907200., -462127./453600., 6755041./3628800., -8641823./3628800., 200029./90720., -1408913./907200., 9449717./7257600., 25713./89600. },
        { -3250433./479001600., 9071219./119750400., -12318413./31933440., 23643791./19958400., -21677723./8870400., 2227571./623700., -33765029./8870400., 12051709./3991680., -296725183./159667200., 164046413./119750400., 26842253./95800320. },
        { 4671./788480., -68928781./958003200., 384709327./958003200., -87064741./63866880., 501289903./159667200., -91910491./17740800., 1007253581./159667200., -102212233./17740800., 36465037./9123840., -99642413./45619200., 1374799219./958003200., 4777223./17418240. },
    };
    
    AlgebraicVector<T> f[ORDER];    // Values of the function 'f' evaluated at previous steps
    AlgebraicVector<T> fn;          // Value  of the function 'f' at the newly computed step
    AlgebraicVector<T> Xp;          // The integration solution at current step during corrector iterations
    AlgebraicVector<T> X = x0;      // The integration solution at the current step
    
    unsigned int steps = 0;         // Counter for the taken integration steps
    double t = t0;                  // Current value of the independent integration variable (usually time)


    // Remember initial conditions as part of the results output
    if( save )
        res.push_back( X );
    


    // Initialization: take the first m-1 steps with another method
    if (ORDER > 1){
        
        while( t < t1 )
        {
            steps++;
            const double dt = min(t1-t, h);

            // Evaluate RHS and store value
            f[ORDER-1-steps] = RHS(X, t);

            // Estimate y(x+dt) using another method of a compatible order (for now use RK8)
            X = RK8 (t, t+dt, dt, tol, X, *RHS, res, false);

            // Update value of independent variables
            t += dt;

            // save step
            if( save )
                res.push_back( X );
            
            // Break loop when initialization is complete
            if (steps == ORDER - 1)
                break;
        }

    }


    // Evaluate RHS at the newly (just) computed integration step
    fn = RHS(X, t);


    // Adams-Bashforth Method with fixed stepsize
    while( t < t1 )
    {
        steps++;
        const double dt = min(t1-t, h);

        // Update (i.e. cycle) the array of function evaluations at previous steps
        for( unsigned int i = ORDER-1; i > 0; --i )
            f[i] = f[i-1];

        f[0] = fn;

        // Estimate y(x+dt) using an m-step Adams-Bashforth-Moulton method (except if the last stepsize needs to be smaller)
        if (dt < h){

            X = RK8 (t, t+dt, dt, tol, X, *RHS, res, false);
        
        }
        else{
        
            // Copy solution at current step
            Xp = X;

            // PREDICT with Adams-Bashforth method
            for( unsigned int j = 0; j < ORDER; j++)
            {
                X += dt * B_Bashforth[ORDER-1][j] * f[ORDER-1-j];
            }

            // EVALUATE RHS at the newly (just) computed integration step
            fn = RHS(X, t+dt);

            // CORRECT (iteratively) with Adams-Moulton method
            for( unsigned int k = 0; k < ITERS; k++){
                
                X = Xp + dt * B_Moulton[ORDER-1][ORDER-1] * fn;

                for( unsigned int j = 0; j < ORDER-1; j++)
                {
                    X += dt * B_Moulton[ORDER-1][j] * f[ORDER-2-j];
                }

                // Re-EVALUATE RHS at the newly (just) computed integration step
                fn = RHS(X, t+dt);
            }
        
        }

        // Update value of independent variables
        t += dt;

        // save step
        if( save )
            res.push_back( X );

        // Display integration progress
        const double progress = 100.00 * (t) / (t1-t0);
        if( steps % 100000 == 0) std::cout << "  Steps: " << steps << "\t" << progress << " \%" << std::endl;
    }

    return X;
};
