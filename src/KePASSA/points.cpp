/*
    Generates a cloud of points to test
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

//static const unsigned int NRANDOM = 1000000;  // number of random points, adjust as needed
static const unsigned int NRANDOM = 10000;  // number of random points, adjust as needed

int main(void)
{
    // set output precision to something useful
    cout.precision(16);
    cout.setf( std::ios::fixed, std:: ios::floatfield );
    // center point
    for(unsigned int i = 0; i < 6; i++)
        cout << 0.0 << endl;
    // iterate all 2^6 corner points
    for(unsigned int j = 0; j < (1<<7); j++)
    {
        unsigned int temp = j;
        for(unsigned int i = 0; i < 6; i++)
        {
            cout << ((temp&1) ? -1.0 : 1.0) << endl;
            temp >>= 1;
        }
    }
    // random points 
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution( -1.0, 1.0 );
    for(unsigned int j = 0; j < NRANDOM; j++)
        for(unsigned int i = 0; i < 6; i++)
            cout << distribution(generator) << endl;

    return 0;
}
