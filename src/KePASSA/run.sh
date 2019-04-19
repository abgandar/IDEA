#!/bin/sh
# Run from the top level of the CMake build directory

# generate points
./src/points > points.dat

# feed points to Hodei's program for reference solutions
epsR=0.001        # in km, as in C++ code
epsV=0.0001
scaleX=`bc -l <<< "$epsR * 1000"`   # convert to m
scaleV=`bc -l <<< "$epsV * 1000"`
echo > reference.dat
while read -r X1
do
read -r X2
read -r X3
read -r X4
read -r X5
read -r X6
# Hodei's code is in m, m/s (this is very slow)
#X1=`bc -l <<< "$X1 * $epsR * 1000"`
#X2=`bc -l <<< "$X2 * $epsR * 1000"`
#X3=`bc -l <<< "$X3 * $epsR * 1000"`
#X4=`bc -l <<< "$X4 * $epsV * 1000"`
#X5=`bc -l <<< "$X5 * $epsV * 1000"`
#X6=`bc -l <<< "$X6 * $epsV * 1000"`

#../src/KePASSA/SS2bSol/bin/Propagator $X1 $X2 $X3 $X4 $X5 $X6 $scaleX $scaleV 1e-13 1e-13 248942.32365024 >> reference.dat
../src/KePASSA/SS2bSol/bin/Propagator $X1 $X2 $X3 $X4 $X5 $X6 $scaleX $scaleV 1e-13 1e-13 497884.64730048 >> reference.dat

done < points.dat

