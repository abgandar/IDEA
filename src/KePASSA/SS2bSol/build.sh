#! /bin/bash

clear


# Get into build folder
cd build


# Select compilers bundle or toolchain. Options: intel, GNU
Bundle=GNU


# Launch CMake
if [ "$Bundle" = "intel" ]; then

	#FC=ifort CC=icc CXX=icpc cmake .
	FC=ifort cmake ..

elif [ "$Bundle" = "GNU" ]; then

	#FC=gfortran CC=gcc CXX=g++ cmake .
	FC=gfortran cmake ..

else
	cmake ..
fi


# Launch Make
make


# Go up a level to execute the binary
cd ..


# Execute compiled binary, run 3 example propagations and store results in 'Results.plt'
#
#   Input arguments are as follows:
#   X       Initial Cartesian position offset along X axis [m]
#   Y       Initial Cartesian position offset along Y axis [m]
#   Z       Initial Cartesian position offset along Z axis [m]
#   Vx      Initial Cartesian velocity offset along X axis [m/s]
#   Vy      Initial Cartesian velocity offset along Z axis [m/s]
#   Vz      Initial Cartesian velocity offset along Y axis [m/s]
#   RelTol  Relative integration tolerance                 [-]
#   AbsTol  Absolute integration tolerance                 [-]
#   Tf      Final time to stop integration                 [s]
#
#   Output is a one-line string containing 6 numbers, respectively:
#   X       Final Cartesian position projected on X axis [m]
#   Y       Final Cartesian position projected on Y axis [m]
#   Z       Final Cartesian position projected on Z axis [m]
#   Vx      Final Cartesian velocity projected on X axis [m/s]
#   Vy      Final Cartesian velocity projected on Y axis [m/s]
#   VZ      Final Cartesian velocity projected on Z axis [m/s]
./bin/Propagator 1e-6 1e-6 1e-6 1e-9 1e-9 1e-9 1e-20 1e-20       >  Results.plt
./bin/Propagator 1e-6 1e-6 1e-6 1e-9 1e-9 1e-9 1e-22 1e-22       >> Results.plt
./bin/Propagator 1e-6 1e-6 1e-6 1e-9 1e-9 1e-9 1e-22 1e-22 10000 >> Results.plt
