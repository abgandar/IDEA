#! /bin/bash


# Cleam using Makefile
cd ./build
make clean
cd ..

# Remove CMakeFiles Folder
#rm -r CMakeFiles

# Remove bin Folder
#rm -r bin

# Remove .MOD files
#rm -f *.mod

# Remove data that the software might have written
#rm -f *.dat
rm -f *.plt

# Remove other CMake generated files
#rm -f cmake_install.cmake
#rm -f CMakeCache.txt
#rm -f Makefile

rm -f -r build/*
