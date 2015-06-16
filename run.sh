#!/bin/bash
# Our run code. Should probably have a bunch of different things here that we comment out as needed.

# NMS/LAPACK compilation, this is my prefered code at the moment
gfortran -c GFMod.f95
gfortran main.f95 GFMod.o -O -L$HOME/lib -llapack -lslatec

# Including MKL, the fucking mess that it is
# gfortran -c GFMod.f95
# gfortran main.f95 GFMod.o -O -L$HOME/lib/NMS -lnms -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -fopenmp 

# Test compilation
# gfortran -c GFMod.f95
# gfortran test.f95 GFMod.o -O -L$HOME/lib -lnms -llapack -lslatec

