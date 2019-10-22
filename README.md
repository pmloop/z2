# Z2 gauge theory in 4d

This is a basic Python adaptation of M. Creutz's Z2 code in C.

## Keep It Simple and Stupid

This code is for learning the basics of Python and Lattice Gauge Theory.

The python code is very slow. Definitely not suitable for production.
There are many ways to speed things up. 
For example, writing the Monte Carlo update routine and measurement routines in c / cpp / fortran 
and run it from python via F2PY. We do not do it here.

## How to Run?

Simply

> gcc z2.c; ./a.out

> g++ z2.cpp; ./a.out

> gfortran z2.f90; ./a.out

> python3 z2.py


## Comments are welcomed.
