# Z2 gauge theory in 4d

## Julia

An improved version of Julia implementation based on the discussion [here](
https://discourse.julialang.org/t/accessing-array-element-splatting-vs-explicit-writing-them-out/64973
).

This is a basic Python (+ Julia) adaptation of M. Creutz's Z2 code in C.

## Keep It Simple and Stupid

This code is for learning the basics of Python and Lattice Gauge Theory.

The python code is very slow. Definitely not suitable for production.
(There are many ways to speed things up. We don't do it here.)

The Julia version does a reasonable job. (now even better! 07.2021)

## How to Run?

Simply

> gcc z2.c; ./a.out

> g++ z2.cpp; ./a.out

> gfortran z2.f90; ./a.out

> python3 z2.py

> julia z2.jl

## Comments are welcomed.
