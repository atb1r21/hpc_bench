# DL_MG

DL_MG is a parallel (MPI+OpenMP) 3D geometric high order finite difference
multigrid solver for Poisson and Poisson-Boltzmann Equations written
in Fortran.

License: BSD 3-Clause License
Copyright (c) 2016, Chris-Kriton Skylaris and the Numerical Algorithms Group Ltd.
All rights reserved.

Authors: Lucian Anton, James Womack and Jacek Dziedzic

Details of the algorithm and subroutines documentations are available
in docs/dl_mg.pdf or for the last release on the web at https://antonl22.bitbucket.io/dlmg

Compilation:

Several make variable, mainly describing the compiler and its flags, need to be initailise for the compilation.
These variables must be provided in a make include file located in platforms/<platform-name>.inc

A good starting point for a template is the file platforms/parallel_laptop.inc which
defines the variables for a local workstation with gcc and MPI installed.

The file platforms/archer.inc is for a Cray XE6 system which uses
several development environments accessible with module command.  Some
other platforms files are includes in the folder 'platforms'.

After the platform file is created the solver library can be build with

make PLATFORM=<platform-name>

Link the code with libdlmg.a and include 'use dl_mg' in the calling subroutine.
