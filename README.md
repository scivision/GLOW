[![Build Status](https://travis-ci.org/scivision/GLOW.svg?branch=master)](https://travis-ci.org/scivision/GLOW)

# GLOW
The GLobal airglOW Model

This directory contains:

   * Modern Fortran source code files,
   * Makefile
   * CMakeLists.txt
   * Example input and output files
   * Example job script
   * `data/` contains input data files,
   * `data/iri90/` contains IRI input data files
   * [Release Notes](ReleaseNotes.rst)


## Python install
This command automatically compiles the Fortran code to access from Python on any platform.

    python setup.py develop

You can then run the self-tests with

    python tests/test.py -v

## Fortran-only mode
While many users use the Python interface to Glow, users on HPCC may want to use MPI directly in Fortran using the Makefiles or CMake. Here's how to compile the ``basic`` example with CMake.

    cd build
    cmake ..
    make
    cd ..

This creates:


executable  |  description
------------|--------------
basic        |   basic GLOW 
driver  | MPI/NetCDF GLOW 

### Fortran examples
With regard to precision, at first try I occasionally see the least digit of precision in the text output files differ. 
This can be related to Stan using Intel Fortran compiler and I'm using Gfortran.
E.g. I get 4.34e-9 and Stan got 4.33e-9.

I also get the message::

    Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
    
Which is a warning sign that different compilers and platforms may yield different results of unpredictable variance.



#### Aurora example (night)
Compare your results vs. Stan's results with `meld`:

    ./basic < in.basic.aur > aur.basic.out

    meld out.basic.aur aur.basic.out


### Aurora Example (night)

    ./basic < in.basic.day > day.basic.out

    meld out.basic.day day.basic.out

### MPI Prereq
Allows parallel execution of GLOW Fortran code for HPCC.

    sudo apt install libopenmpi-dev
    
### NetCDF prereq
allows writing data in NetCDF (optional).

    sudo apt install libnetcdf-dev

### Select Fortran compiler
Simply use the variable `FC`. Example

    FC=ifort make
