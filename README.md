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

### MPI Prereq
Allows parallel execution of GLOW Fortran code for HPCC.

    sudo apt install libopenmpi-dev
    
### NetCDF prereq
allows writing data in NetCDF (optional).

    sudo apt install libnetcdf-dev

### Select Fortran compiler
Simply use the variable `FC`. Example

    FC=ifort make
