|Build Status|

GLOW
====

The GLobal airglOW Model

This directory contains:

-  Modern Fortran source code files,
-  Makefile
-  CMakeLists.txt
-  Example input and output files
-  Example job script
-  ``data/`` contains input data files,
-  ``data/iri90/`` contains IRI input data files
-  `Release Notes <ReleaseNotes.rst>`__

Python install
--------------

This command automatically compiles the Fortran code to access from
Python on any platform::

    pip install -e .

Run self-test with::

    pytest -v

Fortran-only (optional)
-----------------------

While many users use the Python interface to Glow, users on HPCC may
want to use MPI directly in Fortran using the Makefiles or CMake. 
Here's how to compile the ``basic`` example with CMake::

    cd bin
    cmake ..
    make
    cd ..

This creates:

+--------------+-------------------+
| executable   | description       |
+==============+===================+
| basic        | basic GLOW        |
+--------------+-------------------+
| driver       | MPI/NetCDF GLOW   |
+--------------+-------------------+
| testdrv      | Python intfc test |
+--------------+-------------------+

Fortran examples
~~~~~~~~~~~~~~~~

With regard to precision, at first try I occasionally see the least
digit of precision in the text output files differ. This can be related
to Stan using Intel Fortran compiler and I'm using Gfortran. E.g. I get
4.34e-9 and Stan got 4.33e-9.

I also get the message::

::

    Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL

Which is a warning sign that different compilers and platforms may yield
different results of unpredictable variance.

Aurora example (night)
^^^^^^^^^^^^^^^^^^^^^^

Compare your results vs. Stan's results:: with ``meld``::

    ./basic < in.basic.aur > aur.basic.out

    meld reference/aur981.basic.out <(./basic < in.basic.aur)

Dayglow Example (day)
^^^^^^^^^^^^^^^^^^^^^^

::

    ./basic < in.basic.day > day.basic.out


To compare::

    meld reference/out.basic.day <(./basic < in.basic.day)

MPI Prereq
~~~~~~~~~~

Allows parallel execution of GLOW Fortran code for HPC::

    sudo apt install libopenmpi-dev

NetCDF prereq
~~~~~~~~~~~~~

allows reading/writing data in NetCDF (optional)::

    sudo apt install libnetcdf-dev libnetcdff-dev

Select Fortran compiler
~~~~~~~~~~~~~~~~~~~~~~~

Simply use the variable ``FC``::

    cd bin/
    rm -r *

    FC=ifort cmake ..
    make

NOTE: Using the Intel compiler requires that you have `built NetCDF
using Intel
Fortran <https://software.intel.com/en-us/articles/performance-tools-for-software-developers-building-netcdf-with-the-intel-compilers/>`__--this
is an issue ANYTIME you use the Intel Compiler and NetCDF.

.. |Build Status| image:: https://travis-ci.org/scivision/GLOW.svg?branch=master
   :target: https://travis-ci.org/scivision/GLOW
