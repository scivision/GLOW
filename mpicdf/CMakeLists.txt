cmake_minimum_required (VERSION 2.8.12)
project(glow Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

include_directories(../build) # for .mod

add_compile_options(-Warray-bounds -fexceptions)
#---------- NetCDF (optional) ------------
find_package(netCDF)
if(netCDF_DIR)
    list(APPEND CMAKE_MODULE_PATH ${netCDF_DIR})

    include_directories(${netCDF_INCLUDE_DIR})

    add_library(ionetcdf readtgcm.f90 output.f90 )
    target_link_libraries(ionetcdf netcdff) # yes double "ff"
endif()
#---------- MPI (optional) -------------
find_package(MPI)
if(MPI_Fortran_COMPILER)
    # wrapper compiler mpif90
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER}) 

    add_executable(driver glowdriver.f90 tzgrid.f90)
    target_link_libraries(driver glow ionetcdf)
endif()
