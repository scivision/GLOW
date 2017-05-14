cmake_minimum_required (VERSION 2.8.12)
project(glow Fortran)

# note two "ff" on libnetcdff-dev
#
# apt install libnetcdff-dev libopenmpi-dev

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

include_directories(../build) # for .mod

add_compile_options(-Warray-bounds -fexceptions -g)
#---------- NetCDF (optional) ------------
find_package(netCDF)
if(netCDF_DIR)
    include_directories(${netCDF_INCLUDE_DIR})

    add_library(ionetcdf readtgcm.f90 output.f90 )
    target_link_libraries(ionetcdf netcdff) # yes double "ff"
else()
    message(WARNING "NetCDF failed to be included " ${netCDF_INCLUDE_DIR})
endif()
#---------- MPI (optional) -------------
find_package(MPI)
if(MPI_Fortran_COMPILER)
    # wrapper compiler mpif90
    #set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER}) 
    
    add_definitions(${MPI_Fortran_COMPILE_FLAGS})
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    add_compile_options(${MPI_Fortran_LINK_FLAGS} -pthread)

    add_executable(driver glowdriver.f90 tzgrid.f90)
    target_link_libraries(driver glow ionetcdf ${MPI_Fortran_LIBRARIES})
else()
    message(WARNING "MPI failed to be included " ${MPI_Fortran_LIBRARIES})
endif()
