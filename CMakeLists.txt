cmake_minimum_required (VERSION 2.8.12)
project(glow Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ..)

add_compile_options(-O3 -g)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    add_compile_options(-debug full -traceback -fpe-all=0 -check all)
else()
    add_compile_options(-mtune=native -Warray-bounds -fexceptions)
endif()
# --------- GLOW library -------------------
add_library(glow cglow.f90 glow.f90 bands.f90 conduct.f90 egrid.f90 ephoto.f90 etrans.f90 exsect.f fieldm.f gchem.f90 geomag.f90 maxt.f90 mzgrid.f90 qback.f90 rcolum.f90 rout.f90 snoem.f90 snoemint.f90 solzen.f90 ssflux.f90 iri90.f nrlmsise00.f)
#------------- basic example ------------
add_executable(basic glowbasic.f90)
target_link_libraries(basic glow)
#include_directories(../build) # for .mod
#---------- NetCDF (optional) ------------
find_package(netCDF)
if(netCDF_DIR)
    include_directories(${netCDF_INCLUDE_DIR})

    add_library(ionetcdf readtgcm.f90 output.f90 )
    target_link_libraries(ionetcdf netcdff glow) # yes double "ff"
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
    add_compile_options(-pthread)

    add_executable(driver glowdriver.f90 tzgrid.f90)
    target_link_libraries(driver glow ionetcdf ${MPI_Fortran_LIBRARIES})
else()
    message(WARNING "MPI failed to be included " ${MPI_Fortran_LIBRARIES})
endif()

