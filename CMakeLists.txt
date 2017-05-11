cmake_minimum_required (VERSION 2.8.12)
project(glow Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

add_compile_options(-Warray-bounds -fexceptions)
#-------------------------- non - MPI -------
add_executable(basic glowbasic.f90)
target_link_libraries(basic glow)
#---------- NetCDF (optional) ------------
find_package(netCDF)
if(netCDF_DIR)
    list(APPEND CMAKE_MODULE_PATH ${netCDF_DIR})
#    add_definitions(${netCDF_C_COMPILER_FLAGS})
#    message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
endif()
add_library(tcgm readtgcm.f90)
target_link_libraries(tcgm ${netCDF_LIBRARIES})
#---------- MPI (optional) -------------
find_package(MPI)
if(MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER}) #wrapper compiler mpif90
endif()

# --------- GLOW library -------------------
add_library(glow cglow.f90  glow.f90 bands.f90 conduct.f90 egrid.f90 ephoto.f90 etrans.f90 exsect.f fieldm.f gchem.f90 geomag.f90 maxt.f90 mzgrid.f90 qback.f90 rcolum.f90 rout.f90 snoem.f90 snoemint.f90 solzen.f90 ssflux.f90 iri90.f nrlmsise00.f)



# message("netCDF_LIBRARIES = ${netCDF_LIBRARIES}")  # simply 'netcdf'
add_executable(driver glowdriver.f90)
target_link_libraries(driver glow tcgm )

