cmake_minimum_required (VERSION 2.8.12)
project(glow Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

add_compile_options(-Warray-bounds -fexceptions)
#---------- MPI (optional) -------------
find_package(MPI)
if(MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER}) #wrapper compiler mpif90

    add_executable(basicmpi glowbasic.f90)
    target_link_libraries(basicmpi glow)
endif()
# --------- GLOW library -------------------
add_library(glow cglow.f90  glow.f90 bands.f90 conduct.f90 egrid.f90 ephoto.f90 etrans.f90 exsect.f fieldm.f gchem.f90 geomag.f90 maxt.f90 mzgrid.f90 qback.f90 rcolum.f90 rout.f90 snoem.f90 snoemint.f90 solzen.f90 ssflux.f90 iri90.f nrlmsise00.f)

#------ test program: Auroral emissions ----------
add_executable(basic glowbasic.f90)
target_link_libraries(basic glow)

