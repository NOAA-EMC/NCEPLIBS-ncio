# This is the CMake file for the tests directory of NCEPLIBS-ncio.

add_executable(tst_ncio tst_ncio.F90)

target_link_libraries(tst_ncio PRIVATE ncio::ncio NetCDF::NetCDF_Fortran)
add_test(NAME tst_ncio  COMMAND tst_ncio
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

