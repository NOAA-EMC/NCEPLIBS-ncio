# This is the CMake file for the tests directory of NCEPLIBS-ncio.

# Copy some test files from current source dir to out-of-tree build dir.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/dynf000_template.nc.in
    ${CMAKE_CURRENT_BINARY_DIR}/dynf000_template.nc.in)

# Build and run test.
add_executable(tst_ncio tst_ncio.F90)
target_link_libraries(tst_ncio PRIVATE ncio::ncio NetCDF::NetCDF_Fortran)
add_test(NAME tst_ncio  COMMAND tst_ncio)


