# This is the CMake file for the tests directory of NCEPLIBS-ncio.

# Copy some test files from current source dir to out-of-tree build dir.
FILE(GLOB COPY_FILES ${CMAKE_CURRENT_SOURCE_DIR}/dynf000_template.nc.in)
FILE(COPY ${COPY_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Build and run test.
add_executable(tst_ncio tst_ncio.F90)
target_link_libraries(tst_ncio PRIVATE ncio::ncio NetCDF::NetCDF_Fortran)
add_test(NAME tst_ncio  COMMAND tst_ncio)


