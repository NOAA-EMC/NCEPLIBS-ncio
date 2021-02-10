
add_executable(test_ncio test_ncio.f90)

target_link_libraries(test_ncio PRIVATE ncio::ncio NetCDF::NetCDF_Fortran)
add_test(test_ncio test_ncio)

