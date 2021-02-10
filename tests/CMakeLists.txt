
add_executable(test_ncio test_ncio.f90)

target_link_libraries(test_ncio PRIVATE ncio::ncio NetCDF::NetCDF_Fortran)
add_test(NAME test_ncio  COMMAND test_ncio
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

