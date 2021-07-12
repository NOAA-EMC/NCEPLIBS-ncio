# This is the CMake file for the tests directory of NCEPLIBS-ncio.

# Syntax:  add_mpi_test (<TESTNAME>
#                        EXECUTABLE <command>
#                        ARGUMENTS <arg1> <arg2> ...
#                        NUMPROCS <num_procs>
#                        TIMEOUT <timeout>)
function (add_mpi_test TESTNAME)

  # Parse the input arguments
  set (options)
  set (oneValueArgs NUMPROCS TIMEOUT EXECUTABLE)
  set (multiValueArgs ARGUMENTS)
  cmake_parse_arguments (${TESTNAME} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Store parsed arguments for convenience
  set (exec_file ${${TESTNAME}_EXECUTABLE})
  set (exec_args ${${TESTNAME}_ARGUMENTS})
  set (num_procs ${${TESTNAME}_NUMPROCS})
  set (timeout ${${TESTNAME}_TIMEOUT})

  if(MPIEXEC_EXECUTABLE)
    set( MPIEXEC ${MPIEXEC_EXECUTABLE} )
  endif()
  if(NOT MPIEXEC)
    find_program( MPIEXEC NAMES mpiexec mpirun lamexec srun
                  DOC "Executable for running MPI programs." )
  endif()

  # Run tests directly from the command line
  set(EXE_CMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${num_procs}
      ${MPIEXEC_PREFLAGS}  ${exec_file}
      ${MPIEXEC_POSTFLAGS} ${exec_args})

  # Add the test to CTest
  add_test(NAME ${TESTNAME} COMMAND ${EXE_CMD})

  # Adjust the test timeout
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${timeout})

endfunction()

# Copy some test files from current source dir to out-of-tree build dir.
execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/dynf000_template.nc.in
    ${CMAKE_CURRENT_BINARY_DIR}/dynf000_template.nc.in)

execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/dynf000_template.nc.in
    ${CMAKE_CURRENT_BINARY_DIR}/dynf000_par_template.nc.in)

# Build and run test.
add_executable(tst_ncio tst_ncio.F90)
target_link_libraries(tst_ncio PRIVATE ncio::ncio)
add_test(NAME tst_ncio  COMMAND tst_ncio)

# Brian Curtis: MPI Test
add_executable(tst_ncio_mpi tst_ncio_mpi.F90)
target_link_libraries(tst_ncio_mpi PRIVATE ncio::ncio)
add_mpi_test(tst_ncio_mpi
  EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/tst_ncio_mpi
  NUMPROCS 4
  TIMEOUT 120)
