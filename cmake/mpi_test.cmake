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

# Check whether OpenMPI, in which case we need --oversubscribe in GitHub CI
  execute_process(COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE VERSIONOUTPUT)
  if(VERSIONOUTPUT MATCHES ".*OpenRTE.*")
    set(MPIEXEC_PREFLAGS "--oversubscribe" ${MPIEXEC_PREFLAGS})
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
