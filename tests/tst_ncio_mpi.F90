! This is a test for the NCEPLIBS-ncio package.
! Author: Brian Curtis (brian.curtis@noaa.gov) May 2021
program tst_ncio_mpi

  use mpi
  use netcdf
  use module_ncio
  implicit none

  interface
    subroutine check(errorcode)
      integer, intent(in) :: errorcode
    end subroutine check
  end interface

  type(Dataset) :: dset, dsetin, dset_test
  integer, parameter :: YT = 128, XT = 256
  integer, parameter :: HALF_YT = YT/2, HALF_XT = XT/2
  integer, parameter :: MAXDIM = 3
  real(4), allocatable, dimension(:,:,:) :: values_3d, pressfc_check
  ! integer ndim,nvar,ndims,ival,idate(6),icheck(6),ierr,n,nbits
  integer :: my_rank, nprocs
  integer :: mpi_err
  integer :: errcode
  integer :: start(MAXDIM), count(MAXDIM)

  call mpi_init(mpi_err)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, mpi_err)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, mpi_err)

  if (my_rank .eq. 0) print *, '*** Testing NCEPLIBS-ncio with MPI.'

  if (nprocs .ne. 4) then
    print *, 'This test must be run using only 4 processors.'
    stop 4
  endif

  if (my_rank .eq. 0) print *, '*** Testing function open_dataset with paropen=.true.'
  dsetin = open_dataset('dynf000_template.nc.in', errcode=errcode, paropen=.true.)
  call check(errcode)

  if (my_rank .eq. 0) print *,'*** Test creation of new dataset from template...'
  dset = create_dataset('dynf000_par.nc', dsetin, paropen=.true., errcode=errcode)
  call check(errcode)

  if (my_rank .eq. 0) print *,'*** Test that number of variables,dimensions,attributes is read...'

  if (dsetin%nvars .ne. 25) stop 5
  if (dsetin%ndims .ne. 7) stop 6
  if (dsetin%natts .ne. 8) stop 7


  if (my_rank .eq. 0) print *,'*** Test MPI read of 3-D pressfc variable...'

  ! Let's set the starting point of data reads for each pe
  count = (/ HALF_XT, HALF_YT, 1 /)
  if (my_rank .eq. 0) then
     start = (/ 1, 1, 1 /)
  else if (my_rank .eq. 1) then
     start = (/ 1, HALF_YT + 1, 1 /)
  else if (my_rank .eq. 2) then
     start = (/ HALF_XT + 1, 1, 1 /)
  else if (my_rank .eq. 3) then
     start = (/ HALF_XT + 1, HALF_YT + 1, 1 /)
  endif

  call read_vardata(dsetin, 'pressfc', values_3d, ncstart=start, nccount=count, errcode=errcode)
  call check(errcode)

  call close_dataset(dsetin, errcode=errcode)
  call check(errcode)

  if (my_rank .eq. 0) print *,'*** Test write of variable data...'
  call write_vardata(dset,'pressfc', values_3d, ncstart=start, nccount=count, errcode=errcode)
  call check(errcode)

  if (allocated(values_3d)) deallocate(values_3d)

  call close_dataset(dset, errcode=errcode)
  call check(errcode)

  if (my_rank .eq. 0) print *,'*** Reading in new data set and verifying max value'
  dset_test = open_dataset('dynf000_par.nc', errcode=errcode)
  call check(errcode)

  call read_vardata(dset_test, 'pressfc', pressfc_check, errcode=errcode)
  call check(errcode)

  if (maxval(pressfc_check) .ne. 102345.6) stop 33
  call close_dataset(dset_test, errcode=errcode)
  call check(errcode)

  call mpi_finalize(mpi_err)
  if (my_rank .eq. 0) print*, "*** SUCCESS!"

end program tst_ncio_mpi

subroutine check(errorcode)
  use netcdf
  implicit none
  integer, intent(in) :: errorcode

  if(errorcode /= nf90_noerr) then
     print *, 'Error: ', trim(nf90_strerror(errorcode))
     stop 99
  endif
end subroutine check
