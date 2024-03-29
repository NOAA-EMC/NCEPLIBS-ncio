  !> @file
  !! @brief Read 3D var data.
  !! optional return errcode
  !! @author Jeff Whitaker, Cory Martin
  type(Dataset), intent(in) :: dset
  character(len=*), intent(in) :: varname
  integer, intent(in), optional :: nslice
  integer, intent(in), optional :: slicedim
  integer, intent(in), optional :: ncstart(3)
  integer, intent(in), optional :: nccount(3)
  integer, intent(out), optional :: errcode
  integer ncerr, nvar, n, nd, ndim, ncount
  integer, allocatable, dimension(:) :: start, count
  integer :: dimlens(3)
  logical return_errcode
  ! check if use the errcode
  if(present(errcode)) then
     return_errcode=.true.
     errcode = 0
  else
     return_errcode=.false.
  endif
  ! check if count/dims of data are avail
  if (present(nslice)) then
     ncount = nslice
  else
     ncount = 1
  endif
  nvar = get_nvar(dset,varname)
  allocate(start(dset%variables(nvar)%ndims),count(dset%variables(nvar)%ndims))
  start(:) = 1
  count(:) = 1
  dimlens(:) = 1
  if (present(slicedim)) then
     nd = slicedim
  else
     ! slicedim not specified, if data array one dim smaller than
     ! variable slice along last dimension of variable
     nd = dset%variables(nvar)%ndims
  end if
  ndim = 1
  do n=1,dset%variables(nvar)%ndims
     if (present(slicedim) .and. n == nd) then
        start(n) = ncount
        count(n) = 1
     else if (.not. present(slicedim) .and. n == nd .and. dset%variables(nvar)%ndims == 4) then
        start(n) = ncount
        count(n) = 1
     else
        start(n) = 1
        count(n) = dset%variables(nvar)%dimlens(n)
        dimlens(ndim) = dset%variables(nvar)%dimlens(n)
        ndim = ndim + 1
     end if
  end do

  if (dset%variables(nvar)%ndims /= 3 .and. dset%variables(nvar)%ndims /= 4) then
     if (return_errcode) then
        call nccheck(ncerr,halt=.false.)
        errcode=nf90_ebaddim
        return
     else
        print *,'rank of data array != variable ndims (or ndims-1)'
        stop 99
     endif
  endif
  ! allocate/deallocate values
  if (allocated(values)) deallocate(values)
  if (present(ncstart) .and. present(nccount)) then
     allocate(values(nccount(1),nccount(2),nccount(3)))
     start(1)=ncstart(1); count(1)=nccount(1)
     start(2)=ncstart(2); count(2)=nccount(2)
     start(3)=ncstart(3); count(3)=nccount(3)
     if (dset%variables(nvar)%ndims == 4) then
        start(4)=1; count(4)=1
     end if
  else
     if (dset%variables(nvar)%ndims == 4) then
        allocate(values(dimlens(1),dimlens(2),dimlens(3)))
     else
        allocate(values(dset%variables(nvar)%dimlens(1),&
             dset%variables(nvar)%dimlens(2),&
             dset%variables(nvar)%dimlens(3)))
     end if
  end if
  ncerr = nf90_get_var(dset%ncid, dset%variables(nvar)%varid, values,&
                       start=start, count=count)
  ! err check
  if (return_errcode) then
     call nccheck(ncerr,halt=.false.)
     errcode=ncerr
  else
     call nccheck(ncerr)
  endif
