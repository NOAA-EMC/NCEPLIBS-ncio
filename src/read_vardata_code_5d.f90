  !> @file
  !! @brief Read 5D var data.
  !! optional return errcode
  !! @author Jeff Whitaker, Cory Martin
  type(Dataset), intent(in) :: dset
  character(len=*), intent(in) :: varname
  integer, allocatable, dimension(:) :: start, count
  integer, intent(out), optional :: errcode
  integer, intent(in), optional :: nslice
  integer, intent(in), optional :: slicedim
  integer, intent(in), optional :: ncstart(5)
  integer, intent(in), optional :: nccount(5)
  integer ncerr, nvar, n1,n2,n3,n4,n5,nd,ncount
  logical return_errcode
  ! check if use the errcode
  if(present(errcode)) then
     return_errcode=.true.
     errcode = 0
  else
     return_errcode=.false.
  endif
  if (present(nslice)) then
     ncount = nslice
  else
     ncount = 1
  endif
  nvar = get_nvar(dset,varname)
  if (dset%variables(nvar)%ndims /= 5) then
     if (return_errcode) then
        errcode=nf90_ebaddim
        return
     else
        print *,'rank of data array != variable ndims (or ndims-1)'
        stop 99
     endif
  endif
  if (present(slicedim)) then
     nd = slicedim
  else
     nd = dset%variables(nvar)%ndims
  end if
  n1 = dset%variables(nvar)%dimlens(1)
  n2 = dset%variables(nvar)%dimlens(2)
  n3 = dset%variables(nvar)%dimlens(3)
  n4 = dset%variables(nvar)%dimlens(4)
  n5 = dset%variables(nvar)%dimlens(5)
  ! allocate/deallocate values
  if (allocated(values)) deallocate(values)
  allocate(start(dset%variables(nvar)%ndims),count(dset%variables(nvar)%ndims))
  start(:) = 1; count(1)=n1; count(2)=n2; count(3)=n3; count(4)=n4; count(5)=n5
  if (present(ncstart) .and. present(nccount)) then
     start(1)=ncstart(1); count(1)=nccount(1)
     start(2)=ncstart(2); count(2)=nccount(2)
     start(3)=ncstart(3); count(3)=nccount(3)
     start(4)=ncstart(4); count(4)=nccount(4)
     start(5)=ncstart(5); count(5)=nccount(5)
  else if (present(nslice)) then
     start(nd)=ncount; count(nd)=1
  endif
  allocate(values(count(1),count(2),count(3),count(4),count(5)))
  ncerr = nf90_get_var(dset%ncid, dset%variables(nvar)%varid, values,&
                       start,count)
  ! err check
  if (return_errcode) then
     call nccheck(ncerr,halt=.false.)
     errcode=ncerr
  else
     call nccheck(ncerr)
  endif
