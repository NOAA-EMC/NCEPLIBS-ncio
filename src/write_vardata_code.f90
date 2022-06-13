  !> @file
  !! @brief Write variable data.
  !! Optional return errcode
  !! @author Jeff Whitaker, Cory Martin
  type(Dataset), intent(inout) :: dset
  character(len=*), intent(in) :: varname
  integer, intent(in), optional :: nslice
  integer, intent(in), optional :: slicedim
  integer, intent(out), optional :: errcode
  integer ncerr, nvar, ncount, nd, n, ndim
  integer, allocatable, dimension(:) :: start, count, varshape
  logical is_slice
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
     is_slice = .true.
  else
     ncount = 1
     is_slice = .false.
  endif
  ! define variable name and allocate variable
  nvar = get_nvar(dset,varname)
  allocate(start(dset%variables(nvar)%ndims),count(dset%variables(nvar)%ndims))
  allocate(varshape(dset%variables(nvar)%ndims))
  start(:) = 1
  count(:) = 1
  if (is_slice) then
     nd = slicedim
  else
     nd = dset%variables(nvar)%ndims
  end if
  do n=1,dset%variables(nvar)%ndims
     ndim = dset%variables(nvar)%dimids(n)
     if (is_slice .and. n == nd) then
        start(n) = ncount
        count(n) = 1
     else if (n == nd .and. dset%dimensions(ndim)%isunlimited) then
        start(n) = ncount
        varshape = shape(values)
        count(n) = varshape(n)
     else
        start(n) = 1
        count(n) = dset%variables(nvar)%dimlens(n)
     end if
  end do
  ! write operations on a parallel file system are performed collectively
  ncerr = nf90_var_par_access(dset%ncid, dset%variables(nvar)%varid, nf90_collective)
  if (is_slice) then
     if (dset%variables(nvar)%ndims > 1) then
        ncerr = nf90_put_var(dset%ncid, dset%variables(nvar)%varid,values, &
                             start=start,count=count)
     else if (dset%variables(nvar)%ndims == 1) then
        if (return_errcode) then
           errcode = -1
           return
        else
           print *,'cannot write a slice to a 1d variable'
           stop 99
        endif
     endif
  else if (present(ncstart) .and. present(nccount)) then
     ncerr = nf90_put_var(dset%ncid, dset%variables(nvar)%varid,values, &
          start=ncstart, count=nccount)
  else
     ncerr = nf90_put_var(dset%ncid, dset%variables(nvar)%varid, values, &
          start=start, count=count)
  endif
  if (return_errcode) then
     call nccheck(ncerr,halt=.false.)
     errcode=ncerr
     if (ncerr /= 0) return
  else
     call nccheck(ncerr)
  endif
  ! reset unlim dim size for all variables
  if (dset%variables(nvar)%hasunlim) then
     if (return_errcode) then
        call set_varunlimdimlens_(dset,errcode)
     else
        call set_varunlimdimlens_(dset)
     endif
  endif
