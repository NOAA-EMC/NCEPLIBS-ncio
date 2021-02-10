program test_ncio
  use netcdf
  use module_ncio
  implicit none
  
  character(len=72) charatt, time_units
  type(Dataset) :: dset, dsetin
  type(Variable) :: var
  real(4), allocatable, dimension(:) :: values_1d
  real(4), allocatable, dimension(:,:) :: values_2d
  real(4), allocatable, dimension(:,:,:) :: values_3d
  real(4), allocatable, dimension(:,:,:,:) :: values_4d
  real(4) mval,r4val
  integer ndim,nvar,ndims,ival,idate(6),icheck(6),ierr,n
  logical hasit
  
  dsetin = open_dataset('dynf000_template.nc.in')
  ! create a copy of the template file
  print *,'*** Test creation of new dataset from template...'
  dset = create_dataset('dynf000.nc',dsetin)
  print *,'*** Test that number of variables,dimensions,attributes read...'
  if (dsetin%nvars .ne. 23) then
    print *,'***number of variables not correct...'
    stop 99
  endif
  if (dsetin%ndims .ne. 5) then
    print *,'***number of dimensions not correct...'
    stop 99
  endif
  if (dsetin%natts .ne. 8) then
    print *,'***number of attributes not correct...'
    stop 99
  endif
  call close_dataset(dsetin)
  print *,'*** Test read of variable data...'
  call read_vardata(dset, 'pressfc', values_3d)
  call read_vardata(dset, 'vgrd', values_4d)
  values_3d=1.013e5 
  values_4d=99.
  ! populate pressfc and vgrd
  print *,'*** Test write of variable data...'
  call write_vardata(dset,'pressfc', values_3d)
  call write_vardata(dset,'vgrd', values_4d)
  call write_attribute(dset,'foo',1,'ugrd')
  if (allocated(values_1d)) deallocate(values_1d)
  allocate(values_1d(5))
  values_1d = (/1,2,3,4,5/)
  print *,'*** Test write of attributes...'
  call write_attribute(dset,'bar',values_1d,'ugrd')
  call write_attribute(dset,'hello','world')
  print *,'*** Test read of attribute just written...'
  call read_attribute(dset,'hello',charatt)
  if (trim(charatt) .ne. 'world') then
    print *,'***attribute not read or written correctly...'
    stop 99
  endif
  call read_attribute(dset,'missing_value',mval,'pressfc')
  print *,'*** Test read of missing_value attribute...'
  if (mval .ne. -1.e10) then
    print *,'***missing_value not correct...'
    stop 99
  endif
  nvar = get_nvar(dset, 'vgrd')
  print *,'*** Test getting variable id...'
  if (trim(dset%variables(nvar)%name) .ne. 'vgrd') then
    print *,'***variable id not correct...'
    stop 99
  endif
  ! print all info for this variable
  !do n=1,dset%variables(nvar)%ndims
  !   print *,trim(adjustl(dset%variables(nvar)%dimnames(n))),&
  !        dset%variables(nvar)%dimlens(n),&
  !        dset%dimensions(dset%variables(nvar)%dimindxs(n))%len
  !enddo
  idate = get_idate_from_time_units(dset)
  print *,'*** Test getting idate from time units attribute...'
  icheck = (/2016,1,4,6,0,0/)
  if (all(idate .ne. icheck)) then
      print *,'*** idate not correct'
      stop 99
  endif
  time_units = get_time_units_from_idate(idate, time_measure='minutes')
  print *,'*** Test getting time units from idate...'
  if (trim(time_units) .ne. 'minutes since 2016-01-04 06:00:00') then
    print *,'***time units not correct...'
    stop 99
  endif
  call read_attribute(dset,'foo',time_units,errcode=ierr)
  print *,'*** Test the return code is correct for missing attribute...'
  if (ierr .ne. NF90_ENOTATT) then
    print *,'***return code not correct for missing attribute...'
    stop 99
  endif
  call close_dataset(dset)
  
  ! re-open dataset
  dset = open_dataset('dynf000.nc')
  var = get_var(dset,'vgrd') 
  print *,'*** Test number of dimensions for variable is correct...'
  if (var%ndims .ne. 4) then
    print *,'***number of dimensions for variable not correct...'
    stop 99
  endif
  print *,'*** Test reading of data just written...'
  call read_vardata(dset, 'vgrd', values_4d)
  if (minval(values_4d) .ne. 99. .or. maxval(values_4d) .ne. 99.) then
    print *,'***vgrd variable data read as 4d not correct...'
    stop 99
  endif
  ! this should work also, since time dim is singleton
  call read_vardata(dset, 'vgrd', values_3d) 
  if (minval(values_3d) .ne. 99. .or. maxval(values_3d) .ne. 99.) then
    print *,'***vgrd variable data read as 3d not correct...'
    stop 99
  endif
  call read_vardata(dset, 'pressfc', values_3d)
  if (minval(values_3d) .ne. 1.013e5 .or. maxval(values_3d) .ne. 1.013e5) then
    print *,'***pressfc variable data read as 3d not correct...'
    stop 99
  endif
  ! this should work also, since time dim is singleton
  call read_vardata(dset, 'pressfc', values_2d)
  if (minval(values_2d) .ne. 1.013e5 .or. maxval(values_2d) .ne. 1.013e5) then
    print *,'***presssfc variable data read as 2d not correct...'
    stop 99
  endif
  print *, '*** Test has_var function...'
  hasit = has_var(dset,'pressfc')
  if (.not. hasit) then
    print *,'***has_var function failed...'
    stop 99
  endif
  print *,'*** Test has_att function...'
  hasit = has_attr(dset,'max_abs_compression_error','pressfc')
  if (hasit) then
    print *,'***has_att function failed...'
    stop 99
  endif
  print *,'*** Test reading of array attribute...'
  call read_attribute(dset,'ak',values_1d)
  if (values_1d(1) .ne. 20.) then
    print *,'***array attribute read failed...'
    stop 99
  endif
  call close_dataset(dset)
  print *,"SUCCESS!"
end program test_ncio
