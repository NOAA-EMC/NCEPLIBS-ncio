program test_ncio
  use netcdf
  use module_ncio
  implicit none
  
  character(len=72) charatt, time_units
  type(Dataset) :: dset, dsetin
  type(Variable) :: var
  real(4), allocatable, dimension(:) :: values_1d
  real(8), allocatable, dimension(:) :: values8_1d
  real(4), allocatable, dimension(:,:) :: values_2d
  real(4), allocatable, dimension(:,:,:) :: values_3d
  real(4), allocatable, dimension(:,:,:,:) :: values_4d
  real(4) mval,r4val
  integer ndim,nvar,ndims,ival,idate(6),ierr,n
  logical hasit
  
  dsetin = open_dataset('dynf000_template.nc')
  ! create a copy of the template file
  dset = create_dataset('dynf000.nc',dsetin)
  print *,'ncid=',dsetin%ncid
  print *,'nvars=',dsetin%nvars
  print *,'ndims=',dsetin%ndims
  print *,'natts=',dsetin%natts
  call close_dataset(dsetin)
  call read_vardata(dset, 'pressfc', values_3d)
  call read_vardata(dset, 'vgrd', values_4d)
  values_3d=1.013e5 
  values_4d=99.
  ! populate pressfc and vgrd
  call write_vardata(dset,'pressfc', values_3d)
  call write_vardata(dset,'vgrd', values_4d)
  call write_attribute(dset,'foo',1,'ugrd')
  if (allocated(values_1d)) deallocate(values_1d)
  allocate(values_1d(5))
  values_1d = (/1,2,3,4,5/)
  call write_attribute(dset,'bar',values_1d,'ugrd')
  call write_attribute(dset,'hello','world')
  call read_attribute(dset,'hello',charatt)
  print *,trim(charatt)
  call read_attribute(dset,'missing_value',mval,'pressfc')
  print *,'missing_value=',mval
  nvar = get_nvar(dset, 'vgrd')
  print *,trim(dset%variables(nvar)%name)
  do n=1,dset%variables(nvar)%ndims
     print *,trim(adjustl(dset%variables(nvar)%dimnames(n))),&
          dset%variables(nvar)%dimlens(n),&
          dset%dimensions(dset%variables(nvar)%dimindxs(n))%len
  enddo
  idate = get_idate_from_time_units(dset)
  print *,'idate=',idate
  time_units = get_time_units_from_idate(idate, time_measure='minutes')
  print *,'time_units=',trim(time_units)
  call read_attribute(dset,'foo',time_units,errcode=ierr)
  print *,'error code = ',ierr,'should be',NF90_ENOTATT
  call close_dataset(dset)
  
  ! re-open dataset
  dset = open_dataset('dynf000.nc')
  var = get_var(dset,'vgrd') 
  print *,'vgrd has ',var%ndims,' dims'
  call read_vardata(dset, 'vgrd', values_4d)
  print *,'min/max vgrd (4d)'
  print *,minval(values_4d),maxval(values_4d)
  ! this should work also, since time dim is singleton
  call read_vardata(dset, 'vgrd', values_3d) 
  print *,'min/max vgrd (3d)'
  print *,minval(values_3d),maxval(values_3d)
  print *,'min/max diff=',minval(values_4d(:,:,:,1)-values_3d),maxval(values_4d(:,:,:,1)-values_3d)
  ! this should work also, since time dim is singleton
  print *,'min/max pressfc (3d)'
  call read_vardata(dset, 'pressfc', values_3d)
  print *,minval(values_3d),maxval(values_3d)
  print *,'min/max pressfc (2d)'
  call read_vardata(dset, 'pressfc', values_2d)
  print *,minval(values_2d),maxval(values_2d)
  print *,'min/max diff=',minval(values_3d(:,:,1)-values_2d),maxval(values_3d(:,:,1)-values_2d)
  print *,'min/max pfull (1d_r8)'
  call read_vardata(dset, 'pfull', values8_1d)
  print *,minval(values8_1d),maxval(values8_1d)
  hasit = has_var(dset,'pressfc')
  print *,'has var pressfc',hasit
  hasit = has_attr(dset,'max_abs_compression_error','pressfc')
  print *,'pressfc has max_abs_compression_error attribute',hasit
  call read_attribute(dset,'ncnsto',ival)
  print *,'ncnsto =',ival
  call read_attribute(dset,'ak',values_1d)
  print *,'ak =',values_1d
  call close_dataset(dset)
end program test_ncio
