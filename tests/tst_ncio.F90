! This is a test for the NCEPLIBS-ncio package.
program tst_ncio

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
  real(4), allocatable, dimension(:,:,:,:,:) :: values_5d
  real(4), dimension(10,10) :: quantize1, quantize2
  real(4) mval,r4val,qerr
  integer ndim,nvar,ndims,ival,idate(6),icheck(6),ierr,n,nbits
  logical hasit

  print *,'*** Testing NCEPLIBS-ncio.'

  dsetin = open_dataset('dynf000_template.nc.in')
  print *,'*** Test creation of new dataset from template...'
  dset = create_dataset('dynf000.nc',dsetin)
  print *,'*** Test that number of variables,dimensions,attributes is read...'
  if (dsetin%nvars .ne. 25) then
     print *,'***number of variables not correct...'
     stop 99
  endif
  if (dsetin%ndims .ne. 7) then
     print *,'***number of dimensions not correct...'
     stop 99
  endif
  if (dsetin%natts .ne. 8) then
     print *,'***number of attributes not correct...'
     stop 99
  endif

  print *,'*** Test read of variable data...'
  call read_vardata(dsetin, 'pressfc', values_3d)
  call read_vardata(dsetin, 'vgrd', values_4d)
  call read_vardata(dsetin, 'tmp_spread', values_5d)
  if (maxval(values_3d) .ne. 102345.6) then
     print *,'*** read_vardata not working properly...'
     print *, 'maxvalue(pressfc) != 102345.6'
     stop 99
  end if
  if (minval(values_4d) .ne. -5.5) then
     print *,'*** read_vardata not working properly...'
     print *, 'minvalue(vgrd) != -5.5'
     stop 99
  end if
  if ((minval(values_5d) .ne. -1.0) .and. (maxval(values_5d) .ne. 1.0)) then
     print *,'*** read_vardata not working properly...'
     print *, 'maxvalue(tmp_spread) != -1.0 and maxvalue(tmp_spread) != 1.0'
     stop 99
  end if
  print *,'*** Test reading of slice for 5d var...'
  call read_vardata(dsetin, 'tmp_spread', values_5d, 10, 3)
  if ( all(shape(values_5d) .ne. (/256,128,1,2,1/)) ) then
     print *,'***shape of 5d slice incorrect...'
  endif
  call close_dataset(dsetin)
  values_3d=1.013e5
  values_4d=99.

  ! populate pressfc and vgrd
  print *,'*** Test write of variable data...'
  call write_vardata(dset,'pressfc', values_3d)
  call write_vardata(dset,'vgrd', values_4d)

  ! write just a slice along 3rd dimension (index 10)
  values_3d = -99.
  call write_vardata(dset, 'vgrd', values_3d, 10, 3)
  call write_attribute(dset,'foo',1,'ugrd')
  if (allocated(values_1d)) deallocate(values_1d)
  allocate(values_1d(5))
  values_1d = (/1,2,3,4,5/)

  ! the next tests exercise the quantize data functionality
  ! this packs the data so that the last 32-nbits of each float are
  ! zero for efficient zlib compression (note this results in 'lossy'
  ! compression and the qerr value is the max error introduced by
  ! this data packing procedure
  print *,'*** Test quantize data to 8 bit ...'
  nbits = 8 ! eight bit only
  quantize1 = 314.1592653589793
  quantize1(5,:) = 123.456789
  quantize1(10,:) = 1013.254321
  call quantize_data(quantize1, quantize2, nbits, qerr)
  if (abs(quantize2(1,1)-quantize1(1,1)) .gt. qerr) then
     print *,'*** quantize data not working properly...'
     print *,'abs(quantize2(1,1)-quantize1(1,1))==',abs(quantize2(1,1)-quantize1(1,1))
     print *,'error==',qerr
     stop 99
  end if

  print *,'*** Test quantize data to 32 bit ...'
  nbits = 32
  call quantize_data(quantize1, quantize2, nbits, qerr)
  if (abs(quantize2(1,1)-quantize1(1,1)) .gt. qerr) then
     print *,'*** quantize data not working properly...'
     print *,'abs(quantize2(1,1)-quantize1(1,1))==',abs(quantize2(1,1)-quantize1(1,1))
     print *,'error==',qerr
     stop 99
  end if

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

  print *,'*** Test variable dimension info...'
  if (trim(adjustl(dset%variables(nvar)%dimnames(1))) .ne. 'grid_xt' .or. &
       dset%variables(nvar)%dimlens(1) .ne. 256 .or. &
       dset%dimensions(dset%variables(nvar)%dimindxs(1))%len .ne. 256) then
     print *,'***dimension 1 info not correct...'
     stop 99
  endif
  if (trim(adjustl(dset%variables(nvar)%dimnames(2))) .ne. 'grid_yt' .or. &
       dset%variables(nvar)%dimlens(2) .ne. 128 .or. &
       dset%dimensions(dset%variables(nvar)%dimindxs(2))%len .ne. 128) then
     print *,'***dimension 2 info not correct...'
     stop 99
  endif
  if (trim(adjustl(dset%variables(nvar)%dimnames(3))) .ne. 'pfull' .or. &
       dset%variables(nvar)%dimlens(3) .ne. 64 .or. &
       dset%dimensions(dset%variables(nvar)%dimindxs(3))%len .ne. 64) then
     print *,'***dimension 3 info not correct...'
     stop 99
  endif
  if (trim(adjustl(dset%variables(nvar)%dimnames(4))) .ne. 'time' .or. &
       dset%variables(nvar)%dimlens(4) .ne. 1 .or. &
       dset%dimensions(dset%variables(nvar)%dimindxs(4))%len .ne. 1) then
     print *,'***dimension 4 info not correct...'
     stop 99
  endif

  print *,'*** Test getting idate from time units attribute...'
  idate = get_idate_from_time_units(dset)
  icheck = (/2016,1,4,6,0,0/)
  if (all(idate .ne. icheck)) then
     print *,'*** idate not correct'
     stop 99
  endif

  print *,'*** Test getting time units from idate...'
  time_units = get_time_units_from_idate(idate, time_measure='minutes')
  if (trim(time_units) .ne. 'minutes since 2016-01-04 06:00:00') then
     print *,'***time units not correct...'
     stop 99
  endif

  print *,'*** Test the return code is correct for missing attribute...'
  call read_attribute(dset,'foo',time_units,errcode=ierr)
  if (ierr .ne. NF90_ENOTATT) then
     print *,'***return code not correct for missing attribute...'
     stop 99
  endif
  call close_dataset(dset)

  ! re-open dataset
  dset = open_dataset('dynf000.nc')
  print *,'*** Test number of dimensions for variable is correct...'
  var = get_var(dset,'vgrd')
  if (var%ndims .ne. 4) then
     print *,'***number of dimensions for variable not correct...'
     stop 99
  endif

  print *,'*** Test reading of data just written...'
  call read_vardata(dset, 'vgrd', values_4d)
  if (minval(values_4d) .ne. -99. .or. maxval(values_4d) .ne. 99.) then
     print *,'***vgrd variable data read as 4d not correct...'
     stop 99
  endif

  ! this should work also, since time dim is singleton
  call read_vardata(dset, 'vgrd', values_3d)
  if (minval(values_3d) .ne. -99. .or. maxval(values_3d) .ne. 99.) then
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

  print *,'*** Test reading of slice for 3d var...'
  ! read 10th element along 3rd dimension
  call read_vardata(dset, 'vgrd', values_3d,10,3)
  if ( all(shape(values_3d) .ne. (/256,128,1/)) ) then
     print *,'***shape of 3d slice incorrect...'
  endif
  if ( all(values_3d .ne. -99.) ) then
     print *,'***data in 3d slice incorrect...'
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

end program tst_ncio
