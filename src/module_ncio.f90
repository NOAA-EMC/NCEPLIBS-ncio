!> @file
!! @brief Module for reading/writing netcdf files output by the GFS.
!! @author Jeff Whitaker <jeffrey.s.whitaker@noaa.gov> @date 201910

!> This is a module for reading/writing netcdf files output by the
!! GFS. It handles 32 and 64 bit real variables, 8, 16 and 32 bit
!! integer variables and char variables. Variables can have up to 5
!! dimensions.
!!
!! @note Wwriting requires a template file.
!!
!! @author Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
!! @date Oct, 2019
module module_ncio

  use netcdf
  use mpi
  
  implicit none
  private

  type Variable
     integer varid !< NetCDF variable ID.
     integer ndims !< Number of dimensions.
     integer dtype !< NetCDF data type.
     integer natts !< Number of attributes.
     integer deflate_level !< Compression level (if > 0).
     logical shuffle  !< Shuffle filter?
     logical hasunlim !< Has an unlimited dim?
     character(len=nf90_max_name) :: name !< Variable name.
     integer, allocatable, dimension(:) :: dimids !< NetCDF dimension IDs.
     integer, allocatable, dimension(:) :: dimindxs !< Indices into Dataset%dimensions for associated dimensions.
     character(len=nf90_max_name), allocatable, dimension(:) :: dimnames !< Names of associated dimensions.
     integer, allocatable, dimension(:) :: dimlens !< Current dimension lengths (updated after every write_vardata call).
     integer, allocatable, dimension(:) :: chunksizes !< Chunksizes to use with data.
  end type Variable

  type Dimension
     integer dimid !< netCDF dimension ID
     integer len   !< dimension length (updated after every write_vardata call)
     logical isunlimited !< unlimited? 
     character(len=nf90_max_name) :: name !< name name of dimension
  end type Dimension

  type Dataset
     integer :: ncid  !< netCDF ID.
     integer :: nvars !< number of variables in dataset
     integer :: ndims !< number of dimensions in dataset
     integer :: natts !< number of dataset (global) attributes
     integer :: nunlimdim !< dimension ID for unlimited dimension
     logical :: ishdf5 !< is underlying disk format HDF5?
     logical :: isparallel !< was file opened for parallel I/O?
     character(len=500) filename !< netCDF filename
     type(Variable), allocatable, dimension(:) :: variables !< array of Variable instances
     type(Dimension), allocatable, dimension(:) :: dimensions !< array of Dimension instances
  end type Dataset

  !> @defgroup read_vardata_8param Read Variable Data with slice/start/count Arrays.
  !! Read data from variable varname in dataset dset, return in it
  !! array values.
  !!
  !! @param[in] dset Input dataset instance returned by
  !! open_dataset/create_dataset.
  !! @param[in] varname Input string name of variable.
  !! @param[inout] values Array to hold variable data. Must be an allocatable
  !! array with same rank as variable varname (or 1 dimension less).
  !! @param[in] nslice optional index along dimension slicedim
  !! @param[in] slicedim optional, if nslice is set, index of which
  !! dimension to slice with nslice, default is ndims.
  !! @param[in] ncstart optional, if ncstart and nccount are set, manually
  !! specify the start and count of netCDF read.
  !! @param[in] nccount optional, if ncstart and nccount are set, manually
  !! specify the start and count of netCDF read.
  !! @param[out] errcode optional error return code. If not specified,
  !! program will stop if a nonzero error code returned from netcdf
  !! library.
  !! @author Jeff Whitaker
  !!
  !! @defgroup read_vardata_4param Read Variable Data without slice/start/count Arrays.
  !!
  !! When reading 5D arrays, the subsetting options are not present in
  !! NCEPLIBS-ncio.
  !!
  !! Read data from variable varname in dataset dset, return in it
  !! array values.
  !!
  !! @param[in] dset Input dataset instance returned by
  !! open_dataset/create_dataset.
  !! @param[in] varname Input string name of variable.
  !! @param[inout] values Array to hold variable data. Must be an allocatable
  !! array with same rank as variable varname (or 1 dimension less).
  !! @param[out] errcode optional error return code. If not specified,
  !! program will stop if a nonzero error code returned from netcdf
  !! library.
  !!
  !! @author Jeff Whitaker

  !> This interface provides access to the read_vardata
  !! subroutines. The 5D subroutines have a different number of
  !! parameters from the lower-dimension read subroutines.
  !!
  !! @author Jeff Whitaker
  interface read_vardata
     module procedure read_vardata_1d_r4, read_vardata_2d_r4, read_vardata_3d_r4,&
          read_vardata_4d_r4, read_vardata_5d_r4, &
          read_vardata_1d_r8, read_vardata_2d_r8, read_vardata_3d_r8,&
          read_vardata_4d_r8, read_vardata_5d_r8, &
          read_vardata_1d_int, read_vardata_2d_int, &
          read_vardata_3d_int, read_vardata_4d_int, read_vardata_5d_int, &
          read_vardata_1d_short, read_vardata_2d_short, &
          read_vardata_3d_short, read_vardata_4d_short, read_vardata_5d_short , &
          read_vardata_1d_byte, read_vardata_2d_byte, &
          read_vardata_3d_byte, read_vardata_4d_byte, read_vardata_5d_byte, &
          read_vardata_1d_char, read_vardata_2d_char, &
          read_vardata_3d_char, read_vardata_4d_char, read_vardata_5d_char
  end interface read_vardata
  
  !> Write data (in array values) to variable varname in dataset dset.
  !!
  !! @param[in] dset Input dataset instance returned by open_dataset/create_dataset.
  !! @param[in] varname Input string name of variable.
  !! @param[in] values Array with variable data. Must be an
  !! allocatable array with same rank as variable varname (or 1
  !! dimension less).
  !! @param[in] nslice optional index along dimension slicedim.
  !! @param[in] slicedim optional, if nslice is set, index of which
  !! dimension to slice with nslice, default is ndims.
  !! @param[in] ncstart optional, if ncstart and nccount are set,
  !! manually specify the start and count of netCDF write.
  !! @param[in] nccount optional, if ncstart and nccount are set,
  !! manually specify the start and count of netCDF write.
  !! @param[out] errcode optional error return code. If not specified,
  !! program will stop if a nonzero error code returned from netcdf
  !! library.
  !!
  !! @author Jeff Whitaker 
  interface write_vardata
     module procedure write_vardata_1d_r4, write_vardata_2d_r4, write_vardata_3d_r4,&
          write_vardata_4d_r4, write_vardata_1d_r8, write_vardata_2d_r8, write_vardata_3d_r8,&
          write_vardata_4d_r8, write_vardata_1d_int, write_vardata_2d_int, &
          write_vardata_3d_int, write_vardata_4d_int, &
          write_vardata_5d_int, write_vardata_5d_r4, write_vardata_5d_r8, &
          write_vardata_1d_short, write_vardata_2d_short, write_vardata_3d_short, &
          write_vardata_4d_short, write_vardata_5d_short, &
          write_vardata_1d_byte, write_vardata_2d_byte, write_vardata_3d_byte, &
          write_vardata_4d_byte, write_vardata_5d_byte, &
          write_vardata_1d_char, write_vardata_2d_char, write_vardata_3d_char, &
          write_vardata_4d_char, write_vardata_5d_char
  end interface write_vardata
    
  !> Read attribute 'attname' return in 'values'.
  !!
  !! If optional argument 'varname' is given, a variable attribute is
  !! returned. If the attribute is a 1d array, values should be an
  !! allocatable 1d array of the correct type.
  !!
  !! @param[in] dset Input dataset instance returned by open_dataset/create_dataset.
  !! @param[in] attname Input string name of attribute to read.
  !! @param[inout] values Array with attribute data. Must be an allocatable
  !! array with same rank as attribute attname.
  !! @param[in] varname optional, if provided, attribute will be read from variable
  !! varname, otherwise attribute will be assumed to be a global attribute.
  !! @param[out] errcode optional error return code. If not specified,
  !! program will stop if a nonzero error code returned from netcdf
  !! library.
  !!
  !! @author Jeff Whitaker
  interface read_attribute
     module procedure read_attribute_r4_scalar, read_attribute_int_scalar,&
          read_attribute_r8_scalar, read_attribute_r4_1d,&
          read_attribute_int_1d, read_attribute_r8_1d, read_attribute_char, &
          read_attribute_short_scalar, read_attribute_short_1d, &
          read_attribute_byte_scalar, read_attribute_byte_1d
  end interface read_attribute
  
  !> Write attribute 'attname' with data in 'values'. If optional
  !! argument 'varname' is given, a variable attribute is written.
  !! values can be a real(4), real(8), integer, string or 1d array.
  !!
  !! @param[in] dset Input dataset instance returned by open_dataset/create_dataset.
  !! @param[in] attname Input string name of attribute to write.
  !! @param[inout] values Array with attribute data. Must be an allocatable
  !! array with same rank as attribute attname.
  !! @param[in] varname optional, if provided, attribute will be written to variable
  !! varname, otherwise attribute will be assumed to be a global attribute.
  !! @param[out] errcode optional error return code. If not specified,
  !! program will stop if a nonzero error code returned from netcdf
  !! library.
  !!
  !! @author Jeff Whitaker
  interface write_attribute
     module procedure write_attribute_r4_scalar, write_attribute_int_scalar,&
          write_attribute_r8_scalar, write_attribute_r4_1d,&
          write_attribute_int_1d, write_attribute_r8_1d, write_attribute_char, &
          write_attribute_short_scalar, write_attribute_short_1d, &
          write_attribute_byte_scalar, write_attribute_byte_1d
  end interface write_attribute
  
  !> Quantize data.
  !!
  !! @param[in] dataIn Input array of data to quantize.
  !! @param[out] dataOut Output array of quantized data.
  !! @param[in] nbits Number of bits to quantize data by (0-32). This subroutine
  !! will convert data to 32 bit integers in range 0 to 2**nbits-1, then cast
  !! back to 32 bit floats (data is then quantized in steps proportional to
  !! to 2**nbits so last 32-nbits in floating point representation
  !! should be zero for efficient zlib compression).
  !! @param[out] compress_err Maximum absolute error between dataIn and dataOut.
  !!
  !! @author Jeff Whitaker, Cory Martin
  interface quantize_data
     module procedure quantize_data_2d, quantize_data_3d, &
          quantize_data_4d, quantize_data_5d
  end interface quantize_data
  public :: open_dataset, create_dataset, close_dataset, Dataset, Variable, Dimension, &
       read_vardata, read_attribute, write_vardata, write_attribute, get_ndim, &
       get_nvar, get_var, get_dim, get_idate_from_time_units, &
       get_time_units_from_idate, quantize_data, has_var, has_attr

contains

  !> Check return code, print error message.
  !!
  !! @param[in] status the status indicator
  !! @param[in] halt the halt option
  !! @param[in] fname the filename
  !! @author Jeff Whitaker 
  subroutine nccheck(status,halt,fname)
    implicit none
    integer, intent (in) :: status
    logical, intent(in), optional :: halt
    character(len=500), intent(in), optional :: fname
    logical stopit
    if (present(halt)) then
       stopit = halt
    else
       stopit = .true.
    endif
    if (status /= nf90_noerr) then
       write(0,*) status, trim(nf90_strerror(status))
       if (present(fname)) then
          write(0,*) trim(fname)
       end if
       if (stopit) stop 99
    end if
  end subroutine nccheck

  !> Get Dimension object given name.
  !! 
  !! @param dset the dataset
  !! @param dimname the dimension name
  !! @author Jeff Whitaker
  function get_dim(dset, dimname) result(dim)
    type(Dataset) :: dset
    type(Dimension) :: dim
    character(len=*), intent(in) :: dimname
    integer ndim
    ndim = get_ndim(dset, dimname)
    dim = dset%dimensions(ndim)
  end function get_dim

  !> Get Dimension index given name.
  !! Dimension object can then be accessed via Dataset%dimensions(nvar)
  !!
  !! @param dset the dataset
  !! @param dimname the dimension name
  !! 
  !! @return 
  !! @author Jeff Whitaker
  integer function get_ndim(dset, dimname)
    type(Dataset), intent(in) :: dset
    character(len=*), intent(in) :: dimname
    integer ndim
    get_ndim = -1
    do ndim=1,dset%ndims
       if (trim(dset%dimensions(ndim)%name) == trim(dimname)) then
          get_ndim = ndim
          exit
       endif
    enddo
  end function get_ndim

  !> Get Variable object given name.
  !!
  !! @param dset the dataset
  !! @param varname the variable name
  !! 
  !! @return 
  !! @author Jeff Whitaker
  function get_var(dset, varname) result (var)
    type(Dataset) :: dset
    type(Variable) :: var
    character(len=*) :: varname
    integer nvar
    nvar = get_nvar(dset, varname)
    var = dset%variables(nvar)
  end function get_var

  !> @return .true. is varname exists in dset, otherwise .false.
  !!
  !! @param dset the dataset
  !! @param varname the variable name
  !! 
  !! @return 
  logical function has_var(dset, varname)
    type(Dataset) :: dset
    character(len=*) :: varname
    integer nvar
    nvar = get_nvar(dset, varname)
    if (nvar > 0) then
       has_var=.true.
    else
       has_var=.false.
    endif
  end function has_var

  !> @return .true. if attribute exists in dset, otherwise .false.
  !! use optional kwarg varname to check for a variable attribute.
  !!
  !! @param dset 
  !! @param attname the attribute name
  !! @param varname the variable name
  !! 
  !! @return 
  !! @author Jeff Whitaker
  logical function has_attr(dset, attname, varname)
    type(Dataset) :: dset
    character(len=*) :: attname
    character(len=*), optional :: varname
    integer nvar, varid, ncerr
    nvar = get_nvar(dset, varname)
    if(present(varname))then
       nvar = get_nvar(dset,varname)
       if (nvar < 0) then
          has_attr = .false.
          return
       endif
       varid = dset%variables(nvar)%varid
    else
       varid = NF90_GLOBAL
    endif
    ncerr = nf90_inquire_attribute(dset%ncid, varid, attname)
    if (ncerr /= 0) then
       has_attr=.false.
    else
       has_attr=.true.
    endif
  end function has_attr

  !> Get Variable index given name.
  !!
  !! @param dset the dataset
  !! @param varname the variable name
  !! 
  !! @return
  !! @author Jeff Whitaker
  integer function get_nvar(dset,varname)
    type(Dataset), intent(in) :: dset
    character(len=*), intent(in) :: varname
    integer nvar
    get_nvar = -1
    do nvar=1,dset%nvars
       if (trim(dset%variables(nvar)%name) == trim(varname)) then
          get_nvar = nvar
          exit
       endif
    enddo
  end function get_nvar

  !> Reset dimension length (dimlens) for unlim dim for all variables.
  !!
  !! @param dset the dataset
  !! @param errcode the err code
  !! 
  !! @author Jeff Whitaker
  subroutine set_varunlimdimlens_(dset,errcode)
    type(Dataset), intent(inout) :: dset
    integer, intent(out), optional :: errcode
    integer ndim,n,nvar,ncerr
    logical return_errcode
    if(present(errcode)) then
       return_errcode=.true.
       errcode = 0
    else
       return_errcode=.false.
    endif
    ! loop over all vars
    do nvar=1,dset%nvars
       ! does var have unlim dimension?
       if (dset%variables(nvar)%hasunlim) then
          ! loop over all var dimensions
          do ndim=1,dset%variables(nvar)%ndims
             n = dset%variables(nvar)%dimindxs(ndim)
             ! n is the dimension index for this variable dimension
             ! if this dim is unlimited, update dimlens entry
             if (dset%dimensions(n)%isunlimited) then
                ncerr = nf90_inquire_dimension(dset%ncid,&
                     dset%dimensions(n)%dimid, &
                     len=dset%variables(nvar)%dimlens(ndim))
                if (return_errcode) then
                   call nccheck(ncerr,halt=.false.)
                   errcode=ncerr
                   return
                else
                   call nccheck(ncerr)
                endif
                ! also update len attribute of Dimension object
                dset%dimensions(n)%len = dset%variables(nvar)%dimlens(ndim)
             endif
          enddo
       endif
    enddo
  end subroutine set_varunlimdimlens_

  !> Open existing dataset, create dataset object for reading netcdf
  !! file.
  !!
  !! @param filename filename of netCDF Dataset.
  !! @param errcode optional error return code. If not specified
  !!          the program will stop if a nonzero error code returned by the netcdf lib.
  !! @param paropen optional flag to indicate whether to open dataset for parallel
  !!          access (Default .false.)
  !! @param mpicomm optional MPI communicator to use (Default MPI_COMM_WORLD)
  !!          ignored if paropen=F
  !!
  !! @returns Dataset object.
  !! @author Jeff Whitaker
  function open_dataset(filename,errcode,paropen, mpicomm) result(dset)
    implicit none
    character(len=*), intent(in) :: filename
    type(Dataset) :: dset
    integer, intent(out), optional :: errcode
    logical, intent(in), optional :: paropen
    integer, intent(in), optional :: mpicomm
    integer ncerr,nunlimdim,ndim,nvar,n,formatnum
    logical return_errcode
    if(present(errcode)) then
       return_errcode=.true.
       errcode = 0
    else
       return_errcode=.false.
    endif
    if (present(paropen)) then
       if (paropen) then
          dset%isparallel = .true.
       else
          dset%isparallel = .false.
       end if
    else
       dset%isparallel = .false.
    end if
    ! open netcdf file, get info, populate Dataset object.
    if (dset%isparallel) then
       if (present(mpicomm)) then
          ncerr = nf90_open(trim(filename), ior(NF90_NOWRITE, NF90_MPIIO), &
               comm=mpicomm, info = mpi_info_null, ncid=dset%ncid)
       else
          ncerr = nf90_open(trim(filename), ior(NF90_NOWRITE, NF90_MPIIO), &
               comm=mpi_comm_world, info = mpi_info_null, ncid=dset%ncid)
       end if
    else
       ncerr = nf90_open(trim(filename), NF90_NOWRITE, ncid=dset%ncid)
    end if
    if (return_errcode) then
       call nccheck(ncerr,halt=.false.,fname=filename)
       errcode=ncerr
       if (ncerr /= 0) return
    else
       call nccheck(ncerr,fname=filename)
    endif
    ncerr = nf90_inquire(dset%ncid, dset%ndims, dset%nvars, dset%natts, nunlimdim, formatnum=formatnum)
    if (return_errcode) then
       errcode=ncerr
       call nccheck(ncerr,halt=.false.,fname=filename)
       if (ncerr /= 0) return
    else
       call nccheck(ncerr,fname=filename)
    endif
    if (formatnum == nf90_format_netcdf4 .or. formatnum == nf90_format_netcdf4_classic) then
       dset%ishdf5 = .true.
    else
       dset%ishdf5 = .false.
    endif
    dset%filename = trim(filename)
    allocate(dset%variables(dset%nvars))
    allocate(dset%dimensions(dset%ndims))
    do ndim=1,dset%ndims
       dset%dimensions(ndim)%dimid = ndim
       ncerr = nf90_inquire_dimension(dset%ncid, ndim, name=dset%dimensions(ndim)%name, &
            len=dset%dimensions(ndim)%len)
       if (return_errcode) then
          errcode=ncerr
          call nccheck(ncerr,halt=.false.,fname=filename)
          if (ncerr /= 0) return
       else
          call nccheck(ncerr,fname=filename)
       endif
       if (ndim == nunlimdim) then
          dset%dimensions(ndim)%isunlimited = .true.
       else
          dset%dimensions(ndim)%isunlimited = .false.
       endif
    enddo
    do nvar=1,dset%nvars
       dset%variables(nvar)%hasunlim = .false.
       dset%variables(nvar)%varid = nvar
       ncerr = nf90_inquire_variable(dset%ncid, nvar,&
            name=dset%variables(nvar)%name,&
            natts=dset%variables(nvar)%natts,&
            xtype=dset%variables(nvar)%dtype,&
            ndims=dset%variables(nvar)%ndims)
       if (return_errcode) then
          errcode=ncerr
          call nccheck(ncerr,halt=.false.,fname=filename)
          if (ncerr /= 0) return
       else
          call nccheck(ncerr,fname=filename)
       endif
       allocate(dset%variables(nvar)%dimids(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%dimindxs(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%dimlens(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%chunksizes(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%dimnames(dset%variables(nvar)%ndims))
       if (dset%ishdf5) then
          ncerr = nf90_inquire_variable(dset%ncid, nvar,&
               dimids=dset%variables(nvar)%dimids,&
               deflate_level=dset%variables(nvar)%deflate_level,&
               chunksizes=dset%variables(nvar)%chunksizes,&
               shuffle=dset%variables(nvar)%shuffle)
       else
          ncerr = nf90_inquire_variable(dset%ncid, nvar,&
               dimids=dset%variables(nvar)%dimids)
       endif
       if (return_errcode) then
          errcode=ncerr
          call nccheck(ncerr,halt=.false.,fname=filename)
          if (ncerr /= 0) return
       else
          call nccheck(ncerr,fname=filename)
       endif
       do ndim=1,dset%variables(nvar)%ndims
          do n=1,dset%ndims
             if (dset%variables(nvar)%dimids(ndim) == dset%dimensions(n)%dimid) then
                exit
             endif
          enddo
          dset%variables(nvar)%dimindxs(ndim) = n
          dset%variables(nvar)%dimlens(ndim) = dset%dimensions(n)%len
          dset%variables(nvar)%dimnames(ndim) = dset%dimensions(n)%name
          if (dset%dimensions(n)%isunlimited) then
             dset%variables(nvar)%hasunlim = .true.
          endif
       enddo
    enddo
  end function open_dataset

  !> Create new dataset, using an existing dataset object to define.
  !! Variables, dimensions and attributes.
  !!
  !! @param[in] filename filename for netCDF Dataset.
  !! @param[in] dsetin dataset object to use as a template.
  !! @param[in] copy_vardata optional flag to control whether all variable
  !!              data is copied (Default is .false., only coordinate
  !!              variable data is copied).
  !! @param[in] paropen optional flag to indicate whether to open dataset for parallel
  !!          access (Default .false.)
  !! @param[in] nocompress optional flag to disable compression  even if input dataset is
  !!             compressed (Default .false.).
  !! @param[in] mpicomm optional MPI communicator to use (Default MPI_COMM_WORLD)
  !!          ignored if paropen=F
  !! @param[out] errcode optional error return code. If not specified
  !!          the program will stop if a nonzero error code returned by the netcdf lib.
  !!
  !! @returns Dataset object.
  !! @author Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
  function create_dataset(filename, dsetin, copy_vardata, paropen, nocompress, mpicomm, errcode) result(dset)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=nf90_max_name) :: attname, varname
    logical, intent(in), optional :: copy_vardata
    type(Dataset) :: dset
    type(Dataset), intent(in) :: dsetin
    logical, intent(in), optional :: paropen
    integer, intent(in), optional :: mpicomm
    logical, intent(in), optional :: nocompress
    integer, intent(out), optional :: errcode
    integer ncerr,ndim,nvar,n,ishuffle,natt
    logical copyd, coordvar, compress
    real(8), allocatable, dimension(:) :: values_1d
    real(8), allocatable, dimension(:,:) :: values_2d
    real(8), allocatable, dimension(:,:,:) :: values_3d
    real(8), allocatable, dimension(:,:,:,:) :: values_4d
    real(8), allocatable, dimension(:,:,:,:,:) :: values_5d
    integer, allocatable, dimension(:) :: ivalues_1d
    integer, allocatable, dimension(:,:) :: ivalues_2d
    integer, allocatable, dimension(:,:,:) :: ivalues_3d
    integer, allocatable, dimension(:,:,:,:) :: ivalues_4d
    integer, allocatable, dimension(:,:,:,:,:) :: ivalues_5d
    character, allocatable, dimension(:) :: cvalues_1d
    character, allocatable, dimension(:,:) :: cvalues_2d
    character, allocatable, dimension(:,:,:) :: cvalues_3d
    character, allocatable, dimension(:,:,:,:) :: cvalues_4d
    character, allocatable, dimension(:,:,:,:,:) :: cvalues_5d
    logical return_errcode
    if(present(errcode)) then
       return_errcode=.true.
       errcode = 0
    else
       return_errcode=.false.
    endif
    if (present(copy_vardata)) then
       ! copy all variable data
       copyd = .true.  
    else
       ! only copy coordinate variable data
       copyd = .false. 
    endif
    if (present(paropen)) then
       if (paropen) then
          dset%isparallel = .true.
       else
          dset%isparallel = .false.
       end if
    else
       dset%isparallel = .false.
    end if
    compress = .true.
    if (present(nocompress)) then
       if (nocompress) then
          compress = .false.
       end if
    end if
    ! create netcdf file
    if (dsetin%ishdf5) then
       if (dset%isparallel) then
          if (present(mpicomm)) then
             ncerr = nf90_create(trim(filename), &
                  cmode=IOR(NF90_CLOBBER,NF90_NETCDF4), ncid=dset%ncid, &
                  comm = mpicomm, info = mpi_info_null)
          else
             ncerr = nf90_create(trim(filename), &
                  cmode=IOR(NF90_CLOBBER,NF90_NETCDF4), ncid=dset%ncid, &
                  comm = mpi_comm_world, info = mpi_info_null)
          end if
       else
          ncerr = nf90_create(trim(filename), &
               cmode=IOR(NF90_CLOBBER,NF90_NETCDF4), &
               ncid=dset%ncid)
       end if
       dset%ishdf5 = .true.
    else
       ncerr = nf90_create(trim(filename), &
            cmode=IOR(NF90_CLOBBER,NF90_CDF5), &
            ncid=dset%ncid)
       dset%ishdf5 = .false.
    endif
    if (return_errcode) then
       errcode=ncerr
       call nccheck(ncerr,halt=.false.,fname=filename)
       if (ncerr /= 0) return
    else
       call nccheck(ncerr,fname=filename)
    endif
    ! copy global attributes
    do natt=1,dsetin%natts
       ncerr = nf90_inq_attname(dsetin%ncid, NF90_GLOBAL, natt, attname)
       if (return_errcode) then
          errcode=ncerr
          call nccheck(ncerr,halt=.false.)
          if (ncerr /= 0) return
       else
          call nccheck(ncerr)
       endif
       ncerr = nf90_copy_att(dsetin%ncid, NF90_GLOBAL, attname, dset%ncid, NF90_GLOBAL)
       if (return_errcode) then
          errcode=ncerr
          call nccheck(ncerr,halt=.false.)
          if (ncerr /= 0) return
       else
          call nccheck(ncerr)
       endif
    enddo
    dset%natts = dsetin%natts
    dset%filename = trim(filename)
    dset%ndims = dsetin%ndims
    dset%nvars = dsetin%nvars
    allocate(dset%variables(dsetin%nvars))
    allocate(dset%dimensions(dsetin%ndims))
    ! create dimensions
    do ndim=1,dsetin%ndims
       if (dsetin%dimensions(ndim)%isunlimited) then
          ncerr = nf90_def_dim(dset%ncid, trim(dsetin%dimensions(ndim)%name), &
               NF90_UNLIMITED, &
               dset%dimensions(ndim)%dimid)
          if (return_errcode) then
             errcode=ncerr
             call nccheck(ncerr,halt=.false.)
             if (ncerr /= 0) return
          else
             call nccheck(ncerr)
          endif
          dset%dimensions(ndim)%isunlimited = .true.
          dset%nunlimdim = ndim
          dset%dimensions(ndim)%len = 0
          dset%dimensions(ndim)%name = trim(dsetin%dimensions(ndim)%name)
       else
          ncerr = nf90_def_dim(dset%ncid, trim(dsetin%dimensions(ndim)%name),&
               dsetin%dimensions(ndim)%len, &
               dset%dimensions(ndim)%dimid)
          if (return_errcode) then
             errcode=ncerr
             call nccheck(ncerr,halt=.false.)
             if (ncerr /= 0) return
          else
             call nccheck(ncerr)
          endif
          dset%dimensions(ndim)%len = dsetin%dimensions(ndim)%len
          dset%dimensions(ndim)%isunlimited = .false.
          dset%dimensions(ndim)%name = trim(dsetin%dimensions(ndim)%name)
       endif
    enddo
    ! create variables
    do nvar=1,dsetin%nvars
       dset%variables(nvar)%hasunlim = .false.
       dset%variables(nvar)%ndims = dsetin%variables(nvar)%ndims
       allocate(dset%variables(nvar)%dimids(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%dimindxs(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%dimnames(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%dimlens(dset%variables(nvar)%ndims))
       allocate(dset%variables(nvar)%chunksizes(dset%variables(nvar)%ndims))
       dset%variables(nvar)%chunksizes = dsetin%variables(nvar)%chunksizes
       do ndim=1,dset%variables(nvar)%ndims
          do n=1,dset%ndims
             if (trim(dsetin%variables(nvar)%dimnames(ndim)) == &
                  trim(dset%dimensions(n)%name)) then
                exit
             endif
          enddo
          dset%variables(nvar)%dimindxs(ndim) = n
          dset%variables(nvar)%dimids(ndim) = dset%dimensions(n)%dimid
          dset%variables(nvar)%dimlens(ndim) = dset%dimensions(n)%len
          dset%variables(nvar)%dimnames(ndim) = dset%dimensions(n)%name
          if (dset%dimensions(n)%isunlimited) then
             dset%variables(nvar)%hasunlim = .true.
          endif
       enddo
       dset%variables(nvar)%name = dsetin%variables(nvar)%name
       dset%variables(nvar)%dtype = dsetin%variables(nvar)%dtype
       if (maxval(dset%variables(nvar)%chunksizes) > 0 .and. dset%ishdf5) then
          ! workaround for older versions of netcdf-fortran that don't
          ! like zero chunksize to be specified.
          ncerr = nf90_def_var(dset%ncid, &
               trim(dset%variables(nvar)%name),&
               dset%variables(nvar)%dtype, &
               dset%variables(nvar)%dimids, &
               dset%variables(nvar)%varid, &
               chunksizes=dset%variables(nvar)%chunksizes)
       else
          ncerr = nf90_def_var(dset%ncid, &
               trim(dset%variables(nvar)%name),&
               dset%variables(nvar)%dtype, &
               dset%variables(nvar)%dimids, &
               dset%variables(nvar)%varid)
       endif
       if (return_errcode) then
          errcode=ncerr
          call nccheck(ncerr,halt=.false.)
          if (ncerr /= 0) return
       else
          call nccheck(ncerr)
       endif
       if (dsetin%variables(nvar)%deflate_level > 0 .and. dset%ishdf5 .and. compress) then
          if (dsetin%variables(nvar)%shuffle) then
             ishuffle=1
          else
             ishuffle=0
          endif
          ncerr = nf90_def_var_deflate(dset%ncid, dset%variables(nvar)%varid,&
               ishuffle,1,dsetin%variables(nvar)%deflate_level)
          if (return_errcode) then
             errcode=ncerr
             call nccheck(ncerr,halt=.false.)
             if (ncerr /= 0) return
          else
             call nccheck(ncerr)
          endif
          dset%variables(nvar)%shuffle = dsetin%variables(nvar)%shuffle
          dset%variables(nvar)%deflate_level = &
               dsetin%variables(nvar)%deflate_level
       endif
       ! copy variable attributes
       do natt=1,dsetin%variables(nvar)%natts
          ncerr = nf90_inq_attname(dsetin%ncid, dsetin%variables(nvar)%varid, natt, attname)
          if (return_errcode) then
             errcode=ncerr
             call nccheck(ncerr,halt=.false.)
             if (ncerr /= 0) return
          else
             call nccheck(ncerr)
          endif
          if (.not. compress) then
             if (trim(attname) == 'max_abs_compression_error' &
                  .or. trim(attname) == 'nbits') then
                cycle
             end if
          end if
          ncerr = nf90_copy_att(dsetin%ncid, dsetin%variables(nvar)%varid, attname, dset%ncid, dset%variables(nvar)%varid)
          if (return_errcode) then
             errcode=ncerr
             call nccheck(ncerr,halt=.false.)
             if (ncerr /= 0) return
          else
             call nccheck(ncerr)
          endif
       enddo
    enddo
    ncerr = nf90_enddef(dset%ncid)
    if (return_errcode) then
       errcode=ncerr
       call nccheck(ncerr,halt=.false.)
       if (ncerr /= 0) return
    else
       call nccheck(ncerr)
    endif
    ! copy variable data
    ! assumes data is real (32 or 64 bit), or integer (16 or 32 bit) and 1-4d.
    do nvar=1,dsetin%nvars
       varname = trim(dsetin%variables(nvar)%name)
       ! is this variable a coordinate variable?
       coordvar = .false.
       if (trim(varname) == 'lats' .or. trim(varname) == 'lons' .or. &
            trim(varname) == 'lat'  .or. trim(varname) == 'lon') then
          coordvar = .true.
       else
          do ndim=1,dset%ndims
             if (trim(varname) == trim(dset%dimensions(ndim)%name)) then
                coordvar = .true.
             endif
          enddo
       endif
       ! if copy_data flag not given, and not a coordinate var,
       ! skip to next var.
       if (.not. coordvar .and. .not. copyd) cycle
       ! real variable
       if (dsetin%variables(nvar)%dtype == NF90_FLOAT .or.&
            dsetin%variables(nvar)%dtype == NF90_DOUBLE) then
          if (dsetin%variables(nvar)%ndims == 1) then
             call read_vardata(dsetin, varname, values_1d)
             call write_vardata(dset, varname, values_1d)
          else if (dsetin%variables(nvar)%ndims == 2) then
             call read_vardata(dsetin, varname, values_2d)
             call write_vardata(dset, varname, values_2d)
          else if (dsetin%variables(nvar)%ndims == 3) then
             call read_vardata(dsetin, varname, values_3d)
             call write_vardata(dset, varname, values_3d)
          else if (dsetin%variables(nvar)%ndims == 4) then
             call read_vardata(dsetin, varname, values_4d)
             call write_vardata(dset, varname, values_4d)
          else if (dsetin%variables(nvar)%ndims == 5) then
             call read_vardata(dsetin, varname, values_5d)
             call write_vardata(dset, varname, values_5d)
          endif
          ! integer var
       elseif (dsetin%variables(nvar)%dtype == NF90_INT .or.&
            dsetin%variables(nvar)%dtype == NF90_BYTE .or.&
            dsetin%variables(nvar)%dtype == NF90_SHORT) then
          ! TODO:  support NF90_UBYTE, USHORT, UINT, INT64, UINT64
          if (dsetin%variables(nvar)%ndims == 1) then
             call read_vardata(dsetin, varname, ivalues_1d)
             call write_vardata(dset, varname, ivalues_1d)
          else if (dsetin%variables(nvar)%ndims == 2) then
             call read_vardata(dsetin, varname, ivalues_2d)
             call write_vardata(dset, varname, ivalues_2d)
          else if (dsetin%variables(nvar)%ndims == 3) then
             call read_vardata(dsetin, varname, ivalues_3d)
             call write_vardata(dset, varname, ivalues_3d)
          else if (dsetin%variables(nvar)%ndims == 4) then
             call read_vardata(dsetin, varname, ivalues_4d)
             call write_vardata(dset, varname, ivalues_4d)
          else if (dsetin%variables(nvar)%ndims == 5) then
             call read_vardata(dsetin, varname, ivalues_5d)
             call write_vardata(dset, varname, ivalues_5d)
          endif
       elseif (dsetin%variables(nvar)%dtype == NF90_CHAR) then
          if (dsetin%variables(nvar)%ndims == 1) then
             call read_vardata(dsetin, varname, cvalues_1d)
             call write_vardata(dset, varname, cvalues_1d)
          else if (dsetin%variables(nvar)%ndims == 2) then
             call read_vardata(dsetin, varname, cvalues_2d)
             call write_vardata(dset, varname, cvalues_2d)
          else if (dsetin%variables(nvar)%ndims == 3) then
             call read_vardata(dsetin, varname, cvalues_3d)
             call write_vardata(dset, varname, cvalues_3d)
          else if (dsetin%variables(nvar)%ndims == 4) then
             call read_vardata(dsetin, varname, cvalues_4d)
             call write_vardata(dset, varname, cvalues_4d)
          else if (dsetin%variables(nvar)%ndims == 5) then
             call read_vardata(dsetin, varname, cvalues_5d)
             call write_vardata(dset, varname, cvalues_5d)
          endif
       else
          print *,'not copying variable ',trim(adjustl(varname)),&
               ' (unsupported data type or rank)'
       endif
    enddo
  end function create_dataset

  !> Close a netcdf file, deallocate members of dataset object. If
  !! optional error return code errcode is not specified, program will
  !! stop if a nonzero error code returned by the netcdf lib.
  !!
  !! @param dset a Dataset object with the open netCDF file.
  !! @param errcode optional error return code. If not specified the
  !! program will stop if a nonzero error code returned by the
  !! @author Jeff Whitaker
  subroutine close_dataset(dset,errcode)
    type(Dataset), intent(inout) :: dset
    integer, intent(out), optional :: errcode
    integer ncerr, nvar
    logical return_errcode
    if(present(errcode)) then
       return_errcode=.true.
       errcode = 0
    else
       return_errcode=.false.
    endif
    ncerr = nf90_close(ncid=dset%ncid)
    if (return_errcode) then
       errcode=ncerr
       call nccheck(ncerr,halt=.false.)
       if (ncerr /= 0) return
    else
       call nccheck(ncerr)
    endif
    do nvar=1,dset%nvars
       deallocate(dset%variables(nvar)%dimids)
       deallocate(dset%variables(nvar)%dimindxs)
       deallocate(dset%variables(nvar)%dimlens)
       deallocate(dset%variables(nvar)%chunksizes)
       deallocate(dset%variables(nvar)%dimnames)
    enddo
    deallocate(dset%variables,dset%dimensions)
  end subroutine close_dataset

  !> @copydoc read_vardata_8param
  subroutine read_vardata_1d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4), allocatable, dimension(:), intent(inout) :: values
    include "read_vardata_code_1d.f90"
  end subroutine read_vardata_1d_r4

  !> @copydoc read_vardata_8param
  subroutine read_vardata_2d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4), allocatable, dimension(:,:), intent(inout) :: values
    include "read_vardata_code_2d.f90"
  end subroutine read_vardata_2d_r4

  !> @copydoc read_vardata_8param
  subroutine read_vardata_3d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4), allocatable, dimension(:,:,:), intent(inout) :: values
    include "read_vardata_code_3d.f90"
  end subroutine read_vardata_3d_r4

  !> @copydoc read_vardata_8param
  subroutine read_vardata_4d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4), allocatable, dimension(:,:,:,:), intent(inout) :: values
    include "read_vardata_code_4d.f90"
  end subroutine read_vardata_4d_r4

  !> @copydoc read_vardata_4param
  subroutine read_vardata_5d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4), allocatable, dimension(:,:,:,:,:), intent(inout) :: values
    include "read_vardata_code_5d.f90"
  end subroutine read_vardata_5d_r4

  !> @copydoc read_vardata_8param
  subroutine read_vardata_1d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8), allocatable, dimension(:), intent(inout) :: values
    include "read_vardata_code_1d.f90"
  end subroutine read_vardata_1d_r8

  !> @copydoc read_vardata_8param
  subroutine read_vardata_2d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8), allocatable, dimension(:,:), intent(inout) :: values
    include "read_vardata_code_2d.f90"
  end subroutine read_vardata_2d_r8

  !> @copydoc read_vardata_8param
  subroutine read_vardata_3d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8), allocatable, dimension(:,:,:), intent(inout) :: values
    include "read_vardata_code_3d.f90"
  end subroutine read_vardata_3d_r8

  !> @copydoc read_vardata_8param
  subroutine read_vardata_4d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8), allocatable, dimension(:,:,:,:), intent(inout) :: values
    include "read_vardata_code_4d.f90"
  end subroutine read_vardata_4d_r8

  !> @copydoc read_vardata_4param
  subroutine read_vardata_5d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8), allocatable, dimension(:,:,:,:,:), intent(inout) :: values
    include "read_vardata_code_5d.f90"
  end subroutine read_vardata_5d_r8

  !> @copydoc read_vardata_8param
  subroutine read_vardata_1d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer, allocatable, dimension(:), intent(inout) :: values
    include "read_vardata_code_1d.f90"
  end subroutine read_vardata_1d_int

  !> @copydoc read_vardata_8param
  subroutine read_vardata_2d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer, allocatable, dimension(:,:), intent(inout) :: values
    include "read_vardata_code_2d.f90"
  end subroutine read_vardata_2d_int

  !> @copydoc read_vardata_8param
  subroutine read_vardata_3d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer, allocatable, dimension(:,:,:), intent(inout) :: values
    include "read_vardata_code_3d.f90"
  end subroutine read_vardata_3d_int

  !> @copydoc read_vardata_8param
  subroutine read_vardata_4d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer, allocatable, dimension(:,:,:,:), intent(inout) :: values
    include "read_vardata_code_4d.f90"
  end subroutine read_vardata_4d_int

  !> @copydoc read_vardata_4param
  subroutine read_vardata_5d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer, allocatable, dimension(:,:,:,:,:), intent(inout) :: values
    include "read_vardata_code_5d.f90"
  end subroutine read_vardata_5d_int

  !> @copydoc read_vardata_8param
  subroutine read_vardata_1d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2), allocatable, dimension(:), intent(inout) :: values
    include "read_vardata_code_1d.f90"
  end subroutine read_vardata_1d_short

  !> @copydoc read_vardata_8param
  subroutine read_vardata_2d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2), allocatable, dimension(:,:), intent(inout) :: values
    include "read_vardata_code_2d.f90"
  end subroutine read_vardata_2d_short

  !> @copydoc read_vardata_8param
  subroutine read_vardata_3d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2), allocatable, dimension(:,:,:), intent(inout) :: values
    include "read_vardata_code_3d.f90"
  end subroutine read_vardata_3d_short

  !> @copydoc read_vardata_8param
  subroutine read_vardata_4d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2), allocatable, dimension(:,:,:,:), intent(inout) :: values
    include "read_vardata_code_4d.f90"
  end subroutine read_vardata_4d_short

  !> @copydoc read_vardata_4param
  subroutine read_vardata_5d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2), allocatable, dimension(:,:,:,:,:), intent(inout) :: values
    include "read_vardata_code_5d.f90"
  end subroutine read_vardata_5d_short

  !> @copydoc read_vardata_8param
  subroutine read_vardata_1d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1), allocatable, dimension(:), intent(inout) :: values
    include "read_vardata_code_1d.f90"
  end subroutine read_vardata_1d_byte

  !> @copydoc read_vardata_8param
  subroutine read_vardata_2d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1), allocatable, dimension(:,:), intent(inout) :: values
    include "read_vardata_code_2d.f90"
  end subroutine read_vardata_2d_byte

  !> @copydoc read_vardata_8param
  subroutine read_vardata_3d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1), allocatable, dimension(:,:,:), intent(inout) :: values
    include "read_vardata_code_3d.f90"
  end subroutine read_vardata_3d_byte

  !> @copydoc read_vardata_8param
  subroutine read_vardata_4d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1), allocatable, dimension(:,:,:,:), intent(inout) :: values
    include "read_vardata_code_4d.f90"
  end subroutine read_vardata_4d_byte

  !> @copydoc read_vardata_4param
  subroutine read_vardata_5d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1), allocatable, dimension(:,:,:,:,:), intent(inout) :: values
    include "read_vardata_code_5d.f90"
  end subroutine read_vardata_5d_byte

  !> @copydoc read_vardata_8param
  subroutine read_vardata_1d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character, allocatable, dimension(:), intent(inout) :: values
    include "read_vardata_code_1d.f90"
  end subroutine read_vardata_1d_char

  !> @copydoc read_vardata_8param
  subroutine read_vardata_2d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character, allocatable, dimension(:,:), intent(inout) :: values
    include "read_vardata_code_2d.f90"
  end subroutine read_vardata_2d_char

  !> @copydoc read_vardata_8param
  subroutine read_vardata_3d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character, allocatable, dimension(:,:,:), intent(inout) :: values
    include "read_vardata_code_3d.f90"
  end subroutine read_vardata_3d_char

  !> @copydoc read_vardata_8param
  subroutine read_vardata_4d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character, allocatable, dimension(:,:,:,:), intent(inout) :: values
    include "read_vardata_code_4d.f90"
  end subroutine read_vardata_4d_char

  !> @copydoc read_vardata_4param
  subroutine read_vardata_5d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character, allocatable, dimension(:,:,:,:,:), intent(inout) :: values
    include "read_vardata_code_5d.f90"
  end subroutine read_vardata_5d_char

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_1d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4),  dimension(:), intent(in) :: values
    integer, intent(in), optional :: ncstart(1)
    integer, intent(in), optional :: nccount(1)
    include "write_vardata_code.f90"
  end subroutine write_vardata_1d_r4

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_2d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4),  dimension(:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(2)
    integer, intent(in), optional :: nccount(2)
    include "write_vardata_code.f90"
  end subroutine write_vardata_2d_r4

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_3d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4),  dimension(:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(3)
    integer, intent(in), optional :: nccount(3)
    include "write_vardata_code.f90"
  end subroutine write_vardata_3d_r4

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_4d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4),  dimension(:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(4)
    integer, intent(in), optional :: nccount(4)
    include "write_vardata_code.f90"
  end subroutine write_vardata_4d_r4

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_5d_r4(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(4),  dimension(:,:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(5)
    integer, intent(in), optional :: nccount(5)
    include "write_vardata_code.f90"
  end subroutine write_vardata_5d_r4

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_1d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8),  dimension(:), intent(in) :: values
    integer, intent(in), optional :: ncstart(1)
    integer, intent(in), optional :: nccount(1)
    include "write_vardata_code.f90"
  end subroutine write_vardata_1d_r8

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_2d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8),  dimension(:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(2)
    integer, intent(in), optional :: nccount(2)
    include "write_vardata_code.f90"
  end subroutine write_vardata_2d_r8

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_3d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8),  dimension(:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(3)
    integer, intent(in), optional :: nccount(3)
    include "write_vardata_code.f90"
  end subroutine write_vardata_3d_r8

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_4d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8),  dimension(:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(4)
    integer, intent(in), optional :: nccount(4)
    include "write_vardata_code.f90"
  end subroutine write_vardata_4d_r8

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_5d_r8(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    real(8),  dimension(:,:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(5)
    integer, intent(in), optional :: nccount(5)
    include "write_vardata_code.f90"
  end subroutine write_vardata_5d_r8

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_1d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer,  dimension(:), intent(in) :: values
    integer, intent(in), optional :: ncstart(1)
    integer, intent(in), optional :: nccount(1)
    include "write_vardata_code.f90"
  end subroutine write_vardata_1d_int

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_2d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer,  dimension(:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(2)
    integer, intent(in), optional :: nccount(2)
    include "write_vardata_code.f90"
  end subroutine write_vardata_2d_int

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_3d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer,  dimension(:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(3)
    integer, intent(in), optional :: nccount(3)
    include "write_vardata_code.f90"
  end subroutine write_vardata_3d_int

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_4d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer,  dimension(:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(4)
    integer, intent(in), optional :: nccount(4)
    include "write_vardata_code.f90"
  end subroutine write_vardata_4d_int

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_5d_int(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer,  dimension(:,:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(5)
    integer, intent(in), optional :: nccount(5)
    include "write_vardata_code.f90"
  end subroutine write_vardata_5d_int

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_1d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2),  dimension(:), intent(in) :: values
    integer, intent(in), optional :: ncstart(1)
    integer, intent(in), optional :: nccount(1)
    include "write_vardata_code.f90"
  end subroutine write_vardata_1d_short

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_2d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2),  dimension(:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(2)
    integer, intent(in), optional :: nccount(2)
    include "write_vardata_code.f90"
  end subroutine write_vardata_2d_short

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_3d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2),  dimension(:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(3)
    integer, intent(in), optional :: nccount(3)
    include "write_vardata_code.f90"
  end subroutine write_vardata_3d_short

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_4d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2),  dimension(:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(4)
    integer, intent(in), optional :: nccount(4)
    include "write_vardata_code.f90"
  end subroutine write_vardata_4d_short

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_5d_short(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(2),  dimension(:,:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(5)
    integer, intent(in), optional :: nccount(5)
    include "write_vardata_code.f90"
  end subroutine write_vardata_5d_short

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_1d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1),  dimension(:), intent(in) :: values
    integer, intent(in), optional :: ncstart(1)
    integer, intent(in), optional :: nccount(1)
    include "write_vardata_code.f90"
  end subroutine write_vardata_1d_byte

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_2d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1),  dimension(:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(2)
    integer, intent(in), optional :: nccount(2)
    include "write_vardata_code.f90"
  end subroutine write_vardata_2d_byte

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_3d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1),  dimension(:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(3)
    integer, intent(in), optional :: nccount(3)
    include "write_vardata_code.f90"
  end subroutine write_vardata_3d_byte

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_4d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1),  dimension(:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(4)
    integer, intent(in), optional :: nccount(4)
    include "write_vardata_code.f90"
  end subroutine write_vardata_4d_byte

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_5d_byte(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    integer(1),  dimension(:,:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(5)
    integer, intent(in), optional :: nccount(5)
    include "write_vardata_code.f90"
  end subroutine write_vardata_5d_byte

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_1d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character,  dimension(:), intent(in) :: values
    integer, intent(in), optional :: ncstart(1)
    integer, intent(in), optional :: nccount(1)
    include "write_vardata_code.f90"
  end subroutine write_vardata_1d_char

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_2d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character,  dimension(:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(2)
    integer, intent(in), optional :: nccount(2)
    include "write_vardata_code.f90"
  end subroutine write_vardata_2d_char

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_3d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character,  dimension(:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(3)
    integer, intent(in), optional :: nccount(3)
    include "write_vardata_code.f90"
  end subroutine write_vardata_3d_char

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_4d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character,  dimension(:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(4)
    integer, intent(in), optional :: nccount(4)
    include "write_vardata_code.f90"
  end subroutine write_vardata_4d_char

  !> @copydoc module_ncio::write_vardata
  subroutine write_vardata_5d_char(dset, varname, values, nslice, slicedim, ncstart, nccount, errcode)
    character,  dimension(:,:,:,:,:), intent(in) :: values
    integer, intent(in), optional :: ncstart(5)
    integer, intent(in), optional :: nccount(5)
    include "write_vardata_code.f90"
  end subroutine write_vardata_5d_char

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_int_scalar(dset, attname, values, varname, errcode)
    integer, intent(inout) :: values
    include "read_scalar_attribute_code.f90"
  end subroutine read_attribute_int_scalar

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_short_scalar(dset, attname, values, varname, errcode)
    integer(2), intent(inout) :: values
    include "read_scalar_attribute_code.f90"
  end subroutine read_attribute_short_scalar

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_byte_scalar(dset, attname, values, varname, errcode)
    integer(1), intent(inout) :: values
    include "read_scalar_attribute_code.f90"
  end subroutine read_attribute_byte_scalar

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_r4_scalar(dset, attname, values, varname, errcode)
    real(4), intent(inout) :: values
    include "read_scalar_attribute_code.f90"
  end subroutine read_attribute_r4_scalar

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_r8_scalar(dset, attname, values, varname, errcode)
    real(8), intent(inout) :: values
    include "read_scalar_attribute_code.f90"
  end subroutine read_attribute_r8_scalar

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_r4_1d(dset, attname, values, varname, errcode)
    real(4), intent(inout), allocatable, dimension(:) :: values
    include "read_attribute_code.f90"
  end subroutine read_attribute_r4_1d

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_r8_1d(dset, attname, values, varname, errcode)
    real(8), intent(inout), allocatable, dimension(:) :: values
    include "read_attribute_code.f90"
  end subroutine read_attribute_r8_1d

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_int_1d(dset, attname, values, varname, errcode)
    integer, intent(inout), allocatable, dimension(:) :: values
    include "read_attribute_code.f90"
  end subroutine read_attribute_int_1d

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_short_1d(dset, attname, values, varname, errcode)
    integer(2), intent(inout), allocatable, dimension(:) :: values
    include "read_attribute_code.f90"
  end subroutine read_attribute_short_1d

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_byte_1d(dset, attname, values, varname, errcode)
    integer(1), intent(inout), allocatable, dimension(:) :: values
    include "read_attribute_code.f90"
  end subroutine read_attribute_byte_1d

  !> @copydoc module_ncio::read_attribute
  subroutine read_attribute_char(dset, attname, values, varname, errcode)
    character(len=*), intent(inout) :: values
    include "read_scalar_attribute_code.f90"
  end subroutine read_attribute_char

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_int_scalar(dset, attname, values, varname, errcode)
    integer, intent(in) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_int_scalar

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_short_scalar(dset, attname, values, varname, errcode)
    integer(2), intent(in) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_short_scalar

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_byte_scalar(dset, attname, values, varname, errcode)
    integer(1), intent(in) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_byte_scalar

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_r4_scalar(dset, attname, values, varname, errcode)
    real(4), intent(in) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_r4_scalar

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_r8_scalar(dset, attname, values, varname, errcode)
    real(8), intent(in) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_r8_scalar

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_r4_1d(dset, attname, values, varname, errcode)
    real(4), intent(in), allocatable, dimension(:) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_r4_1d

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_r8_1d(dset, attname, values, varname, errcode)
    real(8), intent(in), allocatable, dimension(:) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_r8_1d

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_int_1d(dset, attname, values, varname, errcode)
    integer, intent(in), allocatable, dimension(:) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_int_1d

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_short_1d(dset, attname, values, varname, errcode)
    integer(2), intent(in), allocatable, dimension(:) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_short_1d

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_byte_1d(dset, attname, values, varname, errcode)
    integer(1), intent(in), allocatable, dimension(:) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_byte_1d

  !> @copydoc module_ncio::write_attribute
  subroutine write_attribute_char(dset, attname, values, varname, errcode)
    character(len=*), intent(in) :: values
    include "write_attribute_code.f90"
  end subroutine write_attribute_char

  !> return integer array with year,month,day,hour,minute,second
  !! parsed from time units attribute.
  !!
  !! @param[in] dset Input dataset instance returned by
  !! open_dataset/create_dataset.
  !!
  !! @author Jeff Whitaker
  function get_idate_from_time_units(dset) result(idate)
    type(Dataset), intent(in) :: dset
    integer idate(6)
    character(len=nf90_max_name) :: time_units
    integer ipos1,ipos2
    call read_attribute(dset, 'units', time_units, 'time')
    ipos1 = scan(time_units,"since",back=.true.)+1
    ipos2 = scan(time_units,"-",back=.false.)-1
    read(time_units(ipos1:ipos2),*) idate(1)
    ipos1 = ipos2+2; ipos2=ipos1+1
    read(time_units(ipos1:ipos2),*) idate(2)
    ipos1 = ipos2+2; ipos2=ipos1+1
    read(time_units(ipos1:ipos2),*) idate(3)
    ipos1 = scan(time_units,":")-2
    ipos2 = ipos1+1
    read(time_units(ipos1:ipos2),*) idate(4)
    ipos1 = ipos2+2
    ipos2 = ipos1+1
    read(time_units(ipos1:ipos2),*) idate(5)
    ipos1 = ipos2+2
    ipos2 = ipos1+1
    read(time_units(ipos1:ipos2),*) idate(6)
  end function get_idate_from_time_units

  !> create time units attribute of form 'hours since YYYY-MM-DD HH:MM:SS'
  !! from integer array with year,month,day,hour,minute,second
  !! optional argument 'time_measure' can be used to change 'hours' to
  !! 'days', 'minutes', 'seconds' etc.
  !!
  !! @param[in] idate Input array of integers idate(Year, Month, Day, Hour, Minute, Second)
  !! @param[in] time_measure optional, string to indicate units of time since
  !! idate time ('hours', 'days', 'minutes', 'seconds', etc.).
  !!
  !! @author Jeff Whitaker
  function get_time_units_from_idate(idate, time_measure) result(time_units)
    character(len=*), intent(in), optional :: time_measure
    integer, intent(in) ::  idate(6)
    character(len=12) :: timechar
    character(len=nf90_max_name) :: time_units
    if (present(time_measure)) then
       timechar = trim(time_measure)
    else
       timechar = 'hours'
    endif
    write(time_units,101) idate
101 format(' since ',i4.4,'-',i2.2,'-',i2.2,' ',&
         i2.2,':',i2.2,':',i2.2)
    time_units = trim(adjustl(timechar))//time_units
  end function get_time_units_from_idate

  !> @copydoc module_ncio::quantize_data
  subroutine quantize_data_2d(dataIn, dataOut, nbits, compress_err)
    real(4), intent(in) :: dataIn(:,:)
    real(4), intent(out) :: dataOut(:,:)
    include "quantize_data_code.f90"
  end subroutine quantize_data_2d

  !> @copydoc module_ncio::quantize_data
  subroutine quantize_data_3d(dataIn, dataOut, nbits, compress_err)
    real(4), intent(in) :: dataIn(:,:,:)
    real(4), intent(out) :: dataOut(:,:,:)
    include "quantize_data_code.f90"
  end subroutine quantize_data_3d

  !> @copydoc module_ncio::quantize_data
  subroutine quantize_data_4d(dataIn, dataOut, nbits, compress_err)
    real(4), intent(in) :: dataIn(:,:,:,:)
    real(4), intent(out) :: dataOut(:,:,:,:)
    include "quantize_data_code.f90"
  end subroutine quantize_data_4d

  !> @copydoc module_ncio::quantize_data
  subroutine quantize_data_5d(dataIn, dataOut, nbits, compress_err)
    real(4), intent(in) :: dataIn(:,:,:,:,:)
    real(4), intent(out) :: dataOut(:,:,:,:,:)
    include "quantize_data_code.f90"
  end subroutine quantize_data_5d

end module module_ncio
