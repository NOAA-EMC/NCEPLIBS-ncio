@mainpage

## NCEPLIBS-ncio

NetCDF read/write modules for the NCEP models. This is part of the
[NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) project.

The NCEPLIBS-ncio code here: https://github.com/NOAA-EMC/NCEPLIBS-ncio.

# API

```
  type Variable
     integer varid ! netCDF variable ID
     integer ndims ! number of dimensions
     integer dtype ! netCDF data type
     integer natts ! number of attributes
     integer deflate_level ! compression level (if > 0)
     logical shuffle  ! shuffle filter?
     logical hasunlim ! has an unlimited dim?
     character(len=nf90_max_name) :: name ! variable name
     integer, allocatable, dimension(:) :: dimids ! netCDF dimension IDs
     ! indices into Dataset%dimensions for associated dimensions.
     integer, allocatable, dimension(:) :: dimindxs
     ! names of associated dimensions.
     character(len=nf90_max_name), allocatable, dimension(:) :: dimnames
     ! current dimension lengths (updated after every write_vardata call)
     integer, allocatable, dimension(:) :: dimlens
     integer, allocatable, dimension(:) :: chunksizes
  end type Variable

  type Dimension
     integer dimid ! netCDF dimension ID
     integer len ! dimension length (updated after every write_vardata call)
     logical isunlimited ! unlimited?
     character(len=nf90_max_name) :: name ! name of dimension
  end type Dimension

  type Dataset
     integer :: ncid ! netCDF ID.
     integer :: nvars ! number of variables in dataset
     integer :: ndims ! number of dimensions in dataset
     integer :: natts ! number of dataset (global) attributes
     integer :: nunlimdim ! dimension ID for unlimited dimension
     logical :: ishdf5 ! is underlying disk format HDF5?
     logical :: isparallel ! was file opened for parallel I/O?
     character(len=500) filename ! netCDF filename
     ! array of Variable instances
     type(Variable), allocatable, dimension(:) :: variables
     ! array of Dimension instances
     type(Dimension), allocatable, dimension(:) :: dimensions
  end type Dataset

  logical function has_var(dset, varname)
    ! returns .true. is varname exists in dset, otherwise .false.

  logical function has_attr(dset, attname, varname)
    ! returns .true. if attribute exists in dset, otherwise .false.
    ! use optional kwarg varname to check for a variable attribute.

  function get_dim(dset, dimname) result(dim)
    ! get Dimension object given name

  integer function get_ndim(dset,varname)
    ! get Dimension index (ndim) given name
    ! Dimension object can then be accessed via Dataset%dimensions(ndim)

  function get_var(dset, varname) result (var)
    ! get Variable object given name

  integer function get_nvar(dset,varname)
    ! get variable index (nvar) given name
    ! Variable object can then be accessed via Dataset%variables(nvar)

  function open_dataset(filename,errcode,paropen, mpicomm) result(dset)
    ! open existing dataset, create dataset object for reading netcdf file
    !
    ! filename: filename of netCDF Dataset.
    ! errcode: optional error return code.  If not specified 
    !          the program will stop if a nonzero error code returned by the netcdf lib.
    ! paropen: optional flag to indicate whether to open dataset for parallel
    !          access (Default .false.)
    ! mpicomm: optional MPI communicator to use (Default MPI_COMM_WORLD)
    !          ignored if paropen=F
    !
    ! returns Dataset object.

  function create_dataset(filename, dsetin, copy_vardata, paropen, nocompress, mpicomm, errcode) result(dset)
    ! create new dataset, using an existing dataset object to define
    ! variables, dimensions and attributes.
    !
    ! filename: filename for netCDF Dataset.
    ! dsetin:  dataset object to use as a template.
    ! copyvardata: optional flag to control whether all variable
    !              data is copied (Default is .false., only coordinate
    !              variable data is copied).
    ! errcode: optional error return code.  If not specified 
    !          the program will stop if a nonzero error code returned by the netcdf lib.
    ! paropen: optional flag to indicate whether to open dataset for parallel
    !          access (Default .false.)
    ! nocompress: optional flag to disable compression  even if input dataset is
    !             compressed (Default .false.).
    ! mpicomm: optional MPI communicator to use (Default MPI_COMM_WORLD)
    !          ignored if paropen=F
    !
    ! returns Dataset object.

  subroutine close_dataset(dset,errcode)
    ! close netcdf file, deallocate members of dataset object.
    ! if optional error return code errcode is not specified,
    ! program will stop if a nonzero error code returned by the netcdf lib.

   subroutine read_vardata(dset,varname,values,nslice,slicedim,errcode)
    ! read data from variable varname in dataset dset, return in it array values.
    !
    ! dset:    Input dataset instance returned by open_dataset/create_dataset.
    ! varname: Input string name of variable.
    ! values:  Array to hold variable data.  Must be
    !          an allocatable array with same rank
    !          as variable varname (or 1 dimension less).
    ! nslice:  optional index along dimension slicedim
    ! slicedim: optional, if nslice is set, index of which dimension to slice with
    !          nslice, default is ndims
    ! ncstart: optional, if ncstart and nccount are set, manually specify the
    !          start and count of netCDF read
    ! nccount: optional, if ncstart and nccount are set, manually specify the
    !          start and count of netCDF read
    ! errcode: optional error return code.  If not specified,
    !          program will stop if a nonzero error code returned
    !          from netcdf library.

  subroutine write_vardata(dset,varname,values,nslice,slicedim,errcode)
    ! write data (in array values) to variable varname in dataset dset.
    !
    ! dset:    Input dataset instance returned by open_dataset/create_dataset.
    ! varname: Input string name of variable.
    ! values:  Array with variable data.  Must be
    !          an allocatable array with same rank
    !          as variable varname (or 1 dimension less).
    ! nslice:  optional index along dimension slicedim
    ! slicedim: optional, if nslice is set, index of which dimension to slice with
    !          nslice, default is ndims
    ! ncstart: optional, if ncstart and nccount are set, manually specify the
    !          start and count of netCDF write
    ! nccount: optional, if ncstart and nccount are set, manually specify the
    !          start and count of netCDF write
    ! errcode: optional error return code.  If not specified,
    !          program will stop if a nonzero error code returned
    !          from netcdf library.

  subroutine read_attribute(dset, attname, values, varname, errcode)
    ! read attribute 'attname' return in 'values'.  If optional
    ! argument 'varname' is given, a variable attribute is returned.
    ! if the attribute is a 1d array, values should be an allocatable 1d
    ! array of the correct type.

  subroutine write_attribute(dset, attname, values, varname, errcode)
    ! write attribute 'attname' with data in 'values'.  If optional
    ! argument 'varname' is given, a variable attribute is written.
    ! values can be a real(4), real(8), integer, string or 1d array.
```
