!-------------------------------------------------------------------------------
!
!+  NetCDF module
!
!-------------------------------------------------------------------------------
module mod_netcdf
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !      This module reads/writes NetCDF3 or NetCDF4/HDF5 file for NICAM output.
  ! 
  !++ Current Corresponding Author : C.Kodama
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.90      13-04-18   C.Kodama  : [NEW]
  !
  !      -----------------------------------------------------------------------
  !
  !++ Usage for reading:
  ! 1. netcdf_open_for_read
  !  ( a. netcdf_set_count (if necessary) )
  !  ( b. netcdf_create_grads_ctl (if necessary) )
  ! 2. netcdf_read (necessary times)
  ! 3. netcdf_close
  !
  !++ Usage for writting:
  ! 1. netcdf_open_for_write
  ! 2. netcdf_write (necessary times)
  ! 3. netcdf_close
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  use netcdf  ! need NetCDF4/HDF5
!  use mod_misc, only : &
!       MISC_get_available_fid
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !
  !++ Private parameters
  !
  integer,parameter :: NOT_OPENED     = 0
  integer,parameter :: OPEN_FOR_READ  = 1
  integer,parameter :: OPEN_FOR_WRITE = -1
  integer,parameter :: CLEN           = 1024
  !
  !++ Public procedures
  !
  public :: netcdf_open_for_write
  public :: netcdf_open_for_read
  public :: netcdf_write
  public :: netcdf_read
  public :: netcdf_read_dim
  public :: netcdf_close
  public :: netcdf_create_grads_ctl

  !--- interfaces for multiple array dimensions
  interface netcdf_write
     module procedure netcdf_write_0d
     module procedure netcdf_write_1d
     module procedure netcdf_write_2d
     module procedure netcdf_write_3d
     module procedure netcdf_write_4d
  end interface

  interface netcdf_read
     module procedure netcdf_read_0d
     module procedure netcdf_read_1d
     module procedure netcdf_read_2d
     module procedure netcdf_read_3d
     module procedure netcdf_read_4d
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- Handler for each netCDF file
  !   Do not access each member outside mod_netcdf.f90.
  type,public :: netcdf_handler

     !--- general
     integer             :: id_nc
     integer             :: count(4)
     integer             :: status = NOT_OPENED

     !--- global attribute
     character(CLEN)     :: title
     character(CLEN)     :: history
     character(CLEN)     :: comment

     !--- longitude
     integer             :: id_dim_lon
     integer             :: id_var_lon
     integer             :: imax = -1
     real(8),allocatable :: lon(:)
     character(CLEN)     :: lon_units

     !--- latitude
     integer             :: id_dim_lat
     integer             :: id_var_lat
     integer             :: jmax = -1
     real(8),allocatable :: lat(:)
     character(CLEN)     :: lat_units

     !--- level
     integer             :: id_dim_lev
     integer             :: id_var_lev
     integer             :: kmax = -1
     real(8),allocatable :: lev(:)
     character(CLEN)     :: lev_units

     !--- time
     integer             :: id_dim_time
     integer             :: id_var_time
     integer             :: tmax = -1
     real(8),allocatable :: time(:)
     character(CLEN)     :: time_units

     !--- variable (basic)
     integer             :: id_var
     character(CLEN)     :: var_name
     character(CLEN)     :: var_desc
     character(CLEN)     :: var_units
     real(4)             :: var_missing

     !--- variable (related to compression)
     integer             :: nf90_type
     integer             :: nf90_cmode
     real(4)             :: var_valid_min
     real(4)             :: var_valid_max
     real(4)             :: var_scale
     real(4)             :: var_offset
     integer(2)          :: var_missing_int2
     real(4)             :: var_force_to_set_valid_min    ! maximal value to force to set var_valid_min
     real(4)             :: var_force_to_set_valid_max    ! minimal value to force to set var_valid_max
     integer             :: chunksizes(4)
     logical             :: shuffle       = .true.   ! shuffle filter (affect compression rate)
     logical             :: fletcher32    = .false.  ! checksum
     integer             :: deflate_level = 1        ! level of compression

     !--- statistics monitor (for check)
     real(4)             :: var_org_min  ! actual minimal value of the original data
     real(4)             :: var_org_max  ! actual maximal value of the original data
     real(4)             :: var_out_min  ! actual minimal value of the output data
     real(4)             :: var_out_max  ! actual maximal value of the output data

  end type netcdf_handler

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------  

  !-----------------------------------------------------------------------------
  subroutine netcdf_open_for_write( &
       nc,                         &
       ncfile,                     &
       count,                      &
       title,                      &
       history,                    &
       comment,                    &
       imax,                       &
       jmax,                       &
       kmax,                       &
       tmax,                       &
       lon,                        &
       lat,                        &
       lev,                        &
       time,                       &
       lon_units,                  &
       lat_units,                  &
       lev_units,                  &
       time_units,                 &
       var_name,                   &
       var_desc,                   &
       var_units,                  &
       var_missing,                &
       var_try_comp_2byte,         &
       var_comp_2byte,             &
       var_comp_hdf5,              &
       var_valid_min,              &
       var_valid_max,              &
       var_force_to_set_valid_min, &
       var_force_to_set_valid_max, &
       chunksizes                  &
       )
    type(netcdf_handler),intent(out)          :: nc        ! handler
    character(*),        intent( in)          :: ncfile    ! output netcdf filename
    integer,             intent( in),optional :: count(4)  ! unit array size to write
    character(*),        intent( in),optional :: title
    character(*),        intent( in),optional :: history
    character(*),        intent( in),optional :: comment
    
    integer,             intent( in),optional :: imax
    integer,             intent( in),optional :: jmax
    integer,             intent( in),optional :: kmax
    integer,             intent( in),optional :: tmax
    
    real(8),             intent( in),optional :: lon(:)
    real(8),             intent( in),optional :: lat(:)
    real(8),             intent( in),optional :: lev(:)
    real(8),             intent( in),optional :: time(:)
    
    character(*),        intent( in),optional :: lon_units
    character(*),        intent( in),optional :: lat_units
    character(*),        intent( in),optional :: lev_units
    character(*),        intent( in),optional :: time_units
    
    character(*),        intent( in),optional :: var_name  ! e.g. 'ms_tem', 'sa_u10m'
    character(*),        intent( in),optional :: var_desc  ! description of the variable
    character(*),        intent( in),optional :: var_units
    real(4),             intent( in),optional :: var_missing
    
    logical,             intent( in),optional :: var_try_comp_2byte ! just try to
    logical,             intent( in),optional :: var_comp_2byte     ! force to
    logical,             intent( in),optional :: var_comp_hdf5      ! force to
    real(4),             intent( in),optional :: var_valid_min
    real(4),             intent( in),optional :: var_valid_max
    real(4),             intent( in),optional :: var_force_to_set_valid_min
    real(4),             intent( in),optional :: var_force_to_set_valid_max
    integer,             intent( in),optional :: chunksizes(4)
    
    integer :: i, j, k, t

    !--- set mode
    if( nc%status /= NOT_OPENED ) then
       write(*,*) 'error: file is already opened.'
       stop 1
    end if
    nc%status = OPEN_FOR_WRITE

    !--- global attribute
    nc%title = 'N/A'
    if( present( title ) ) nc%title = trim( title )
    
    nc%history = 'Generated by mod_netcdf.f90.'
    if( present( history ) ) nc%history = trim( history )

    nc%comment = 'Be careful that definition of time coordinate depends on the dataset and time mode (snapshot/average).'
    if( present( comment ) ) nc%comment = trim(nc%comment) // trim( comment )

    !--- longitude
    nc%imax = 1
    if( present( imax ) ) nc%imax = imax

    allocate( nc%lon( nc%imax ) )
    nc%lon(1:nc%imax) = (/ ( real(i-1), i=1, nc%imax ) /)
    if( present( lon ) ) nc%lon(1:nc%imax) = lon(1:nc%imax)

    nc%lon_units = 'degrees_east'
    if( present( lon_units ) ) nc%lon_units = trim( lon_units )

    !--- latitude
    nc%jmax = 1
    if( present( jmax ) ) nc%jmax = jmax

    allocate( nc%lat( nc%jmax ) )
    nc%lat(1:nc%jmax) = (/ ( real(j-1), j=1, nc%jmax ) /)
    if( present( lat ) ) nc%lat(1:nc%jmax) = lat(1:nc%jmax)

    nc%lat_units = 'degrees_north'
    if( present( lat_units ) ) nc%lat_units = trim( lat_units )

    !--- level
    nc%kmax = 1
    if( present( kmax ) ) nc%kmax = kmax

    allocate( nc%lev( nc%kmax ) )
    nc%lev(1:nc%kmax) = (/ ( real(k-1), k=1, nc%kmax ) /)
    if( present( lev ) ) nc%lev(1:nc%kmax) = lev(1:nc%kmax)

    nc%lev_units = 'm'
    if( present( lev_units ) ) nc%lev_units = trim( lev_units )

    !--- time
    nc%tmax = 1
    if( present( tmax ) ) nc%tmax = tmax

    allocate( nc%time( nc%tmax ) )
    nc%time(1:nc%tmax) = (/ ( real(t-1), t=1, nc%tmax ) /)
    if( present( time ) ) nc%time(1:nc%tmax) = time(1:nc%tmax)

    nc%time_units = 'minutes since 2000-01-01 00:00'
    if( present( time_units ) ) nc%time_units = trim( time_units )

    !--- unit dimension
    nc%count(1) = nc%imax
    nc%count(2) = nc%jmax
    nc%count(3) = 1
    nc%count(4) = 1
    if( present( count ) ) nc%count(1:4) = count(1:4)

    !--- variable (basic)
    nc%var_name = 'none'
    if( present( var_name ) ) nc%var_name = trim( var_name )
    write(*,*) 'var_name     = ', trim( nc%var_name )
    
    nc%var_desc = 'N/A'
    if( present( var_desc ) ) nc%var_desc = trim( var_desc )

    nc%var_units = 'N/A'
    if( present( var_units ) ) nc%var_units = trim( var_units )

    nc%var_missing = -0.99900E+35
    if( present( var_missing ) ) nc%var_missing = var_missing

    !--- variable (related to compression)

    ! default table
    call var_comp_def_table( nc, trim(nc%var_name) )

    ! NetCDF3 or NetCDF4/HDF5
    !nc%nf90_cmode = NF90_HDF5
    if( present( var_comp_hdf5 ) ) then
       if( var_comp_hdf5 ) then
          nc%nf90_cmode = NF90_HDF5
       else
          nc%nf90_cmode = NF90_NOCLOBBER
       end if
    end if
    if( nc%nf90_cmode == NF90_HDF5      ) write(*,*) '  using NetCDF4/HDF5'
    if( nc%nf90_cmode == NF90_NOCLOBBER ) write(*,*) '  using NetCDF3'

    ! valid range
    !nc%var_valid_min = nc%var_missing
    if( present( var_valid_min ) ) then
       if( var_valid_min /= nc%var_missing ) then
          nc%var_valid_min = var_valid_min
       end if
    end if

    !nc%var_valid_max = nc%var_missing
    if( present( var_valid_max ) ) then
       if( var_valid_max /= nc%var_missing ) then
          nc%var_valid_max = var_valid_max
       end if
    end if

    ! whether to force to limit values within valid_min/max for high compression rate
    !nc%var_force_to_set_valid_min = nc%var_missing
    if( present( var_force_to_set_valid_min ) ) then
       if( var_force_to_set_valid_min /= nc%var_missing ) then
          nc%var_force_to_set_valid_min = var_force_to_set_valid_min
       end if
    end if

    !nc%var_force_to_set_valid_max = nc%var_missing
    if( present( var_force_to_set_valid_max ) ) then
       if( var_force_to_set_valid_max /= nc%var_missing ) then
          nc%var_force_to_set_valid_max = var_force_to_set_valid_max
       end if
    end if

    if( nc%var_force_to_set_valid_min /= nc%var_missing ) then
       write(*,*) '  forced to set ', nc%var_valid_min, &
            '(for ', trim(nc%var_name) , ' <', nc%var_force_to_set_valid_min, ')'
    end if
    
    if( nc%var_force_to_set_valid_max /= nc%var_missing ) then
       write(*,*) '  forced to set ', nc%var_valid_max, &
            '(for ', trim(nc%var_name) , ' >', nc%var_force_to_set_valid_max, ')'
    end if

    ! 4-byte float or 2-byte integer
    if( nc%var_valid_min == nc%var_missing &
         .or. nc%var_valid_max == nc%var_missing ) then
          nc%nf90_type = NF90_FLOAT
    end if

    if( present( var_try_comp_2byte ) ) then
       if( .not. var_try_comp_2byte ) then
          nc%nf90_type = NF90_FLOAT
       end if
    end if

    if( present( var_comp_2byte ) ) then
       if( .not. var_comp_2byte ) then
          nc%nf90_type = NF90_FLOAT
       else
          nc%nf90_type = NF90_SHORT
       end if
    end if

    nc%var_scale  = nc%var_missing
    nc%var_offset = nc%var_missing
    if( nc%nf90_type == NF90_SHORT ) then
       write(*,*) '  type = 2-byte integer'

       if( nc%var_valid_min /= nc%var_missing &
            .and. nc%var_valid_max /= nc%var_missing ) then
          nc%var_missing_int2 = -2**15
          nc%var_scale = ( nc%var_valid_max - nc%var_valid_min ) / (2**15)
          nc%var_offset = ( nc%var_valid_max + nc%var_valid_min ) * 0.5
          write(*,*) '    scale_factor = ', nc%var_scale, &
                     '(possible err =', nc%var_scale/2, ')'
          write(*,*) '    add_offset   = ', nc%var_offset
       else
          write(*,*) 'error: Both var_valid_min/max should be set when var_comp_2byte=.true.'
          stop 1
       end if
    else
       write(*,*) '  type = 4-byte float'
    end if

    ! chunksize
    nc%chunksizes(:) = (/ nc%imax, nc%jmax, 1, 1 /)
    if( present( chunksizes ) ) nc%chunksizes(:) = chunksizes(:)


    !---------- start netCDF define mode

    !--- create a new netCDF and enter define mode
    ! define mode: to prepare/write metadata
    ! id_dim_*: define as dimension
    ! id_var_*: define as variable
    call check( nf90_create( trim(ncfile), nc%nf90_cmode, nc%id_nc ) )

    !--- longitude
    call check( nf90_def_dim( &
         nc%id_nc, 'lon', nc%imax, nc%id_dim_lon ) )
    call check( nf90_def_var( &
         nc%id_nc, 'lon', NF90_DOUBLE, (/ nc%id_dim_lon /), nc%id_var_lon ) )

    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_lon, 'units', trim(nc%lon_units) ) )
    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_lon, 'long_name', 'Longitude' ) )

    !--- latitude
    call check( nf90_def_dim( &
         nc%id_nc, 'lat', nc%jmax, nc%id_dim_lat ) )
    call check( nf90_def_var( &
         nc%id_nc, 'lat', NF90_DOUBLE, (/ nc%id_dim_lat /), nc%id_var_lat ) )

    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_lat, 'units', trim(nc%lat_units) ) )
    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_lat, 'long_name', 'Latitude' ) )
 
    !--- level
    call check( nf90_def_dim( &
         nc%id_nc, 'lev', nc%kmax, nc%id_dim_lev ) )
    call check( nf90_def_var( &
         nc%id_nc, 'lev', NF90_DOUBLE, (/ nc%id_dim_lev /), nc%id_var_lev ) )

    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_lev, 'units', trim(nc%lev_units) ) )
    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_lev, 'long_name', 'Level' ) )

    !--- time
    call check( nf90_def_dim( &
         nc%id_nc, 'time', nc%tmax, nc%id_dim_time ) )
    call check( nf90_def_var( &
         nc%id_nc, 'time', NF90_DOUBLE, (/ nc%id_dim_time /), nc%id_var_time ) )

    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_time, 'units', trim(nc%time_units) ) )
    call check( nf90_put_att( &
         nc%id_nc, nc%id_var_time, 'long_name', 'Time' ) )

    !--- variable
    if( nc%nf90_cmode == NF90_HDF5 ) then
       call check( nf90_def_var( &
            nc%id_nc, trim(nc%var_name), nc%nf90_type, &
            (/ nc%id_dim_lon, nc%id_dim_lat, nc%id_dim_lev, nc%id_dim_time /), &
            nc%id_var, &
            chunksizes    = nc%chunksizes, &
            shuffle       = nc%shuffle, &
            fletcher32    = nc%fletcher32, &
            deflate_level = nc%deflate_level ) )
    else
       call check( nf90_def_var( &
            nc%id_nc, trim(nc%var_name), nc%nf90_type, &
            (/ nc%id_dim_lon, nc%id_dim_lat, nc%id_dim_lev, nc%id_dim_time /), &
            nc%id_var ) )
    end if

    call check( nf90_put_att( &
         nc%id_nc, nc%id_var, "units", trim(nc%var_units) ) )
    call check( nf90_put_att( &
         nc%id_nc, nc%id_var, "long_name", trim(nc%var_desc) ) )

    if( nc%nf90_type == NF90_SHORT ) then
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "_FillValue", nc%var_missing_int2 ) )
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "scale_factor", nc%var_scale ) )
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "add_offset", nc%var_offset ) )
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "valid_min", int(-2**15+1,2) ) )
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "valid_max", int(2**15-1,2) ) )
    else
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "_FillValue", nc%var_missing ) )
       if( nc%var_valid_min /= nc%var_missing ) then
          call check( nf90_put_att( &
               nc%id_nc, nc%id_var, "valid_min", nc%var_valid_min ) )
       end if
       if( nc%var_valid_max /= nc%var_missing ) then
          call check( nf90_put_att( &
               nc%id_nc, nc%id_var, "valid_max", nc%var_valid_max ) )
       end if
    end if

    if( nc%var_force_to_set_valid_min /= nc%var_missing ) then
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "var_force_to_set_valid_min", nc%var_force_to_set_valid_min ) )
    end if
    if( nc%var_force_to_set_valid_max /= nc%var_missing ) then
       call check( nf90_put_att( &
            nc%id_nc, nc%id_var, "var_force_to_set_valid_max", nc%var_force_to_set_valid_max ) )
    end if

    !--- global attribute
    call check( nf90_put_att( &
         nc%id_nc, NF90_GLOBAL, "title", nc%title ) )
    call check( nf90_put_att( &
         nc%id_nc, NF90_GLOBAL, "history", trim(nc%history) ) )
    call check( nf90_put_att( &
         nc%id_nc, NF90_GLOBAL, "comment", trim(nc%comment) ) )

    !--- write metadata to the disk and leave define mode
    call check( nf90_enddef( nc%id_nc ) )

    !---------- write coordinate values

    !--- longitude
    call check( nf90_put_var( &
         nc%id_nc, nc%id_var_lon, nc%lon, start=(/ 1 /), count=(/ nc%imax /) ) )

    !--- latitude
    call check( nf90_put_var( &
         nc%id_nc, nc%id_var_lat, nc%lat, start=(/ 1 /), count=(/ nc%jmax /) ) )

    !--- level
    call check( nf90_put_var( &
         nc%id_nc, nc%id_var_lev, nc%lev, start=(/ 1 /), count=(/ nc%kmax /) ) )

    !--- time
    call check( nf90_put_var( &
         nc%id_nc, nc%id_var_time, nc%time, start=(/ 1 /), count=(/ nc%tmax /) ) )

    !--- statistics monitor
    nc%var_org_min = nc%var_missing
    nc%var_org_max = nc%var_missing
    nc%var_out_min = nc%var_missing
    nc%var_out_max = nc%var_missing

  end subroutine netcdf_open_for_write


  subroutine netcdf_open_for_read( &
       nc,         &
       ncfile,     &
       count,      &
       title,      &
       history,    &
       comment,    &
       imax,       &
       jmax,       &
       kmax,       &
       tmax,       &
!       lon,        &
!       lat,        &
!       lev,        &
!       time,       &
       lon_units,  &
       lat_units,  &
       lev_units,  &
       time_units, &
       var_name,   &
       var_desc,   &
       var_units,  &
       var_missing &
       )
    type(netcdf_handler),intent(out)          :: nc
    character(*),        intent( in)          :: ncfile
    integer,             intent(out),optional :: count(4)
    character(*),        intent(out),optional :: title
    character(*),        intent(out),optional :: history
    character(*),        intent(out),optional :: comment

    integer,             intent(out),optional :: imax
    integer,             intent(out),optional :: jmax
    integer,             intent(out),optional :: kmax
    integer,             intent(out),optional :: tmax

!    real(8),allocatable, intent(out),optional :: lon(:) ! strictly, illegal!
!    real(8),allocatable, intent(out),optional :: lat(:)
!    real(8),allocatable, intent(out),optional :: lev(:)
!    real(8),allocatable, intent(out),optional :: time(:)

    character(*),        intent(out),optional :: lon_units
    character(*),        intent(out),optional :: lat_units
    character(*),        intent(out),optional :: lev_units
    character(*),        intent(out),optional :: time_units

    character(*),        intent(out),optional :: var_name
    character(*),        intent(out),optional :: var_desc
    character(*),        intent(out),optional :: var_units
    real(4),             intent(out),optional :: var_missing

    integer :: i
    integer             :: wrk_nvar
    integer,allocatable :: wrk_ids(:)
    character(CLEN)     :: wrk_name
    integer :: wrk_type
    logical :: wrk_lon_flag
    logical :: wrk_lat_flag
    logical :: wrk_lev_flag
    logical :: wrk_time_flag
    logical :: wrk_var_flag

    !--- set mode
    if( nc%status /= NOT_OPENED ) then
       write(*,*) 'error: file is already opened.'
       stop 1
    end if
    nc%status = OPEN_FOR_READ

    !--- open
    !
    call check( nf90_open( &
         path=trim(ncfile), mode=NF90_NOWRITE, ncid=nc%id_nc ) )

    !--- get all the variable IDs
    call check( nf90_inquire( nc%id_nc, nvariables=wrk_nvar ) )
    if( wrk_nvar /= 5 ) then
       write(*,*) 'error: wrk_nvar is not equal to 5 (=', wrk_nvar, ')'
       write(*,*) 'It seems that ncfile=', trim(ncfile), &
            'has more than one field variables or some dimension(s) does not exist. ', &
            'mod_netcdf.f90 does not support such a file.'
       stop 1
    end if

    allocate( wrk_ids( wrk_nvar ) )
    call check( nf90_inq_varids( nc%id_nc, wrk_nvar, wrk_ids ) )
    wrk_lon_flag  = .false.
    wrk_lat_flag  = .false.
    wrk_lev_flag  = .false.
    wrk_time_flag = .false.
    wrk_var_flag  = .false.
    do i=1, wrk_nvar
       call check( nf90_inquire_variable( &
            nc%id_nc, wrk_ids(i), name=wrk_name, xtype=wrk_type ) )

       if( trim(wrk_name) == 'lon' .and. .not. wrk_lon_flag ) then
          nc%id_var_lon = wrk_ids(i)
          wrk_lon_flag = .true.
       else if( trim(wrk_name) == 'lat' .and. .not. wrk_lat_flag ) then
          nc%id_var_lat = wrk_ids(i)
          wrk_lat_flag = .true.
       else if( trim(wrk_name) == 'lev' .and. .not. wrk_lev_flag ) then
          nc%id_var_lev = wrk_ids(i)
          wrk_lev_flag = .true.
       else if( trim(wrk_name) == 'time' .and. .not. wrk_time_flag ) then
          nc%id_var_time = wrk_ids(i)
          wrk_time_flag = .true.
       else if( .not. wrk_var_flag ) then
          nc%id_var = wrk_ids(i)
          wrk_var_flag = .true.
          nc%nf90_type = wrk_type
          nc%var_name = trim( wrk_name )
       else
          write(*,*) 'error in reading variables.'
          stop 1
       end if
    end do
    deallocate( wrk_ids )

    !--- global attribute
    call check( nf90_get_att( &
         nc%id_nc, NF90_GLOBAL, "title", nc%title ) )
    if( present( title ) ) title = nc%title

    call check( nf90_get_att( &
         nc%id_nc, NF90_GLOBAL, "history", nc%history ) )
    if( present( history ) ) history = nc%history

    call check( nf90_get_att( &
         nc%id_nc, NF90_GLOBAL, "comment", nc%comment ) )
    if( present( comment ) ) comment = nc%comment

    !--- longitude
    call check( nf90_inq_dimid( &
         nc%id_nc, "lon", nc%id_dim_lon ) )
    call check( nf90_inquire_dimension( &
         nc%id_nc, nc%id_dim_lon, len=nc%imax ) )
    if( present( imax ) ) imax = nc%imax

    allocate( nc%lon( nc%imax ) )
    call check( nf90_get_var( &
         nc%id_nc, nc%id_var_lon, start=(/ 1 /), count=(/ nc%imax /), values=nc%lon ) )
!    if( present( lon ) ) then
!       if( allocated( lon ) ) deallocate( lon )
!       allocate( lon( nc%imax ) )
!       lon(1:nc%imax) = nc%lon(1:nc%imax)
!    end if

    call check( nf90_get_att( &
         nc%id_nc, nc%id_var_lon, 'units', nc%lon_units ) )
    if( present( lon_units ) ) lon_units = trim( nc%lon_units )

    !--- latitude
    call check( nf90_inq_dimid( &
         nc%id_nc, "lat", nc%id_dim_lat ) )
    call check( nf90_inquire_dimension( &
         nc%id_nc, nc%id_dim_lat, len=nc%jmax ) )
    if( present( jmax ) ) jmax = nc%jmax

    allocate( nc%lat( nc%jmax ) )
    call check( nf90_get_var( &
         nc%id_nc, nc%id_var_lat, start=(/ 1 /), count=(/ nc%jmax /), values=nc%lat ) )
!    if( present( lat ) ) then
!       if( allocated( lat ) ) deallocate( lat )
!       allocate( lat( nc%jmax ) )
!       lat(1:nc%jmax) = nc%lat(1:nc%jmax)
!    end if

    call check( nf90_get_att( &
         nc%id_nc, nc%id_var_lat, 'units', nc%lat_units ) )
    if( present( lat_units ) ) lat_units = trim( nc%lat_units )

    !--- level
    call check( nf90_inq_dimid( &
         nc%id_nc, "lev", nc%id_dim_lev ) )
    call check( nf90_inquire_dimension( &
         nc%id_nc, nc%id_dim_lev, len=nc%kmax ) )
    if( present( kmax ) ) kmax = nc%kmax

    allocate( nc%lev( nc%kmax ) )
    call check( nf90_get_var( &
         nc%id_nc, nc%id_var_lev, start=(/ 1 /), count=(/ nc%kmax /), values=nc%lev ) )
!    if( present( lev ) ) then
!       if( allocated( lev ) ) deallocate( lev )
!       allocate( lev( nc%kmax ) )
!       lev(1:nc%kmax) = nc%lev(1:nc%kmax)
!    end if

    call check( nf90_get_att( &
         nc%id_nc, nc%id_var_lev, 'units', nc%lev_units ) )
    if( present( lev_units ) ) lev_units = trim( nc%lev_units )

    !--- time
    call check( nf90_inq_dimid( &
         nc%id_nc, "time", nc%id_dim_time ) )
    call check( nf90_inquire_dimension( &
         nc%id_nc, nc%id_dim_time, len=nc%tmax ) )
    if( present( tmax ) ) tmax = nc%tmax

    allocate( nc%time( nc%tmax ) )
    call check( nf90_get_var( &
         nc%id_nc, nc%id_var_time, start=(/ 1 /), count=(/ nc%tmax /), values=nc%time ) )
!    if( present( time ) ) then
!       if( allocated( time ) ) deallocate( time )
!       allocate( time( nc%tmax ) )
!       time(1:nc%tmax) = nc%time(1:nc%tmax)
!    end if

    call check( nf90_get_att( &
         nc%id_nc, nc%id_var_time, 'units', nc%time_units ) )
    if( present( time_units ) ) time_units = trim( nc%time_units )

    !--- count
    nc%count(1) = nc%imax
    nc%count(2) = nc%jmax
    nc%count(3) = 1
    nc%count(4) = 1
    if( present(count) ) count(1:4) = nc%count(1:4)

    !--- variable (basic)
    if( present( var_name ) ) var_name = trim( nc%var_name )
    !write(*,*) 'var_name     = ', trim( nc%var_name )

    call check( nf90_get_att( &
         nc%id_nc, nc%id_var, 'long_name', nc%var_desc ) )
    if( present( var_desc ) ) var_desc = trim( nc%var_desc )

    call check( nf90_get_att( &
         nc%id_nc, nc%id_var, 'units', nc%var_units ) )
    if( present( var_units ) ) var_units = trim( nc%var_units )

    !--- variable (related to compression)
    if( nc%nf90_type == NF90_SHORT ) then
       call check( nf90_get_att( &
            nc%id_nc, nc%id_var, "_FillValue", nc%var_missing_int2 ) )
       call check( nf90_get_att( &
            nc%id_nc, nc%id_var, "scale_factor", nc%var_scale ) )
       call check( nf90_get_att( &
            nc%id_nc, nc%id_var, "add_offset", nc%var_offset ) )
       nc%var_missing = -0.99900E+35  ! default
    else
       call check( nf90_get_att( &
            nc%id_nc, nc%id_var, "_FillValue", nc%var_missing ) )
    end if

    if( present( var_missing ) ) var_missing = nc%var_missing

  end subroutine netcdf_open_for_read


  subroutine netcdf_set_count( nc, count )
    type(netcdf_handler),intent(inout) :: nc
    integer,             intent(in)    :: count(4)

    !--- check mode
    if( nc%status /= OPEN_FOR_READ ) then
       write(*,*) 'error: file is not opened for reading.'
       stop 1
    end if

    !--- overwrite nc%count
    nc%count(:) = count(:)
  end subroutine netcdf_set_count


  subroutine netcdf_write_0d( nc, var, i, j, k, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var
    integer,             intent(   in) :: i, j, k, t
    real(4),allocatable                :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( 1, 1, 1, 1 ) )
    wrk_var_real4(1,1,1,1) = var
    call netcdf_write_main( nc, wrk_var_real4, i=1, j=j, k=k, t=t )
    deallocate( wrk_var_real4 )
    !call netcdf_write_main( nc, spread(spread(spread(spread(var,1,1),2,1),3,1),4,1), i=1, j=1, k=k, t=t )
  end subroutine netcdf_write_0d

  subroutine netcdf_write_1d( nc, var, j, k, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:)
    integer,             intent(   in) :: j, k, t
    real(4),allocatable                :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( &
         lbound(var,1):ubound(var,1), &
         1, 1, 1 ) )
    wrk_var_real4(:,1,1,1) = var(:)
    call netcdf_write_main( nc, wrk_var_real4, j=j, k=k, t=t )
    deallocate( wrk_var_real4 )
    !call netcdf_write_main( nc, spread(spread(spread(var,2,1),3,1),4,1), j=1, k=k, t=t )
  end subroutine netcdf_write_1d

  subroutine netcdf_write_2d( nc, var, k, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:,:)
    integer,             intent(   in) :: k, t
    real(4),allocatable                :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( &
         lbound(var,1):ubound(var,1), &
         lbound(var,2):ubound(var,2),&
         1, 1 ) )
    wrk_var_real4(:,:,1,1) = var(:,:)
    call netcdf_write_main( nc, wrk_var_real4, k=k, t=t )
    deallocate( wrk_var_real4 )
    !call netcdf_write_main( nc, spread(spread(var,3,1),4,1), k=k, t=t )
  end subroutine netcdf_write_2d

  subroutine netcdf_write_3d( nc, var, t )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:,:,:)
    integer,             intent(   in) :: t
    real(4),allocatable                :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( &
         lbound(var,1):ubound(var,1), &
         lbound(var,2):ubound(var,2),&
         lbound(var,3):ubound(var,3),&
         1 ) )
    wrk_var_real4(:,:,:,1) = var(:,:,:)
    call netcdf_write_main( nc, wrk_var_real4, t=t )
    deallocate( wrk_var_real4 )
    !call netcdf_write_main( nc, spread(var,4,1), t=t )
  end subroutine netcdf_write_3d

  subroutine netcdf_write_4d( nc, var )
    type(netcdf_handler),intent(inout) :: nc
    real(4),             intent(   in) :: var(:,:,:,:)
    real(4),allocatable                :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( &
         lbound(var,1):ubound(var,1), &
         lbound(var,2):ubound(var,2),&
         lbound(var,3):ubound(var,3),&
         lbound(var,4):ubound(var,4) ) )
    wrk_var_real4(:,:,:,:) = var(:,:,:,:)
    call netcdf_write_main( nc, wrk_var_real4 )
    deallocate( wrk_var_real4 )
    !call netcdf_write_main( nc, var )
  end subroutine netcdf_write_4d

  subroutine netcdf_write_main( nc, wrk_var_real4, i, j, k, t )
    type(netcdf_handler),intent(inout)          :: nc
    real(4),             intent(inout)          :: wrk_var_real4(:,:,:,:)
    integer,             intent(   in),optional :: i, j, k, t  ! record position to write

    integer(2),allocatable :: wrk_var_int2(:,:,:,:)

    integer    :: start(4)
    integer    :: out_of_range = 0
    real(4)    :: wrk_min, wrk_max
    integer    :: p, q, r, s

    start(1:4) = (/ 1, 1, 1, 1 /)

    !--- check mode
    if( nc%status /= OPEN_FOR_WRITE ) then
       write(*,*) 'error: file is not opened for writting.'
       stop 1
    end if

    !--- set record position
    if( present(i) ) then
       if( i > 1 .and. nc%count(1) > 1 ) stop 'error: both i>1 and nc%count(1)>1.'
       start(1) = i
    end if
    if( present(j) ) then
       if( j > 1 .and. nc%count(2) > 1 ) stop 'error: both j>1 and nc%count(2)>1.'
       start(2) = j
    end if
    if( present(k) ) then
       if( k > 1 .and. nc%count(3) > 1 ) stop 'error: both k>1 and nc%count(3)>1.'
       start(3) = k
    end if
    if( present(t) ) then
       if( t > 1 .and. nc%count(4) > 1 ) stop 'error: both t>1 and nc%count(4)>1.'
       start(4) = t
    end if

    !--- monitor for input data
    if( any( ubound(wrk_var_real4)-lbound(wrk_var_real4)+1 /= nc%count ) ) then
       write(*,*) 'error: array shapes differ.'
       write(*,'(A,4I8)') 'netcdf_open_for_write(), count = ', nc%count
       write(*,'(A,4I8)') 'netcdf_write_main(), var shape = ', ubound(wrk_var_real4)-lbound(wrk_var_real4)+1
       stop 1
    end if
    wrk_min = minval( wrk_var_real4(:,:,:,:), mask=(wrk_var_real4(:,:,:,:)/=nc%var_missing) )
    if( nc%var_org_min == nc%var_missing ) then
       nc%var_org_min = wrk_min
    else
       nc%var_org_min = min( nc%var_org_min, wrk_min )
    end if
       
    wrk_max = maxval( wrk_var_real4(:,:,:,:), mask=(wrk_var_real4(:,:,:,:)/=nc%var_missing) )
    if( nc%var_org_max == nc%var_missing ) then
       nc%var_org_max = wrk_max
    else
       nc%var_org_max = max( nc%var_org_max, wrk_max )
    end if

    !--- force to limit value range if necessary
    if( nc%var_force_to_set_valid_min /= nc%var_missing ) then
       where( wrk_var_real4 /= nc%var_missing &
            .and. wrk_var_real4 < nc%var_force_to_set_valid_min )
          wrk_var_real4 = nc%var_valid_min
       end where
    end if

    if( nc%var_force_to_set_valid_max /= nc%var_missing ) then
       where( wrk_var_real4 /= nc%var_missing &
            .and. wrk_var_real4 > nc%var_force_to_set_valid_max )
          wrk_var_real4 = nc%var_valid_max
       end where
    end if

    !--- monitor for output data
    wrk_min = minval( wrk_var_real4(:,:,:,:), mask=(wrk_var_real4(:,:,:,:)/=nc%var_missing) )
    if( nc%var_out_min == nc%var_missing ) then
       nc%var_out_min = wrk_min
    else
       nc%var_out_min = min( nc%var_out_min, wrk_min )
    end if
       
    wrk_max = maxval( wrk_var_real4(:,:,:,:), mask=(wrk_var_real4(:,:,:,:)/=nc%var_missing) )
    if( nc%var_out_max == nc%var_missing ) then
       nc%var_out_max = wrk_max
    else
       nc%var_out_max = max( nc%var_out_max, wrk_max )
    end if

    if( nc%var_valid_min /= nc%var_missing ) then
       out_of_range = count( &
            wrk_var_real4(:,:,:,:) /= nc%var_missing &
            .and. wrk_var_real4(:,:,:,:) < nc%var_valid_min )
    end if
    if( nc%var_valid_max /= nc%var_missing ) then
       out_of_range = out_of_range + count( &
            wrk_var_real4(:,:,:,:) /= nc%var_missing &
            .and. wrk_var_real4(:,:,:,:) > nc%var_valid_max )
    end if

    if( out_of_range > 0 ) then
       write(*,*) 'error: some output values are below var_valid_min (#=', out_of_range, ')'

       write(*,'(A)',advance='no') ' valid range  : ['
       if( nc%var_valid_min /= nc%var_missing ) then
          write(*,'(E15.7)',advance='no') nc%var_valid_min
       else
          write(*,'(A15)',advance='no') ' '
       end if
       write(*,'(A)',advance='no') ':'
       if( nc%var_valid_max /= nc%var_missing ) then
          write(*,'(E15.7)',advance='no') nc%var_valid_max
       else
          write(*,'(A15)',advance='no') ' '
       end if
       write(*,'(A)') ' ]'

       write(*,'(A)',advance='no') ' actual range : ['
       write(*,'(E15.7)',advance='no') minval( wrk_var_real4(:,:,:,:), mask=(wrk_var_real4(:,:,:,:)/=nc%var_missing) )
       write(*,'(A)',advance='no') ':'
          write(*,'(E15.7)',advance='no') maxval( wrk_var_real4(:,:,:,:), mask=(wrk_var_real4(:,:,:,:)/=nc%var_missing) )
       write(*,'(A)') ' ]'

       write(*,'(A,4I8)') '(i,j,k,t)   :', start(:)
       stop 1
    end if

    !--- write
    if( nc%nf90_type == NF90_FLOAT ) then ! as 4-byte float
       call check( nf90_put_var( &
            nc%id_nc, nc%id_var, wrk_var_real4, &
            start=start, count=nc%count ) )
    else                                  ! as 2-byte integer
       allocate( wrk_var_int2(  nc%count(1), nc%count(2), nc%count(3), nc%count(4) ) )

!      NOTE by C.Kodama: "where" statement does not correctly work with huge imax and jmax.
!      Probably because of the need of stack memory.
!       where( wrk_var_real4(1:nc%count(1),1:nc%count(2),1:nc%count(3),1:nc%count(4)) == nc%var_missing  )
!          wrk_var_int2 = nc%var_missing_int2
!       elsewhere
!          wrk_var_int2 = nint( ( wrk_var_real4 - nc%var_offset ) / nc%var_scale )
!       end where
       do s=1, nc%count(4)
          do r=1, nc%count(3)
             do q=1, nc%count(2)
                do p=1, nc%count(1)

                   if( wrk_var_real4(p,q,r,s) /= nc%var_missing ) then
                      wrk_var_int2(p,q,r,s) = nint( ( wrk_var_real4(p,q,r,s) - nc%var_offset ) / nc%var_scale )
                   else
                      wrk_var_int2(p,q,r,s) = nc%var_missing_int2
                   end if
                   
                end do
             end do
          end do
       end do

       call check( nf90_put_var( &
            nc%id_nc, nc%id_var, wrk_var_int2, &
            start=start, count=nc%count ) )

       deallocate( wrk_var_int2 )
    end if

  end subroutine netcdf_write_main


  subroutine netcdf_read_0d( nc, var, i, j, k, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var
    integer,             intent(in)  :: i, j, k, t
    real(4),allocatable              :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( 1, 1, 1, 1 ) )
    call netcdf_read_main( nc, wrk_var_real4, i=i, j=j, k=k, t=t )
    var = wrk_var_real4(1,1,1,1)
    deallocate( wrk_var_real4 )
  end subroutine netcdf_read_0d

  subroutine netcdf_read_1d( nc, var, j, k, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:)
    integer,             intent(in)  :: j, k, t
    real(4),allocatable              :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( &
         lbound(var,1):ubound(var,1), &
         1, 1, 1 ) )
    call netcdf_read_main( nc, wrk_var_real4, j=j, k=k, t=t )
    var(:) = wrk_var_real4(:,1,1,1)
    deallocate( wrk_var_real4 )
  end subroutine netcdf_read_1d

  subroutine netcdf_read_2d( nc, var, k, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:,:)
    integer,             intent(in)  :: k, t
    real(4),allocatable              :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( &
         lbound(var,1):ubound(var,1), &
         lbound(var,2):ubound(var,2),&
         1, 1 ) )
    call netcdf_read_main( nc, wrk_var_real4, k=k, t=t )
    var(:,:) = wrk_var_real4(:,:,1,1)
    deallocate( wrk_var_real4 )
  end subroutine netcdf_read_2d

  subroutine netcdf_read_3d( nc, var, t )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:,:,:)
    integer,             intent(in)  :: t
    real(4),allocatable              :: wrk_var_real4(:,:,:,:)

    allocate( wrk_var_real4( &
         lbound(var,1):ubound(var,1), &
         lbound(var,2):ubound(var,2),&
         lbound(var,3):ubound(var,3),&
         1 ) )
    call netcdf_read_main( nc, wrk_var_real4, t=t )
    var(:,:,:) = wrk_var_real4(:,:,:,1)
    deallocate( wrk_var_real4 )
  end subroutine netcdf_read_3d

  subroutine netcdf_read_4d( nc, var )
    type(netcdf_handler),intent(in)  :: nc
    real(4),             intent(out) :: var(:,:,:,:)
    call netcdf_read_main( nc, var )
  end subroutine netcdf_read_4d

  subroutine netcdf_read_main( nc, wrk_var_real4, i, j, k, t )
    type(netcdf_handler),intent( in)          :: nc
    real(4),             intent(out)          :: wrk_var_real4(:,:,:,:)
    integer,             intent( in),optional :: i, j, k, t  ! record position to read

    integer(2),allocatable :: wrk_var_int2(:,:,:,:)

    integer :: start(4)
    integer :: p, q, r, s

    start(1:4) = (/ 1, 1, 1, 1 /)

    !--- check mode
    if( nc%status /= OPEN_FOR_READ ) then
       write(*,*) 'error: file is not opened for reading.'
       stop 1
    end if

    !--- set record position
    if( present(i) ) then
       if( i > 1 .and. nc%count(1) > 1 ) stop 'error: both i>1 and nc%count(1)>1.'
       start(1) = i
    end if
    if( present(j) ) then
       if( j > 1 .and. nc%count(2) > 1 ) stop 'error: both j>1 and nc%count(2)>1.'
       start(2) = j
    end if
    if( present(k) ) then
       if( k > 1 .and. nc%count(3) > 1 ) stop 'error: both k>1 and nc%count(3)>1.'
       start(3) = k
    end if
    if( present(t) ) then
       if( t > 1 .and. nc%count(4) > 1 ) stop 'error: both t>1 and nc%count(4)>1.'
       start(4) = t
    end if

    !--- read
    if( nc%nf90_type == NF90_FLOAT ) then ! as 4-byte float
       call check( nf90_get_var( &
            nc%id_nc, nc%id_var, wrk_var_real4, &
            start=start, count=nc%count ) )

    else                                  ! as 2-byte integer
       allocate( wrk_var_int2(  nc%count(1), nc%count(2), nc%count(3), nc%count(4) ) )

       call check( nf90_get_var( &
            nc%id_nc, nc%id_var, wrk_var_int2, &
            start=start, count=nc%count ) )

!       where( wrk_var_int2 == nc%var_missing_int2 )
!          var = nc%var_missing
!       elsewhere
!          var = wrk_var_int2 * nc%var_scale + nc%var_offset
!       end where
       do s=1, nc%count(4)
          do r=1, nc%count(3)
             do q=1, nc%count(2)
                do p=1, nc%count(1)

                   if( wrk_var_int2(p,q,r,s) /= nc%var_missing_int2 ) then
                      wrk_var_real4(p,q,r,s) = wrk_var_int2(p,q,r,s) * nc%var_scale + nc%var_offset
                   else
                      wrk_var_real4(p,q,r,s) = nc%var_missing
                   end if

                end do
             end do
          end do
       end do

       deallocate( wrk_var_int2 )    
    end if

  end subroutine netcdf_read_main


  subroutine netcdf_read_dim( &
       nc,         &
       lon,        &
       lat,        &
       lev,        &
       time        &
       )
    type(netcdf_handler),intent(in) :: nc
    real(8),intent(out),optional    :: lon(1:nc%imax)
    real(8),intent(out),optional    :: lat(1:nc%jmax)
    real(8),intent(out),optional    :: lev(1:nc%kmax)
    real(8),intent(out),optional    :: time(1:nc%tmax)

    !--- set mode
    if( nc%status /= OPEN_FOR_READ ) then
       write(*,*) 'error: file is not opened.'
       stop 1
    end if

    if( present( lon ) ) then
       lon(1:nc%imax) = nc%lon(1:nc%imax)
    end if
    if( present( lat ) ) then
       lat(1:nc%jmax) = nc%lat(1:nc%jmax)
    end if
    if( present( lev ) ) then
       lev(1:nc%kmax) = nc%lev(1:nc%kmax)
    end if
    if( present( time ) ) then
       time(1:nc%tmax) = nc%time(1:nc%tmax)
    end if
  end subroutine netcdf_read_dim
  

  subroutine netcdf_close( nc )
    type(netcdf_handler),intent(inout) :: nc

    if( nc%status == OPEN_FOR_WRITE ) then
       write(*,*) '  statistics monitor:'
       write(*,*) '    original data range = [', nc%var_org_min, ':', nc%var_org_max, ']'
       write(*,*) '    output   data range = [', nc%var_out_min, ':', nc%var_out_max, ']'
    end if
    
    call check( nf90_close( nc%id_nc ) )
    
    if( allocated( nc%lon  ) ) deallocate( nc%lon  )
    if( allocated( nc%lat  ) ) deallocate( nc%lat  )
    if( allocated( nc%lev  ) ) deallocate( nc%lev  )
    if( allocated( nc%time ) ) deallocate( nc%time )
    nc%status = NOT_OPENED
  end subroutine netcdf_close


  ! based on mod_grads.f90 and prg_ico2ll.f90 in NICAM
  subroutine netcdf_create_grads_ctl( nc, fid_ctl, endian )
    type(netcdf_handler),intent( in)          :: nc
    integer,             intent( in)          :: fid_ctl
    character(*),        intent( in),optional :: endian

    real(8)         :: lon_start, lon_int
    real(8)         :: lon_int_min, lon_int_max

    real(8)         :: lat_start, lat_int
    real(8)         :: lat_int_min, lat_int_max

    real(8)         :: lev_start, lev_int
    real(8)         :: lev_int_min, lev_int_max

    character(CLEN) :: time_start, time_int
    integer         :: time_int_min, time_int_max
    character(CLEN) :: wrk(4)
    character(CLEN) :: c_year, c_mon, c_day
    integer         :: mon
    character(3)    :: c_mon_ref(12) = (/ &
         'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
         'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

    integer :: i, j, k, t
    character(CLEN) :: wrk_endian

    !--- check mode
    if( nc%status == NOT_OPENED ) then
       write(*,*) 'error: file is not opened.'
       stop 1
    end if

    !--- endian
    wrk_endian = 'BIG_ENDIAN'
    if( present( endian ) ) wrk_endian = trim( endian )

    !--- longitude
    lon_start   = nc%var_missing
    lon_int     = nc%var_missing
    lon_int_min = minval( nc%lon(2:nc%imax) - nc%lon(1:nc%imax-1) )
    lon_int_max = maxval( nc%lon(2:nc%imax) - nc%lon(1:nc%imax-1) )
    if( lon_int_min == lon_int_max ) then
       lon_start = nc%lon(1)
       lon_int   = lon_int_min
    end if

    !--- latitude
    lat_start   = nc%var_missing
    lat_int     = nc%var_missing
    lat_int_min = minval( nc%lat(2:nc%jmax) - nc%lat(1:nc%jmax-1) )
    lat_int_max = maxval( nc%lat(2:nc%jmax) - nc%lat(1:nc%jmax-1) )
    if( lat_int_min == lat_int_max ) then
       lat_start = nc%lat(1)
       lat_int   = lat_int_min
    end if

    !--- level
    lev_start   = nc%var_missing
    lev_int     = nc%var_missing
    lev_int_min = minval( nc%lev(2:nc%kmax) - nc%lev(1:nc%kmax-1) )
    lev_int_max = maxval( nc%lev(2:nc%kmax) - nc%lev(1:nc%kmax-1) )
    if( lev_int_min == lev_int_max ) then
       lev_start = nc%lev(1)
       lev_int   = lev_int_min
    end if

    !--- time
    read( nc%time_units, * ) wrk(1:4)

    c_year = wrk(3)
    c_year = c_year(1:4)

    c_mon = wrk(3)
    c_mon = c_mon(6:7)
    read( c_mon, * ) mon
    c_mon = trim( c_mon_ref(mon) )

    c_day = wrk(3)
    c_day = c_day(9:10)

    time_start = trim(wrk(4)(1:5)) // 'Z' &
         // trim(c_day) // trim(c_mon) // trim(c_year)

    time_int_min = int( minval( nc%time(2:nc%tmax) - nc%time(1:nc%tmax-1) ) )
    time_int_max = int( maxval( nc%time(2:nc%tmax) - nc%time(1:nc%tmax-1) ) )
    if( time_int_min /= time_int_max ) then
       write(*,*) 'error: time intervals depend on time.'
       stop 1
    end if

    write( time_int, * ) time_int_min
    time_int = adjustl( time_int )

    if( trim(wrk(1)) == 'minutes' ) then
       time_int = trim(time_int) // 'mn'
    else if( trim(wrk(1)) == 'hours' ) then
       time_int = trim(time_int) // 'hr'
    else if( trim(wrk(1)) == 'days' ) then
       time_int = trim(time_int) // 'dy'
    else if( trim(wrk(1)) == 'months' ) then
       time_int = trim(time_int) // 'mn'
    else if( trim(wrk(1)) == 'years' ) then
       time_int = trim(time_int) // 'yr'
    else
       write(*,*) 'error: time increment = ', trim(wrk(1)), 'is not supported.'
       stop 1
    end if

    !--- output
    write(fid_ctl,fmt='(2a)')        'DSET ','^'//trim(nc%var_name)//'.grd'
    write(fid_ctl,fmt='(2a)')        'TITLE ',trim(nc%title)
    write(fid_ctl,fmt='(2a)')        'OPTIONS ' // trim(wrk_endian)
    write(fid_ctl,fmt='(a,e12.5)')   'UNDEF ',real(nc%var_missing,4)

    if( lon_start == nc%var_missing .or. lon_int == nc%var_missing ) then
       write(fid_ctl,fmt='(a,i5,a)')    'XDEF ',nc%imax, ' LEVELS'
       !write(fid_ctl,fmt='(5x,5f10.3)')    (lon(i)*180.0D0/pi,i=1,imax)
       !write(fid_ctl,fmt='(5x,5f10.3)') (nc%lon(i),i=1,nc%imax)
       write(fid_ctl,fmt='(5x,5f10.4)') (nc%lon(i),i=1,nc%imax)  ! [mod] C.Kodama: follow prg_ico2ll.f90.
    else
       !write(fid_ctl,fmt='(a,i5,a,2f10.3)')    'XDEF ',nc%imax, ' LINEAR', &
       write(fid_ctl,fmt='(a,i5,a,2f10.4)')    'XDEF ',nc%imax, ' LINEAR', &
            lon_start, lon_int
    end if

    if( lat_start == nc%var_missing .or. lat_int == nc%var_missing ) then
       write(fid_ctl,fmt='(a,i5,a)')    'YDEF ',nc%jmax, ' LEVELS'
       !write(fid_ctl,fmt='(5x,5f10.3)')(lat(j)*180.0D0/pi,j=1,jmax)
       !write(fid_ctl,fmt='(5x,5f10.3)') (nc%lat(j),j=1,nc%jmax)
       write(fid_ctl,fmt='(5x,5f10.4)') (nc%lat(j),j=1,nc%jmax)
    else
       !write(fid_ctl,fmt='(a,i5,a,2f10.3)')  'YDEF ',nc%jmax, ' LINEAR', &
       write(fid_ctl,fmt='(a,i5,a,2f10.4)')  'YDEF ',nc%jmax, ' LINEAR', &
            lat_start, lat_int
    end if

    if( lev_start == nc%var_missing .or. lev_int == nc%var_missing ) then
!       write(fid_ctl,fmt='(a,i5,a)')    'ZDEF ',nc%kmax, ' LEVELS'
!       write(fid_ctl,fmt='(5f10.3)')((nc%lev(k)),k=1,nc%kmax)
!       write(fid_ctl,fmt='(a,i5,a,5f10.3)')    'ZDEF ',nc%kmax, ' LEVELS', ((nc%lev(k)),k=1,nc%kmax)
       write(fid_ctl,fmt='(a,i5,a)',advance='no')    'ZDEF ',nc%kmax, ' LEVELS'
       !write(fid_ctl,fmt='(5f10.3)') ((nc%lev(k)),k=1,nc%kmax)
       write(fid_ctl,fmt='(5f11.4)') ((nc%lev(k)),k=1,nc%kmax)  ! [mod] C.Kodama: to support dfq_isccp2.nc.
    else
       !write(fid_ctl,fmt='(a,i5,a,2f10.3)')    'ZDEF ',nc%kmax, ' LINEAR', &
       write(fid_ctl,fmt='(a,i5,a,2f11.4)')    'ZDEF ',nc%kmax, ' LINEAR', &
            lev_start, lev_int
    end if

    write(fid_ctl,fmt='(a,i5,2a,1x,a)')  'TDEF ',nc%tmax,' LINEAR ', &
         trim(time_start),trim(time_int)
    write(fid_ctl,fmt='(a,i5)')     'VARS ',1
    if( nc%kmax == 1 ) then
       write(fid_ctl,fmt='(a,2i5,1x,4a)')    trim(nc%var_name),0, 99,  &
            trim(nc%var_desc),'[',trim(nc%var_units),']'
    else
       write(fid_ctl,fmt='(a,2i5,1x,4a)')    trim(nc%var_name),nc%kmax, 99,  &
            trim(nc%var_desc),'[',trim(nc%var_units),']'
    end if
    write(fid_ctl,fmt='(a)')        'ENDVARS '

  end subroutine netcdf_create_grads_ctl


  subroutine check( status )
    implicit none
    integer,intent( in) :: status
    
    if( status == NF90_NOERR ) return
    
    write(*,*) 'Netcdf error number = ', status
    write(*,*) 'Error Message: ', trim( NF90_STRERROR(status) )
    stop 1
  end subroutine check
  
  subroutine handle_error( status )
    implicit none
    integer, intent (in) :: status
    
    if( status == nf90_NoErr ) return
    
    write(*,*) 'Netcdf error number = ', status
    write(*,*) 'Error Message: ', trim( NF90_STRERROR(status) )
    stop 1
  end subroutine handle_error


  subroutine var_comp_def_table( nc, var_name )
    type(netcdf_handler),intent(inout) :: nc
    character(*),        intent(   in) :: var_name

    !--- default of default: loss-less compression
    nc%nf90_type                  = NF90_FLOAT
    nc%nf90_cmode                 = NF90_HDF5
    nc%var_valid_min              = nc%var_missing
    nc%var_valid_max              = nc%var_missing
    nc%var_force_to_set_valid_min = nc%var_missing
    nc%var_force_to_set_valid_max = nc%var_missing

    !--- default for each variable
    select case( trim(var_name) )
    case('dfq_isccp2')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 0.0
       nc%var_valid_max              = 1.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('la_lai')
       continue

    case('la_rof')
       continue

    case('la_rofl')
       continue

    case('la_snw')
       continue

    case('la_tg')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              =  50.0
       nc%var_valid_max              = 350.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('la_wg')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 0.0
       nc%var_valid_max              = 1.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ms_dh')
       continue

    case('ms_pres')
       continue

    case('ms_qc')
       nc%var_valid_min              = 0.0
       nc%var_force_to_set_valid_min = nc%var_valid_min

    case('ms_qg')
       nc%var_valid_min              = 0.0
       nc%var_force_to_set_valid_min = 1.0e-12

    case('ms_qi')
       nc%var_valid_min              = 0.0
       nc%var_force_to_set_valid_min = nc%var_valid_min

    case('ms_qr')
       nc%var_valid_min              = 0.0
       nc%var_force_to_set_valid_min = 1.0e-12

    case('ms_qs')
       nc%var_valid_min              = 0.0
       nc%var_force_to_set_valid_min = 1.0e-12

    case('ms_qv')
       nc%var_valid_min              = 0.0
       nc%var_force_to_set_valid_min = nc%var_valid_min

    case('ms_rh')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 0.0
       nc%var_valid_max              = 3.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ms_rho')
       continue

    case('ms_tem')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 50.0
       nc%var_valid_max              = 700.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ms_u')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = -300.0
       nc%var_valid_max              = 300.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ms_v')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = -300.0
       nc%var_valid_max              = 300.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ms_w')
       continue

    case('oa_ice')
       continue

    case('oa_icr')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 0.0
       nc%var_valid_max              = 1.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('oa_ist')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              =  50.0
       nc%var_valid_max              = 300.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('oa_snw')
       continue

    case('oa_sst')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 250.0
       nc%var_valid_max              = 350.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_albedo')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 0.0
       nc%var_valid_max              = 1.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_cld_frac')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = 0.0
       nc%var_valid_max              = 1.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_cldi')
       continue

    case('ss_cldi')
       continue

    case('sa_cldw')
       continue

    case('ss_cldw')
       continue

    case('sa_evap')
       continue

    case('ss_evap')
       continue

    case('sa_evap_energy')
       continue

    case('sa_lh_sfc')
       continue

    case('sa_lwd_sfc')
       continue

    case('sa_lwd_toa')
       continue

    case('sa_lwu_sfc')
       continue

    case('sa_lwu_toa')
       continue

    case('ss_lwu_toa')
       continue

    case('sa_lwu_toa_c')
       continue

    case('sa_q2m')
       continue

    case('sa_sh_sfc')
       continue

    case('sa_slp')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              =  70000.0
       nc%var_valid_max              = 130000.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ss_slp')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              =  70000.0
       nc%var_valid_max              = 130000.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_swd_sfc')
       continue

    case('sa_swd_toa')
       continue

    case('sa_swu_sfc')
       continue

    case('sa_swu_toa')
       continue

    case('sa_swu_toa_c')
       continue

    case('sa_t2m')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              =  50.0
       nc%var_valid_max              = 700.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_tauu')
       continue

    case('sa_tauv')
       continue

    case('sa_tem_atm')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              =  50.0
       nc%var_valid_max              = 700.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_tem_sfc')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              =  50.0
       nc%var_valid_max              = 700.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_tppn')
       continue

    case('ss_tppn')
       continue

    case('sa_energy')
       continue

    case('sa_u10m')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = -300.0
       nc%var_valid_max              =  300.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ss_u10m')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = -300.0
       nc%var_valid_max              =  300.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_v10m')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = -300.0
       nc%var_valid_max              =  300.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('ss_v10m')
       nc%nf90_type                  = NF90_SHORT
       nc%var_valid_min              = -300.0
       nc%var_valid_max              =  300.0
       nc%var_force_to_set_valid_min = nc%var_valid_min
       nc%var_force_to_set_valid_max = nc%var_valid_max

    case('sa_vap_atm')
       continue

!    case default
!       write(*,*) 'warning:'
    end select

  end subroutine var_comp_def_table

end module mod_netcdf
