!=========================================================
  program grads2netcdf
!=========================================================

   use mod_netcdf
   implicit none

   character(len=50) :: gfile
   character(len=50) :: nfile
   character(len=50) :: rfile
   namelist/fileinfo/ gfile,nfile,rfile

   integer, parameter :: CLEN=1024

   ! for grads file
   integer           :: nx,ny
   character(CLEN)   :: varn
   character(CLEN)   :: vard
   character(CLEN)   :: varu  
   namelist/gridinfo/ nx,ny,varn,vard,varu
   real(4),allocatable  :: gdata(:,:)

   type(netcdf_handler) :: hd_r    ! handler
   type(netcdf_handler) :: hd_w    ! handler

    !--- general
    integer              :: count(4)  ! unit array size to write

    !--- global attribute
    character(CLEN)   :: title
    character(CLEN)   :: history
    character(CLEN)   :: comment
    
    !--- long, lati, level, time
    integer              :: imax
    integer              :: jmax
    integer              :: kmax 
    integer              :: tmax    
    real(8),allocatable  :: lon(:)
    real(8),allocatable  :: lat(:)
    real(8),allocatable  :: lev(:)
    real(8),allocatable  :: time(:)
    
    character(CLEN)       :: lon_units
    character(CLEN)       :: lat_units
    character(CLEN)       :: lev_units
    character(CLEN)       :: time_units
    
    !---variable (basic)
    character(CLEN) :: var_name  ! e.g. 'ms_tem', 'sa_u10m'
    character(CLEN) :: var_desc  ! description of the variable
    character(CLEN) :: var_units
    real(4) :: var_missing

    !--- variable (related to compression)    
    logical :: var_try_comp_2byte ! just try to
    logical :: var_comp_2byte     ! force to
    logical :: var_comp_hdf5      ! force to
    real(4) :: var_valid_min
    real(4) :: var_valid_max
    real(4) :: var_force_to_set_valid_min
    real(4) :: var_force_to_set_valid_max
    integer :: chunksizes(4)

!----------------------------------------------
!  read namelist
!----------------------------------------------
   open(10,file='conv.in',status='old')
   read(10,nml=fileinfo)
   write(6,nml=fileinfo)
   read(10,nml=gridinfo)
   write(6,nml=gridinfo)
   close(10)

!----------------------------------------------
!  read grads file
!----------------------------------------------

   write(6,*) "...open & read grads data ","'",trim(gfile),"'..."
   allocate( gdata(nx,ny) )
   open(50,file=trim(gfile),access='direct',form='unformatted', &
        status='unknown',recl=nx*ny*4)
   read(50,rec=1) gdata
   close(50)

!----------------------------------------------
!  read header info from reference file
!----------------------------------------------
    write(6,*) "...open & read header information from ","'",trim(rfile),"'..."
    call netcdf_open_for_read ( &
           hd_r,                &   ! out
           trim(rfile)          &   ! in
           )

    allocate( lon( hd_r%imax ) )
    allocate( lat( hd_r%jmax ) )
    allocate( lev( 1 ) )
    allocate( time( 1 ) )

    write(6,*) "...setting parameters..."
    title       = trim(hd_r%title)
    history     = trim(hd_r%history)
    comment     = trim(hd_r%comment)
    imax        = hd_r%imax
    jmax        = hd_r%jmax
    kmax        = 1
    tmax        = 1
    lon(1:imax) = hd_r%lon(1:imax)
    lat(1:jmax) = hd_r%lat(1:jmax)
    lev(1)      = 0
    time(1)     = 0
    lon_units   = hd_r%lon_units
    lat_units   = hd_r%lat_units
    lev_units   = hd_r%lev_units
    time_units  = 'minutes since 1999-05-01 00:00'
   
    var_name   = trim(varn)
    var_desc   = trim(vard)  ! description of the variable
    var_units  = trim(varu)
    var_missing = -0.99900E+35
    
    if((nx.ne.imax).or.(ny.ne.jmax))then
      write(6,*) "mismatch (nx,ny) = (",nx,ny,") ; (imax,jmax) = (",imax,jmax,")"
      stop
    endif

    call netcdf_close(hd_r)

!----------------------------------------------
!  write with netcdf format
!----------------------------------------------

    write(6,*) "...open file ","'",trim(nfile),"'..."

    call netcdf_open_for_write (   &
           hd_w,                   &   ! out
           trim(nfile),            &   ! in
           !count=count,            &
           title=title,            &
           history=history,        &
           comment=comment,        &
           imax=imax,              &
           jmax=jmax,              &
           kmax=kmax,              &
           tmax=tmax,              &
           lon=lon,                &
           lat=lat,                &
           lev=lev,                &
           time=time,              &
           lon_units=lon_units,    &
           lat_units=lat_units,    &
           lev_units=lev_units,    &
           time_units=time_units,  &
           var_name=var_name,      &
           var_desc=var_desc,      &
           var_units=var_units,    &
           var_missing=var_missing &
           )
   
    write(6,*) "...write data to ","'",trim(var_name),"'..."
    call netcdf_write( hd_w, gdata, k=kmax, t=tmax)

    write(6,*) "...close...","'",trim(nfile),"'..."
    call netcdf_close( hd_w )

  end program grads2netcdf
