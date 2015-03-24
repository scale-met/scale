program netcdf2grads_h
  !-----------------------------------------------------------------------------
  !> post-process for scale high-performance
  !> convert from netcdf to grads format (with combine and slice)
  !>
  !> original: 2015/02/03
  !-----------------------------------------------------------------------------

  use netcdf
  use mpi

  implicit none

  !--- parameters
  integer,  parameter :: CSHT       = 16
  integer,  parameter :: CMID       = 64
  integer,  parameter :: CLNG       = 128
  integer,  parameter :: SP         = 4
  integer,  parameter :: DP         = 8
  integer,  parameter :: max_vcount = 500
  integer,  parameter :: max_tcount = 500
  integer,  parameter :: max_zcount = 100
  integer,  parameter :: master     = 0
  integer,  parameter :: FID_STD    = 6
  integer,  parameter :: FID_CONF   = 20
  integer,  parameter :: FID_RCNF   = 21
  integer,  parameter :: FID_LOGF   = 22
  integer,  parameter :: FID_CTL    = 23
  integer,  parameter :: FID_DAT    = 24
  integer             :: FID_LOG    = 22

  integer,  parameter :: err_internal = 0
  integer,  parameter :: err_netcdf   = -1
  real(SP), parameter :: UNDEF_SP     = -9.9999D7

  integer,  parameter :: vt_2d     = 0  ! vtype index
  integer,  parameter :: vt_3d     = 1  ! vtype index
  integer,  parameter :: vt_height = 2  ! vtype index
  integer,  parameter :: vt_land   = 3  ! vtype index
  integer,  parameter :: vt_urban  = 4  ! vtype index
  integer,  parameter :: vt_tpmsk  = 5  ! vtype index

  integer,  parameter :: a_slice   = 0  ! atype index
  integer,  parameter :: a_max     = 1  ! atype index
  integer,  parameter :: a_min     = 2  ! atype index
  integer,  parameter :: a_sum     = 3  ! atype index
  integer,  parameter :: a_ave     = 4  ! atype index

  !--- setting for scale running
  integer         :: START_TSTEP    = 1
  integer         :: END_TSTEP      = 1
  integer         :: INC_TSTEP      = 1
  integer         :: DOMAIN_NUM     = 1
  integer         :: VCOUNT         = 1
  integer         :: ZCOUNT         = 1
  integer         :: ZSTART         = 3
  integer         :: TARGET_ZLEV(max_zcount) = 3
  real(DP)        :: EXTRA_TINTERVAL = -9.999
  character(5)    :: EXTRA_TUNIT    = ""
  character(CLNG) :: IDIR           = "./data"
  character(CLNG) :: ODIR           = "."
  character(CLNG) :: CONFFILE       = "./run.conf"
  character(CSHT) :: VNAME(max_vcount) = ""
  character(CSHT) :: ANALYSIS       = "SLICE"
  character(CSHT) :: Z_LEV_TYPE     = "GRID"
  character(5)    :: DELT           = "1mn"
  character(15)   :: STIME          = "00:00Z01JAN2000"
  character(15)   :: FTIME          = "2000010100"
  character(CLNG) :: LOG_BASENAME   = "LOG"
  logical         :: LOG_ALL_OUTPUT = .false.
  logical         :: Z_LEV_LIST     = .true.
  logical         :: Z_MERGE_OUT    = .false.  ! only for slice

  integer         :: PRC_NUM_X
  integer         :: PRC_NUM_Y
  integer         :: TIME_STARTDATE(6)
  real(DP)        :: HISTORY_DEFAULT_TINTERVAL
  character(CMID) :: HISTORY_DEFAULT_BASENAME
  character(5)    :: HISTORY_DEFAULT_TUNIT
  logical         :: HISTORY_DEFAULT_ZINTERP = .true.
  logical         :: HIST_BND = .true.

  !--- variables for work
  real(SP),allocatable :: cz(:), cdz(:)
  real(SP),allocatable :: cx(:), cdx(:)
  real(SP),allocatable :: cy(:), cdy(:)
  real(SP),allocatable :: var_2d(:,:)
  real(SP),allocatable :: p_var(:,:)

  real(DP),allocatable :: p_cx(:), p_cdx(:)
  real(DP),allocatable :: p_cy(:), p_cdy(:)
  real(DP),allocatable :: p_3d(:,:,:,:)
  real(DP),allocatable :: p_3d_urban(:,:,:,:)
  real(DP),allocatable :: p_3d_land(:,:,:,:)
  real(DP),allocatable :: p_2d(:,:,:)
  real(DP),allocatable :: p_2dt(:,:)

  real(SP),allocatable :: sendbuf(:,:)
  real(SP),allocatable :: sendbuf_gx(:)
  real(SP),allocatable :: sendbuf_gy(:)
  real(SP),allocatable :: recvbuf(:,:)
  real(SP),allocatable :: cx_gather(:)
  real(SP),allocatable :: cy_gather(:)
  real(SP),allocatable :: cdx_gather(:)
  real(SP),allocatable :: cdy_gather(:)
  real(SP),allocatable :: vgrid(:)

  integer :: start_3d(4)
  integer :: start_2d(3)
  integer :: start_2dt(2)
  integer :: count_3d(4)
  integer :: count_2d(3)
  integer :: count_urban(4)
  integer :: count_land(4)
  integer :: count_height(3)
  integer :: count_tpmsk(2)

  integer :: vtype
  integer :: atype
  integer :: nt
  integer :: nx   ! num of x-dimension in the combined file
  integer :: ny   ! num of y-dimension in the combined file
  integer :: nz   ! num of z-dimension in the combined file
  integer :: uz   ! num of z-dimension in the combined file for urban
  integer :: lz   ! num of z-dimension in the combined file for land
  integer :: nzg  ! num of z-dimension in the combined file for grid
  integer :: nxh  ! num of halo grids for x-dimension
  integer :: nyh  ! num of halo grids for y-dimension
  integer :: nzh  ! num of halo grids for z-dimension
  integer :: ks, ke
  integer :: nst, nen

  integer :: xproc, yproc, tproc
  integer :: it, ix, jy, iz, iv
  integer :: zz
  integer :: ndim
  integer :: idom
  integer :: irank
  integer :: isize
  integer :: ierr

  integer :: nxp   ! num of x-dimension in partial file
  integer :: nyp   ! num of y-dimension in partial file
  integer :: nxgp  ! num of x-dimension in partial file for grid
  integer :: nygp  ! num of y-dimension in partial file for grid
  integer :: mnxp  ! maximum num of x-grid in partial file among all files
  integer :: mnyp  ! maximum num of y-grid in partial file among all files
  integer :: nxg_tproc
  integer :: nyg_tproc
  integer :: mnx   ! maximum num of x-dimension in partial file among all files
  integer :: mny   ! maximum num of y-dimension in partial file among all files

  integer :: nm, nmnge
  integer :: work
  integer :: yy, mm, dd, hh, mn, sc
  integer,allocatable :: rk_mnge(:)

  character(6)    :: num
  character(3)    :: cmm(12)
  character(CLNG) :: ncfile
  character(CMID) :: varname
  character(CMID) :: fconf
  character(CLNG) :: fname_save

  logical :: flag_bnd       = .false.
  logical :: open_file      = .false.
  logical :: LOUT           = .false.
  !-----------------------------------------------------------------------------

  data cmm / "JAN", "FEB", "MAR", "APL", "MAY", "JUN", &
             "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" /

  namelist /LOGOUT/            &
    LOG_BASENAME,              &
    LOG_ALL_OUTPUT
  namelist /INFO/              &
    START_TSTEP,               &
    END_TSTEP,                 &
    INC_TSTEP,                 &
    DOMAIN_NUM,                &
    VCOUNT,                    &
    ZCOUNT,                    &
    ZSTART,                    &
    ANALYSIS,                  &
    CONFFILE,                  &
    IDIR,                      &
    ODIR,                      &
    EXTRA_TINTERVAL,           &
    EXTRA_TUNIT,               &
    Z_LEV_LIST,                &
    Z_LEV_TYPE,                &
    Z_MERGE_OUT
  namelist /VARI/              &
    VNAME,                     &
    TARGET_ZLEV
  namelist /GRADS/             &
    DELT,                      &
    STIME
  namelist  / PARAM_TIME /     &
    TIME_STARTDATE
  namelist  / PARAM_PRC /      &
    PRC_NUM_X,                 &
    PRC_NUM_Y
  namelist  / PARAM_HISTORY /  &
    HISTORY_DEFAULT_BASENAME,  &
    HISTORY_DEFAULT_TINTERVAL, &
    HISTORY_DEFAULT_TUNIT,     &
    HISTORY_DEFAULT_ZINTERP
  namelist  /PARAM_HIST/       &
    HIST_BND
  !-----------------------------------------------------------------------------

  !### initialization
  call MPI_INIT( ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, isize, ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, irank, ierr )

  !--- Read from argument
  if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
     write(*, *) "ERROR: Program needs a config file!"
     call err_abort( 1, __LINE__ )
  else
     call get_command_argument(1,fconf)
  endif

  call read_conf_logout
  call logio_init( irank, LOUT )

  if ( LOUT ) write( FID_LOG, '(1X,A)') "-----------------------------------"
  if ( LOUT ) write( FID_LOG, '(1X,A)') "      START netcdf to grads"
  if ( LOUT ) write( FID_LOG, '(1X,A)') "-----------------------------------"
  if ( LOUT ) write( FID_LOG, '(1X,"isize:",I6,2X,"irank:",I6)') isize, irank
  call read_conf

  idom  = DOMAIN_NUM
  xproc = PRC_NUM_X
  yproc = PRC_NUM_Y
  tproc = xproc * yproc
  nmnge = tproc / isize
  work = mod( tproc, isize )
  if ( work /= 0 ) then
     if ( LOUT ) write (*, *) "ERROR: specified num of mpi processes is not adequate."
     if ( LOUT ) write (*, *) "*** specify the num; PRC_X*PRC_Y shuold be divisable by the num."
     call err_abort( 1, __LINE__ )
  elseif ( isize > tproc ) then
     if ( LOUT ) write (*, *) "ERROR: num of mpi processes is larger than that of the scale-les run."
     call err_abort( 1, __LINE__ )
  endif

  allocate ( rk_mnge (nmnge) )
  call set_rank_manage( irank, nmnge, rk_mnge )

  !### read and combine
  write ( num,'(I6.6)' ) rk_mnge(1)
  ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
  if ( LOUT ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File: ", trim(ncfile)

  call netcdf_retrieve_dims( ncfile )
  call set_array_size
  call allocation( irank ) ! "irank" is right, not "rk_mnge"

  do nm = 1, nmnge
     call netcdf_read_grid( rk_mnge(nm), nm )
  enddo

  call gather_grid

  if ( irank == master ) call combine_grid

  if ( Z_MERGE_OUT ) call make_vgrid

  !--- read data and combine
  nst = START_TSTEP
  nen = END_TSTEP
  nt  = (nen - nst + 1) / INC_TSTEP
  if ( nt > max_tcount ) then
     if ( LOUT ) write (*, *) "ERROR: overflow maximum of tcount"
     call err_abort( 1, __LINE__ )
  endif

  call cal_init

  do it = nst, nen, INC_TSTEP !--- time loop

     call set_calender
     if ( LOUT ) write( FID_LOG, '(1X,A,I3)' ) "+++ TIME STEP:   ", it
     if ( LOUT ) write( FID_LOG, * ) "+++ Date: ", STIME

     do iv = 1, vcount !--- var loop
        varname = trim( VNAME(iv) )
        call set_vtype( ncfile, varname, vtype, ndim )
        call set_atype
        if ( LOUT ) write( FID_LOG, '(1X,A,A,A,I1,A)' ) &
        "+++ VARIABLE: ", trim(varname), " (vtype = ", vtype, ")"

        select case( vtype )
        case ( vt_urban, vt_land, vt_height, vt_3d )

           do iz = 1, ZCOUNT        !--- level loop
              if ( atype == a_slice ) then
                 zz = TARGET_ZLEV(iz)
                 if ( LOUT ) write( FID_LOG, '(1X,A,I3)' ) "+++ Z LEVEL: ", zz
                 call check_targ_zlev( zz )
              endif

              do nm = 1, nmnge
                 call netcdf_read_var( rk_mnge(nm), nm, it, zz, varname, vtype, p_var )
              enddo

              call gather_vars( p_var, recvbuf )

              if ( irank == master ) then
                 call combine_vars_2d( recvbuf, var_2d )
                 if ( iz /= 1 .and. Z_MERGE_OUT ) then
                    call write_vars_zmerge( var_2d, iz )
                 else
                    call create_ctl( varname, idom, it, zz )
                    call write_vars( var_2d, varname, idom, it, zz )
                 endif
              endif
           enddo !--- level loop

        case ( vt_2d, vt_tpmsk )

           !if ( atype == a_slice ) then
           !   zz = TARGET_ZLEV(1)
           !   call check_targ_zlev( zz )
           !endif

           do nm = 1, nmnge
              call netcdf_read_var( rk_mnge(nm), nm, it, zz, varname, vtype, p_var )
           enddo

           call gather_vars( p_var, recvbuf )

           if ( irank == master ) then
              call combine_vars_2d( recvbuf, var_2d )
              call create_ctl( varname, idom, it, zz )
              call write_vars( var_2d, varname, idom, it, zz )
           endif

        case default
           call err_abort( 0, __LINE__ )
        end select

     enddo !--- var loop

     call cal_increment
  enddo !--- time loop

  ! finalization
  if ( open_file .and. LOUT ) close ( FID_LOG )
  call MPI_FINALIZE( ierr )

  !-----------------------------------------------------------------------------
  !> END MAIN ROUTINE
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Individual Procedures
  !-----------------------------------------------------------------------------

  !> open logfile
  !---------------------------------------------------------------------------
  subroutine logio_init( &
      irank,   & ! [in   ]
      LOUT     ) ! [inout]
    implicit none

    integer, intent(in)    :: irank
    logical, intent(inout) :: LOUT

    character(6) :: num
    !---------------------------------------------------------------------------

    LOUT = .false.
    if ( LOG_ALL_OUTPUT ) then
       LOUT = .true.
    else
       if ( irank == master ) LOUT = .true.
    endif

    if ( trim(LOG_BASENAME) .eq. "STDOUT" ) then
       FID_LOG = FID_STD
       open_file = .false.
    else
       FID_LOG = FID_LOGF
       open_file = .true.
    endif

    if ( open_file .and. LOUT ) then
       write( num,'(I6.6)' ) irank
       open ( FID_LOG, file=trim(LOG_BASENAME)//".pe"//num, &
              status='replace', form='formatted' )
    endif

    return
  end subroutine logio_init


  !> read from netcdf: retrieve dimension size
  !---------------------------------------------------------------------------
  subroutine netcdf_retrieve_dims( &
      ncfile  ) ! [in]
    implicit none

    character(CLNG), intent(in) :: ncfile

    integer :: ncid, dimid
    integer :: istat
    !---------------------------------------------------------------------------

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_dimid ( ncid,  'x',  dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nxp )
    istat = nf90_inq_dimid ( ncid,  'CX', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nxgp )
    nxh = ( nxgp - nxp ) / 2

    istat = nf90_inq_dimid ( ncid,  'y',  dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nyp )
    istat = nf90_inq_dimid ( ncid,  'CY', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nygp )
    nyh = ( nygp - nyp ) / 2

    istat = nf90_inq_dimid ( ncid,  'z',  dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nz )
    istat = nf90_inq_dimid ( ncid,  'CZ', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=nzg )
    nzh = ( nzg - nz ) / 2
    ks = 1
    ke = nz

    istat = nf90_inq_dimid ( ncid,  'uz', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=uz )

    istat = nf90_inq_dimid ( ncid,  'lz', dimid )
    istat = nf90_inquire_dimension( ncid, dimid, len=lz )

    istat = nf90_close(ncid)
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
    return
  end subroutine netcdf_retrieve_dims


  !> read from netcdf: grid data
  !---------------------------------------------------------------------------
  subroutine netcdf_read_grid( &
      imnge,   & ! [in]
      nm       ) ! [in]
    implicit none

    integer, intent(in) :: imnge
    integer, intent(in) :: nm

    integer :: ncid, varid
    integer :: istat
    integer :: is, ie, js, je
    character(CLNG) :: ncfile
    character(6)    :: num
    !---------------------------------------------------------------------------

    write ( num,'(I6.6)' ) imnge
    ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
    if ( LOUT ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File (grd): ", trim(ncfile)

    call set_index_readbuf_grid( nm, is, ie, js, je )

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_varid( ncid, 'CZ', varid )
    istat = nf90_get_var( ncid, varid, cz )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CDZ', varid )
    istat = nf90_get_var( ncid, varid, cdz )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CX', varid )
    istat = nf90_get_var( ncid, varid, p_cx(is:ie) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CDX', varid )
    istat = nf90_get_var( ncid, varid, p_cdx(is:ie) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CY', varid )
    istat = nf90_get_var( ncid, varid, p_cy(js:je) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_inq_varid( ncid, 'CDY', varid )
    istat = nf90_get_var( ncid, varid, p_cdy(js:je) )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_read_grid


  !> read from netcdf: variables data
  !---------------------------------------------------------------------------
  subroutine netcdf_read_var( &
      imnge,    & ! [in ]
      nm,       & ! [in ]
      it,       & ! [in ]
      zz,       & ! [in ]
      varname,  & ! [in ]
      vtype,    & ! [in ]
      p_var     ) ! [out]
    implicit none

    integer, intent(in) :: imnge
    integer, intent(in) :: nm, it, zz
    character(CMID), intent(in)  :: varname
    integer,         intent(in)  :: vtype
    real(SP),        intent(out) :: p_var(:,:)

    integer :: ncid, varid
    integer :: istat
    integer :: is, ie, js, je
    integer :: isn, ien, jsn, jen, nzn
    character(CLNG) :: ncfile
    character(6)    :: num
    !---------------------------------------------------------------------------

    write ( num,'(I6.6)' ) imnge
    ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
    if ( LOUT ) write( FID_LOG, '(1x,A,A)' ) "+++ Target File (var): ", trim(ncfile)

    call irank2ixjy( imnge, ix, jy )
    call set_flag_bnd( ix, jy, flag_bnd )
    call set_index_readbuf( nm, is, ie, js, je, im_bnd=flag_bnd )
    call set_index_netcdf( atype, vtype, ix, jy, it, zz, isn, ien, jsn, jen, nzn )

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_varid( ncid, trim(varname), varid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    select case( vtype )
    case ( vt_urban )
       istat = nf90_get_var( ncid, varid, p_3d_urban(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_urban )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call simple_analysis( is,  ie,  js,  je,  &
                             isn, jsn, nzn, p_3d_urban, p_var )

    case ( vt_land )
       istat = nf90_get_var( ncid, varid, p_3d_land(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_land )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call simple_analysis( is,  ie,  js,  je,  &
                             isn, jsn, nzn, p_3d_land, p_var )

    case ( vt_height )
       istat = nf90_get_var( ncid, varid, p_2d(isn:ien,jsn:jen,1), start=start_2d, count=count_height )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2d(isn:ien,jsn:jen,1) )

    case ( vt_tpmsk )
       istat = nf90_get_var( ncid, varid, p_2dt(isn:ien,jsn:jen), start=start_2dt, count=count_tpmsk )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2dt(isn:ien,jsn:jen) )

    case ( vt_3d )
       istat = nf90_get_var( ncid, varid, p_3d(isn:ien,jsn:jen,1:nzn,1), start=start_3d, count=count_3d )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       call simple_analysis( is,  ie,  js,  je,  &
                             isn, jsn, nzn, p_3d, p_var )

    case ( vt_2d )
       istat = nf90_get_var( ncid, varid, p_2d(isn:ien,jsn:jen,1), start=start_2d, count=count_2d )
       if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
       p_var(is:ie,js:je) = real( p_2d(isn:ien,jsn:jen,1) )

    case default
       call err_abort( 0, __LINE__ )

    end select

    istat = nf90_close(ncid)
    return
  end subroutine netcdf_read_var


  !> combine region divided data: [2D]
  !---------------------------------------------------------------------------
  subroutine simple_analysis( &
      is,       & ! [in ]
      ie,       & ! [in ]
      js,       & ! [in ]
      je,       & ! [in ]
      isn,      & ! [in ]
      jsn,      & ! [in ]
      nzn,      & ! [in ]
      indata,   & ! [in ]
      outdata   ) ! [out]
    implicit none

    integer,  intent(in)  :: is, ie, js, je
    integer,  intent(in)  :: isn, jsn, nzn
    real(DP), intent(in)  :: indata(:,:,:,:)
    real(SP), intent(out) :: outdata(:,:)

    integer :: i, j, ni, nj, k
    real(SP) :: work
    !---------------------------------------------------------------------------

    select case( atype )
    case ( a_slice )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          outdata(i,j) = real( indata(ni,nj,1,1) )
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_max )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          outdata(i,j) = real( maxval(indata(ni,nj,:,1)) )
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_min )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          outdata(i,j) = real( minval(indata(ni,nj,:,1)) )
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_sum )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          work = 0.0D0
          do k = 1, nzn
             work = work + real( indata(ni,nj,k,1) )
          enddo
          outdata(i,j) = work
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case ( a_ave )
       nj = jsn
       do j = js, je
       ni = isn
       do i = is, ie
          work = 0.0D0
          do k = 1, nzn
             work = work + real( indata(ni,nj,k,1) )
          enddo
          outdata(i,j) = work / real(nzn)
       ni = ni + 1
       enddo
       nj = nj + 1
       enddo
    case default
       call err_abort( 0, __LINE__ )
    end select

    return
  end subroutine simple_analysis


  !> combine region divided data: grid data
  !---------------------------------------------------------------------------
  subroutine combine_grid()
    implicit none

    integer :: ix, jy, iix, jjy
    integer :: is, ie, js, je
    integer :: gix, gjy
    integer :: gis, gie, gjs, gje
    !---------------------------------------------------------------------------

    jy = 1
    do ix = 1, xproc
       call set_index     ( ix, jy, is,  ie,  js,  je  )
       call set_index_grid( ix, jy, gis, gie, gjs, gje )

       gix = gis
       do iix = is, ie
          cdx(iix) = cdx_gather(gix+nxh)
          cx(iix)  = cx_gather( gix+nxh)
          gix = gix + 1
       enddo
    enddo

    ix = 1
    do jy = 1, yproc
       call set_index     ( ix, jy, is,  ie,  js,  je  )
       call set_index_grid( ix, jy, gis, gie, gjs, gje )

       gjy = gjs
       do jjy = js, je
          cdy(jjy) = cdy_gather(gjy+nyh)
          cy(jjy)  = cy_gather( gjy+nyh)
          gjy = gjy + 1
       enddo
    enddo

    return
  end subroutine combine_grid

  !> combine region divided data: [2D]
  !---------------------------------------------------------------------------
  subroutine combine_vars_2d( &
      gathered,   & ! [in ]
      combined    ) ! [out]
    implicit none

    real(SP), intent(in)  :: gathered(:,:)
    real(SP), intent(out) :: combined(:,:)

    integer :: ix, jy, iix, jjy, gix, gjy
    integer :: is, ie, js, je
    integer :: gis, gie, gjs, gje
    !---------------------------------------------------------------------------

    do jy = 1, yproc
    do ix = 1, xproc
       call set_index( ix, jy, is, ie, js, je )
       call set_index_gathered( ix, jy, gis, gie, gjs, gje )

       gjy = gjs
       do jjy = js, je
          gix = gis
          do iix = is, ie
             combined(iix,jjy) = gathered( gix, gjy )
             gix = gix + 1
          enddo
          gjy = gjy + 1
       enddo
    enddo
    enddo

    return
  end subroutine combine_vars_2d


  !> communication: gather [1D]
  !---------------------------------------------------------------------------
  subroutine gather_grid()
    implicit none

    integer :: sendcounts
    integer :: recvcounts
    integer :: iix, jjy
    integer :: ierr
    !---------------------------------------------------------------------------

    ! grids for x-direction
    sendcounts = nxgp*nmnge
    recvcounts = nxgp*nmnge

    sendbuf_gx(:) = UNDEF_SP
    do iix = 1, nxgp*nmnge
       sendbuf_gx(iix) = real( p_cx(iix) )
    enddo
    call MPI_GATHER( sendbuf_gx(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cx_gather(:),   &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    sendbuf_gx(:) = UNDEF_SP
    do iix = 1, nxgp*nmnge
       sendbuf_gx(iix) = real( p_cdx(iix) )
    enddo
    call MPI_GATHER( sendbuf_gx(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cdx_gather(:),  &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    ! grids for y-direction
    sendcounts = nygp*nmnge
    recvcounts = nygp*nmnge

    sendbuf_gy(:) = UNDEF_SP
    do jjy = 1, nygp*nmnge
       sendbuf_gy(jjy) = real( p_cy(jjy) )
    enddo
    call MPI_GATHER( sendbuf_gy(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cy_gather(:),   &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    sendbuf_gy(:) = UNDEF_SP
    do jjy = 1, nygp*nmnge
       sendbuf_gy(jjy) = real( p_cdy(jjy) )
    enddo
    call MPI_GATHER( sendbuf_gy(:),  &
                     sendcounts,     &
                     MPI_REAL,       &
                     cdy_gather(:),  &
                     recvcounts,     &
                     MPI_REAL,       &
                     master,         &
                     MPI_COMM_WORLD, &
                     ierr            )

    return
  end subroutine gather_grid

  !> communication: gather [2D]
  !---------------------------------------------------------------------------
  subroutine gather_vars( &
      invar,     & ! [in ]
      recvbuf    ) ! [out]
    implicit none

    real(SP), intent(in)  :: invar (:,:)
    real(SP), intent(out) :: recvbuf(:,:)

    integer :: sendcounts
    integer :: recvcounts
    integer :: iix, jjy
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = mnxp * mnyp * nmnge 
    recvcounts = mnxp * mnyp * nmnge 

    sendbuf(:,:) = UNDEF_SP
    do jjy = 1, mnyp*nmnge 
    do iix = 1, mnxp
       sendbuf(iix,jjy) = invar(iix,jjy)
    enddo
    enddo

    call MPI_GATHER( sendbuf(:,:), &
                     sendcounts,      &
                     MPI_REAL,        &
                     recvbuf(:,:), &
                     recvcounts,      &
                     MPI_REAL,        &
                     master,          &
                     MPI_COMM_WORLD,  &
                     ierr             )

    return
  end subroutine gather_vars


  !> create control file
  !---------------------------------------------------------------------------
  subroutine create_ctl( &
      varname,  & ! [in]
      idom,     & ! [in]
      it,       & ! [in]
      zz        ) ! [in]
    implicit none

    character(CMID), intent(in) :: varname
    integer,         intent(in) :: idom, it, zz

    character(2)    :: cdom
    character(3)    :: clev
    character(CLNG) :: fname
    character(CLNG) :: fname2
    !---------------------------------------------------------------------------

     write(cdom,'(i2.2)') idom
     select case( atype )
     case ( a_slice )
        if ( Z_MERGE_OUT ) then
           clev = "-3d"
        else
           write(clev,'(i3.3)') zz
        endif
     case ( a_max )
        clev = "max"
     case ( a_min )
        clev = "min"
     case ( a_sum )
        clev = "sum"
     case ( a_ave )
        clev = "ave"
     end select
     fname = trim(ODIR)//'/'//trim(varname)//'_d'//cdom//'z'//clev//'_'//trim(FTIME)
     fname2 = trim(varname)//'_d'//cdom//'z'//clev//'_'//trim(FTIME)
     if ( LOUT ) write( FID_LOG, '(1X,A,A)') "Create ctl file: ", trim(fname)
     open ( FID_CTL, file=trim(fname)//".ctl", form="formatted", access="sequential" )

     write( FID_CTL, '(a,1x,a)') "DSET", "^"//trim(fname2)//".grd"
     write( FID_CTL, '(a)') "TITLE SCALE3 data output"
     write( FID_CTL, '(a)') "OPTIONS BIG_ENDIAN"
     write( FID_CTL, '(a,1x,e15.7)') "UNDEF", -9.9999001E+30
     write( FID_CTL, '(a,3x,i7,1x,a)') "XDEF", nx, "LEVELS"
     write( FID_CTL, '(5(1x,e15.7))') cx(1:nx)*1.d-3
     write( FID_CTL, '(a,3x,i7,1x,a)') "YDEF", ny, "LEVELS"
     write( FID_CTL, '(5(1x,e15.7))') cy(1:ny)*1.d-3

     select case( vtype )
     case ( vt_urban, vt_land, vt_tpmsk )
        write( FID_CTL, '(a,3x,i7,1x,a,1x,a)') "ZDEF", 1, "linear", "1 1"
     case ( vt_height, vt_2d )
        write( FID_CTL, '(a,3x,i7,1x,a,1x,e15.7)') "ZDEF", 1, "LEVELS", cz(zz)*1.d-3
     case ( vt_3d )
        if ( Z_MERGE_OUT ) then
           write( FID_CTL, '(a,3x,i7,1x,a)') "ZDEF", ZCOUNT, "LEVELS"
           write( FID_CTL, '(5(1x,e15.7))') vgrid(1:ZCOUNT)*1.d-3
        else
           write( FID_CTL, '(a,3x,i7,1x,a,1x,e15.7)') "ZDEF", 1, "LEVELS", cz(zz)*1.d-3
        endif

     end select

     write( FID_CTL, '(a,3x,i5,1x,a,1x,a,3x,a)') "TDEF", 1, "LINEAR", trim(STIME), trim(DELT)
     write( FID_CTL, '(a,3x,i2)') "VARS", 1
     if ( Z_MERGE_OUT ) then
        write( FID_CTL, '(a,1x,i7,1x,i2,1x,a)') trim(varname), ZCOUNT, 99, "NONE"
     else
        write( FID_CTL, '(a,1x,i7,1x,i2,1x,a)') trim(varname), 0, 99, "NONE"
     endif
     write( FID_CTL, '(a)') "ENDVARS"

     close( FID_CTL )

    return
  end subroutine create_ctl


  !> write data file
  !---------------------------------------------------------------------------
  subroutine write_vars( &
      var_2d,   & ! [in]
      varname,  & ! [in]
      idom,     & ! [in]
      it,       & ! [in]
      zz        ) ! [in]
    implicit none

    real(SP),        intent(in) :: var_2d(:,:)
    character(CMID), intent(in) :: varname
    integer,         intent(in) :: idom, it, zz

    integer         :: irec
    integer(8)      :: irecl
    character(2)    :: cdom
    character(3)    :: clev
    character(CLNG) :: fname
    !---------------------------------------------------------------------------

    irec  = 1
    irecl = int(nx,kind=8) * int(ny,kind=8) * 4_8

    write(cdom,'(i2.2)') idom
    select case( atype )
    case ( a_slice )
       write(clev,'(i3.3)') zz
    case ( a_max )
       clev = "max"
    case ( a_min )
       clev = "min"
    case ( a_sum )
       clev = "sum"
    case ( a_ave )
       clev = "ave"
    end select
    fname = trim(ODIR)//'/'//trim(varname)//'_d'//cdom//'z'//clev//'_'//trim(FTIME)//".grd"
    if ( Z_MERGE_OUT ) fname_save = fname
    if ( LOUT ) write( FID_LOG, '(1X,A,A)') "+++ Output data file: ", trim(fname)
    if ( LOUT ) write( FID_LOG, *) "+++ Check data range: ", maxval(var_2d), minval(var_2d)

    open( FID_DAT, file=trim(fname), form="unformatted", access="direct", recl=irecl)
    write( FID_DAT, rec=irec ) var_2d(:,:)
    close( FID_DAT )

    return
  end subroutine write_vars


  !> write data file
  !---------------------------------------------------------------------------
  subroutine write_vars_zmerge( &
      var_2d,   & ! [in]
      irec      ) ! [in]
    implicit none

    real(SP),        intent(in) :: var_2d(:,:)
    integer,         intent(in) :: irec

    integer(8)      :: irecl
    !---------------------------------------------------------------------------

    irecl = int(nx,kind=8) * int(ny,kind=8) * 4_8
    if ( LOUT ) write( FID_LOG, '(1X,A,A)') "+++ Merged output to: ", trim(fname_save)
    if ( LOUT ) write( FID_LOG, *) "+++ Check data range: ", maxval(var_2d), minval(var_2d)

    open( FID_DAT, file=trim(fname_save), form="unformatted", access="direct", recl=irecl)
    write( FID_DAT, rec=irec ) var_2d(:,:)
    close( FID_DAT )

    return
  end subroutine write_vars_zmerge


  !> read configulation namelists
  !---------------------------------------------------------------------------
  subroutine read_conf_logout()
    implicit none
    !---------------------------------------------------------------------------

    open ( FID_CONF, file=trim(fconf), status='old', &
           form="formatted", delim='apostrophe', iostat=ierr )
    if ( ierr /= 0 ) then
       write (*, *) "ERROR: fail to open net2g.conf file"
       call err_abort( 1, __LINE__ )
    endif

    read  ( FID_CONF, nml=LOGOUT, iostat=ierr )
    if ( ierr /= 0 ) then
       write (*, *) "ERROR: fail to read LOGOUT"
       call err_abort( 1, __LINE__ )
    endif

    close ( FID_CONF )

    return
  end subroutine read_conf_logout


  !> read configulation namelists
  !---------------------------------------------------------------------------
  subroutine read_conf()
    implicit none

    integer :: n, m
    !---------------------------------------------------------------------------

    !--- read namelist file
    open ( FID_CONF, file=trim(fconf), status='old', &
           form="formatted", delim='apostrophe', iostat=ierr )
    if ( ierr /= 0 ) then
       if ( LOUT ) write (*, *) "ERROR: fail to open net2g.conf file"
       call err_abort( 1, __LINE__ )
    endif

    read  ( FID_CONF, nml=INFO, iostat=ierr )
    if ( LOUT ) write ( FID_LOG, nml=INFO )

    if ( ZCOUNT > max_zcount ) then
       if ( LOUT ) write (*, *) "ERROR: overflow maximum of zcount"
       call err_abort( 1, __LINE__ )
    endif

    if ( VCOUNT > max_vcount ) then
       if ( LOUT ) write (*, *) "ERROR: overflow maximum of vcount"
       call err_abort( 1, __LINE__ )
    elseif( VCOUNT < 1 ) then
       if ( LOUT ) write (*, *) "ERROR: specify at least one target variable"
       call err_abort( 1, __LINE__ )
    endif

    rewind( FID_CONF )
    read  ( FID_CONF, nml=VARI, iostat=ierr )
    if ( LOUT ) write ( FID_LOG, nml=VARI )
    if ( LOUT ) write( FID_LOG,* ) ""
    do n=1, VCOUNT
       VNAME(n) = trim(VNAME(n))
       if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,A)' ) "+++ Listing Vars: (", n, ") ", VNAME(n)
    enddo
    if ( LOUT ) write( FID_LOG,* ) ""

    if ( Z_LEV_LIST ) then
       do n=1, ZCOUNT
          if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5)' ) "+++ Listing Levs: (", n, ") ", TARGET_ZLEV(n)
       enddo
    else
       m = ZSTART
       do n=1, ZCOUNT
          TARGET_ZLEV(n) = m
          if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5)' ) "+++ Listing Levs: (", n, ") ", TARGET_ZLEV(n)
          m = m + 1
       enddo
    endif
    if ( LOUT ) write( FID_LOG,* ) ""

    rewind( FID_CONF )
    read  ( FID_CONF, nml=GRADS, iostat=ierr )
    if ( LOUT ) write ( FID_LOG, nml=GRADS )

    close ( FID_CONF )


    !--- read run.conf file
    open ( FID_RCNF, file=trim(CONFFILE), status='old', &
           form="formatted", delim='apostrophe', iostat=ierr )
    if ( ierr /= 0 ) then
       if ( LOUT ) write (*, *) "ERROR: fail to open running *.conf file"
       call err_abort( 1, __LINE__ )
    endif

    read  ( FID_RCNF, nml=PARAM_TIME, iostat=ierr )
    if ( LOUT ) write ( FID_LOG, nml=PARAM_TIME )

    read  ( FID_RCNF, nml=PARAM_PRC, iostat=ierr )
    if ( LOUT ) write ( FID_LOG, nml=PARAM_PRC )

    rewind( FID_RCNF )
    read  ( FID_RCNF, nml=PARAM_HIST, iostat=ierr )
    if ( LOUT ) write ( FID_LOG, nml=PARAM_HIST )

    rewind( FID_RCNF )
    read  ( FID_RCNF, nml=PARAM_HISTORY, iostat=ierr )
    if ( LOUT ) write ( FID_LOG, nml=PARAM_HISTORY )

    close(FID_RCNF)

    if ( EXTRA_TINTERVAL > 0 ) then
       HISTORY_DEFAULT_TINTERVAL = EXTRA_TINTERVAL
       HISTORY_DEFAULT_TUNIT     = EXTRA_TUNIT
       if ( LOUT ) write( FID_LOG, '(1X,A,F5.1,A)' ) "+++ USE EXTRA TIME: ", &
                        HISTORY_DEFAULT_TINTERVAL, trim(HISTORY_DEFAULT_TUNIT)
    endif

    !--- tentative
    if ( HIST_BND ) then
       if ( LOUT ) write (*, *) "HIST_BND is currently unsupported"
       call err_abort( 1, __LINE__ )
    endif

    return
  end subroutine read_conf


  !> setting of indices
  !---------------------------------------------------------------------------
  subroutine irank2ixjy( &
      irank,    & ! [in ]
      ix,       & ! [out]
      jy        ) ! [out]
    implicit none

    integer, intent(in)  :: irank
    integer, intent(out) :: ix, jy
    !---------------------------------------------------------------------------

    ix = mod( irank, xproc )  + 1
    jy = int( irank / xproc ) + 1

    return
  end subroutine irank2ixjy


  !> setting of flag boundary
  !---------------------------------------------------------------------------
  subroutine set_vtype( &
      ncfile,    & ! [in ]
      varname,   & ! [in ]
      vtype,     & ! [out]
      ndim       ) ! [out]
    implicit none

    character(CLNG), intent(in) :: ncfile
    character(CMID), intent(in) :: varname
    integer, intent(out)        :: vtype
    integer, intent(out)        :: ndim

    integer :: ncid, varid
    integer :: istat
    !---------------------------------------------------------------------------

    istat = nf90_open( trim(ncfile), nf90_nowrite, ncid )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)
 
    istat = nf90_inq_varid( ncid, trim(varname), varid )
    istat = nf90_inquire_variable(ncid, varid, ndims=ndim )
    if (istat .ne. nf90_noerr) call handle_err(istat, __LINE__)

    istat = nf90_close(ncid)

    select case( trim(varname) )
    case ( "TRL_URB", "TBL_URB", "TGL_URB" )
       vtype = vt_urban
    case ( "LAND_TEMP", "LAND_WATER" )
       vtype = vt_land
    case ( "height" )
       vtype = vt_height
    case ( "topo", "lsmask" )
       vtype = vt_tpmsk
       Z_MERGE_OUT = .false.
    case default
       if( ndim == 4 ) then
          vtype = vt_3d
       elseif( ndim == 3 ) then
          vtype = vt_2d
          Z_MERGE_OUT = .false.
       else
          call err_abort( 0, __LINE__ )
       end if
    end select

    return
  end subroutine set_vtype


  !> setting of flag analysis
  !---------------------------------------------------------------------------
  subroutine set_atype()
    implicit none
    !---------------------------------------------------------------------------

    select case( trim(ANALYSIS) )
    case ( "SLICE", "slice" )
       atype = a_slice
    case ( "MAX", "max", "MAXIMUM", "maximum" )
       atype = a_max
       ZCOUNT  = 1
       Z_MERGE_OUT = .false.
    case ( "MIN", "min", "MINIMUM", "minimum" )
       atype = a_min
       ZCOUNT  = 1
       Z_MERGE_OUT = .false.
    case ( "SUM", "sum", "SUMMATION", "summation" )
       atype = a_sum
       ZCOUNT  = 1
       Z_MERGE_OUT = .false.
    case ( "AVE", "ave", "AVERAGE", "average" )
       atype = a_ave
       ZCOUNT  = 1
       Z_MERGE_OUT = .false.
    case default
       if ( LOUT ) write (*, *) "ERROR: specified analysis type is not appropiate"
       if ( LOUT ) write (*, *) "***** ", trim(ANALYSIS)
       call err_abort( 1, __LINE__ )
    end select

    return
  end subroutine set_atype


  !> setting of flag boundary
  !---------------------------------------------------------------------------
  subroutine set_flag_bnd( &
      ix,       & ! [in ]
      jy,       & ! [in ]
      flag_bnd  ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    logical, intent(out) :: flag_bnd
    !---------------------------------------------------------------------------

    flag_bnd = .false.

    if ( HIST_BND ) then
       if ( ix == 1 .or. ix == xproc ) then
          flag_bnd = .true.
       endif
       if ( jy == 1 .or. jy == yproc ) then
          flag_bnd = .true.
       endif
    endif

    return
  end subroutine set_flag_bnd


  !> setting of indices
  !---------------------------------------------------------------------------
  subroutine set_index( &
      ix,      & ! [in ]
      jy,      & ! [in ]
      is,      & ! [out]
      ie,      & ! [out]
      js,      & ! [out]
      je       ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    integer, intent(out) :: is, ie
    integer, intent(out) :: js, je
    !---------------------------------------------------------------------------

    if ( HIST_BND ) then
       if ( ix == 1 ) then
          is = 1
          ie = nxp
       elseif ( ix == xproc ) then
          is = (ix-1)*(nxp-2)+2 + 1
          ie = (ix-1)*(nxp-2)+2 + nxp
       else
          is = (ix-1)*nxp+2 + 1
          ie = (ix-1)*nxp+2 + nxp
       endif

       if ( jy == 1 ) then
          js = 1
          je = nyp
       elseif ( jy == yproc ) then
          js = (jy-1)*(nyp-2)+2 + 1
          je = (jy-1)*(nyp-2)+2 + nyp
       else
          js = (jy-1)*nyp+2 + 1
          je = (jy-1)*nyp+2 + nyp
       endif
    else
       is = (ix-1)*nxp + 1
       ie = (ix-1)*nxp + nxp
       js = (jy-1)*nyp + 1
       je = (jy-1)*nyp + nyp
    endif

    return
  end subroutine set_index


  !> setting of indices for grid data
  !---------------------------------------------------------------------------
  subroutine set_index_grid( &
      ix,      & ! [in ]
      jy,      & ! [in ]
      is,      & ! [out]
      ie,      & ! [out]
      js,      & ! [out]
      je       ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    integer, intent(out) :: is, ie
    integer, intent(out) :: js, je
    !---------------------------------------------------------------------------

    is = (ix-1)*nxgp + 1
    ie = (ix-1)*nxgp + nxgp
    js = (jy-1)*nygp*xproc + 1
    je = (jy-1)*nygp*xproc + nygp

    return
  end subroutine set_index_grid


  !> setting of indices for grid data
  !---------------------------------------------------------------------------
  subroutine set_index_gathered( &
      ix,      & ! [in ]
      jy,      & ! [in ]
      is,      & ! [out]
      ie,      & ! [out]
      js,      & ! [out]
      je       ) ! [out]
    implicit none

    integer, intent(in)  :: ix, jy
    integer, intent(out) :: is, ie
    integer, intent(out) :: js, je
    !---------------------------------------------------------------------------

    if ( HIST_BND ) then
       if ( ix == 1 .or. ix == xproc ) then
          is = 1
          ie = mnxp
       else
          is = 1
          ie = mnxp - 2
       endif

       if ( jy == 1 ) then
          js = 1
          je = mnyp
       elseif ( jy == yproc ) then
          js = (ix-1)*mnyp + (jy-1)*mnyp*xproc + 1
          je = (ix-1)*mnyp + (jy-1)*mnyp*xproc + mnxp
       else
          js = (ix-1)*mnyp + (jy-1)*mnyp*xproc + 1
          je = (ix-1)*mnyp + (jy-1)*mnyp*xproc + (mnxp-2)
       endif
    else
       is = 1
       ie = (mnxp-2)

       js = (ix-1)*mnyp + (jy-1)*mnyp*xproc + 1
       je = (ix-1)*mnyp + (jy-1)*mnyp*xproc + (mnyp-2)
    endif

    return
  end subroutine set_index_gathered


  !> setting of indices for reading buffer
  !---------------------------------------------------------------------------
  subroutine set_index_readbuf( &
      nm,      & ! [in ]
      is,      & ! [out]
      ie,      & ! [out]
      js,      & ! [out]
      je,      & ! [out]
      im_bnd   ) ! [in ] optional
    implicit none

    integer, intent(in)  :: nm               ! num of loop for manage ranks
    integer, intent(out) :: is, ie           ! start index, end index
    integer, intent(out) :: js, je           ! start index, end index
    logical, intent(in), optional :: im_bnd  ! flag of boundary (edge tile)

    logical :: not_bnd
    !---------------------------------------------------------------------------

    not_bnd = .true.
    if ( present(im_bnd) ) then
       if ( im_bnd ) not_bnd = .false.
    endif

    if ( HIST_BND ) then
       if ( not_bnd ) then
          is = 1
          ie = nxp
          js = (nm-1)*mnyp + 1
          je = (nm-1)*mnyp + nyp
       else
          is = 1
          ie = mnxp
          js = (nm-1)*mnyp + 1
          je = (nm-1)*mnyp + mnyp
       endif
    else
       is = 1
       ie = nxp
       js = (nm-1)*nyp + 1
       je = (nm-1)*nyp + nyp
    endif

    return
  end subroutine set_index_readbuf


  !> setting of indices for reading buffer for grid
  !---------------------------------------------------------------------------
  subroutine set_index_readbuf_grid( &
      nm,      & ! [in ]
      is,      & ! [out]
      ie,      & ! [out]
      js,      & ! [out]
      je       ) ! [out]
    implicit none

    integer, intent(in)  :: nm               ! num of loop for manage ranks
    integer, intent(out) :: is, ie           ! start index, end index
    integer, intent(out) :: js, je           ! start index, end index
    !---------------------------------------------------------------------------

    is = (nm-1)*nxgp + 1
    ie = (nm-1)*nxgp + nxgp
    js = (nm-1)*nygp + 1
    je = (nm-1)*nygp + nygp

    return
  end subroutine set_index_readbuf_grid


  !> setting of indices for reading buffer (netcdf) with counts
  !---------------------------------------------------------------------------
  subroutine set_index_netcdf( &
      atype,  & ! [in ]
      vtype,  & ! [in ]
      ix,     & ! [in ]
      jy,     & ! [in ]
      it,     & ! [in ]
      zz,     & ! [in ]
      is,     & ! [out]
      ie,     & ! [out]
      js,     & ! [out]
      je,     & ! [out]
      nzn     ) ! [out]
    implicit none

    integer, intent(in)  :: atype, vtype
    integer, intent(in)  :: ix, jy
    integer, intent(in)  :: it, zz
    integer, intent(out) :: is, ie  ! start index, end index
    integer, intent(out) :: js, je  ! start index, end index
    integer, intent(out) :: nzn     ! number of z-dimension in netcdf
    !---------------------------------------------------------------------------

    is = 1
    js = 1

    if ( HIST_BND ) then
       if ( ix == 1 .or. ix == xproc ) then
          ie = mnxp
       else
          ie = mnxp - 2
       endif
       if ( jy == 1 .or. jy == yproc ) then
          je = mnyp
       else
          je = mnyp - 2
       endif
    else
       ie = nxp
       je = nyp
    endif

    count_tpmsk(1:2) = (/ ie, je    /)
    start_2d   (1:3) = (/ 1,  1, it /)
    start_2dt  (1:2) = (/ 1,  1     /)

    select case( atype )
    case ( a_slice )
       count_3d    (1:4) = (/ ie, je, 1,  1  /)
       count_2d    (1:3) = (/ ie, je,     1  /)
       count_urban (1:4) = (/ ie, je, 1,  1  /)
       count_land  (1:4) = (/ ie, je, 1,  1  /)
       count_height(1:3) = (/ ie, je, 1      /)
       start_3d    (1:4) = (/ 1,  1,  zz, it /)
       nzn = 1
    case ( a_max, a_min, a_sum, a_ave )
       count_3d    (1:4) = (/ ie, je, nz, 1  /)
       count_2d    (1:3) = (/ ie, je,     1  /)
       count_urban (1:4) = (/ ie, je, uz, 1  /)
       count_land  (1:4) = (/ ie, je, lz, 1  /)
       count_height(1:3) = (/ ie, je, nz     /)
       start_3d    (1:4) = (/ 1,  1,  1,  it /)

       select case( vtype )
       case ( vt_urban )
          nzn = uz
       case ( vt_land )
          nzn = lz
       case ( vt_3d )
          nzn = nz
       case ( vt_2d, vt_height, vt_tpmsk )
          if ( LOUT ) write (*, *) "ERROR: specified anal-type is not appropiate for the var"
          if ( LOUT ) write (*, *) "***** ", trim(ANALYSIS), trim(varname)
          call err_abort( 1, __LINE__ )
       case default
          call err_abort( 0, __LINE__ )
       end select
  
    case default
       call err_abort( 0, __LINE__ )
    end select

    return
  end subroutine set_index_netcdf


  !> setting of ranks should be managed in the process
  !---------------------------------------------------------------------------
  subroutine set_rank_manage( &
      irank,     & ! [in ]
      nmnge,     & ! [in ]
      rk_mnge    ) ! [out]
    implicit none

    integer, intent(in)  :: irank
    integer, intent(in)  :: nmnge
    integer, intent(out) :: rk_mnge(:)

    integer :: n
    !---------------------------------------------------------------------------

    if ( LOUT ) write( FID_LOG, '(1X,A,I5)') &
                "+++ number of ranks to manage: ", nmnge
    do n = 1, nmnge
       rk_mnge(n) = irank * nmnge + (n-1)
       if ( LOUT ) write( FID_LOG, '(1X,A,I5,A,I5)') &
                   "+++ myrank: ", irank, " -->  manage: ", rk_mnge(n)
    enddo

    return
  end subroutine set_rank_manage


  !> setting of indices
  !---------------------------------------------------------------------------
  subroutine set_array_size()
    implicit none
    !---------------------------------------------------------------------------

    if ( HIST_BND ) then
       mnxp = nxgp - 2
       mnyp = nygp - 2
       nx = (mnxp-2) * xproc + 4
       ny = (mnyp-2) * yproc + 4
    else
       mnxp = nxp
       mnyp = nyp
       nx = mnxp * xproc
       ny = mnyp * yproc
    endif

    mnx = mnxp * xproc
    mny = mnyp * yproc
    nxg_tproc = nxgp * tproc
    nyg_tproc = nygp * tproc

    if ( LOUT ) write( FID_LOG, '(1X,"+++ nx:",I7,2X,"ny:",I7)') nx, ny

    return
  end subroutine set_array_size


  !> setting of calender indices
  !---------------------------------------------------------------------------
  subroutine set_calender()
    implicit none

    character(2) :: cmn, chh, cdd, cmm2
    character(4) :: cyy
    !---------------------------------------------------------------------------

    write(cmn, '(I2.2)') mn
    write(chh, '(I2.2)') hh
    write(cdd, '(I2.2)') dd
    write(cmm2,'(I2.2)') mm
    write(cyy, '(I4.4)') yy
    STIME = chh//':'//cmn//'Z'//cdd//cmm(mm)//cyy
    FTIME = cyy//cmm2//cdd//chh//cmn

    return
  end subroutine set_calender


  !> setting of indices
  !---------------------------------------------------------------------------
  subroutine check_targ_zlev( &
      zz   ) ! [in]
    implicit none

    integer, intent(in)  :: zz
    !---------------------------------------------------------------------------

    select case( vtype )
    case ( vt_urban )
       if ( zz < 1 .or. zz > uz ) then
          if ( LOUT ) write (*, *) "ERROR: requested level is in K-HALO [urban]"
          call err_abort( 1, __LINE__ )
       endif
    case ( vt_land )
       if ( zz < 1 .or. zz > lz ) then
          if ( LOUT ) write (*, *) "ERROR: requested level is in K-HALO [land]"
          call err_abort( 1, __LINE__ )
       endif
    case ( vt_height )
       if ( zz < ks .or. zz > ke ) then
          if ( LOUT ) write (*, *) "ERROR: requested level is in K-HALO [height]"
          call err_abort( 1, __LINE__ )
       endif
    case ( vt_3d )
       if ( zz < ks .or. zz > ke ) then
          if ( LOUT ) write (*, *) "ERROR: requested level is in K-HALO [3D data]"
          call err_abort( 1, __LINE__ )
       endif
    case ( vt_2d, vt_tpmsk )
       if ( zz /= 1 ) then
          if ( LOUT ) write (*, *) "ERROR: requested level is not exist [2D data]"
          call err_abort( 1, __LINE__ )
       endif
    case default
        call err_abort( 0, __LINE__ )
    end select

    return
  end subroutine check_targ_zlev


  !> make vgrid for merged-z output
  !---------------------------------------------------------------------------
  subroutine make_vgrid()
    implicit none

    integer :: iz, ik, k
    real(DP) :: diff, mini
    !---------------------------------------------------------------------------

    allocate( vgrid( ZCOUNT ) )

    select case( Z_LEV_TYPE )
    case ( "GRID", "grid" )
       do iz = 1, ZCOUNT
          vgrid(iz) = cz( TARGET_ZLEV(iz) )
          if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,F8.2)' ) "+++ Target Height: (", iz, ") ", vgrid(iz)
       enddo

    case ( "HEIGHT", "height", "HGT", "hgt" )
       do iz = 1, ZCOUNT
          mini = 9.9999D+10
          do k=nzh+1, nz+nzh+1
             diff = abs(real(TARGET_ZLEV(iz)) - cz(k))
             if ( diff < mini ) then
                mini = diff
                ik = k
             endif
          enddo
          if ( LOUT ) write( FID_LOG, '(1X,A,I5,A,F8.3)' ) &
          "+++ Search Nearest - request: ", TARGET_ZLEV(iz), "  diff: ", mini

          TARGET_ZLEV(iz) = ik
          vgrid(iz) = cz( TARGET_ZLEV(iz) )
          if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,F8.2)' ) "+++ Target Height: (", iz, ") ", vgrid(iz)
       enddo

    case default
       if ( LOUT ) write (*, *) "ERROR: requested Z_LEV_TYPE is not supported"
       call err_abort( 1, __LINE__ )
    end select

    return
  end subroutine make_vgrid


  !> calender initialization
  !---------------------------------------------------------------------------
  subroutine cal_init()
    implicit none

    real(DP) :: inc
    integer  :: i, tint
    character(2) :: tunit
    character(3) :: cint
    !---------------------------------------------------------------------------

    yy = TIME_STARTDATE(1)
    mm = TIME_STARTDATE(2)
    dd = TIME_STARTDATE(3)
    hh = TIME_STARTDATE(4)
    mn = TIME_STARTDATE(5)
    sc = TIME_STARTDATE(6)

    if ( START_TSTEP > 1 ) then
       do i=1, START_TSTEP-1
          call cal_increment
       enddo
    endif

    inc = HISTORY_DEFAULT_TINTERVAL * dble(INC_TSTEP)
    select case ( HISTORY_DEFAULT_TUNIT )
    case ( "SEC", "sec" )
       if ( inc < 60.0D0 ) then
          if ( LOUT ) write( FID_LOG, '(1X,A)') &
                      "*** WARNING: HISTORY_DEFAULT_TINTERVAL is not compatible!"
          if ( LOUT ) write( FID_LOG, '(1X,A,I7,A)') &
                      "*** ", int(HISTORY_DEFAULT_TINTERVAL), " is too short for Grads"
          tint  = 1     !tentative
          tunit = "mn"  !tentative
       else
          inc = inc / 60.0D0
          if ( inc < 60.0D0 ) then
             tint  = int(inc)
             tunit = "mn"
          else
             inc = inc / 60.0D0
             tint  = int(inc)
             tunit = "hr"
          endif
       endif
    case ( "MIN", "min" )
       if ( inc < 60.0D0 ) then
          tint  = int(inc)
          tunit = "mn"
       else
          inc = inc / 60.0D0
          tint  = int(inc)
          tunit = "hr"
       endif
    case ( "HOUR", "hour" )
       tint  = int(inc)
       tunit = "hr"
    case default
        call err_abort( 0, __LINE__ )
    end select

    if ( tint <= 0 ) then
       tint = 1  ! avoid zero
    endif
    write(cint,'(I3)') tint
    DELT = trim(cint)//tunit

    return
  end subroutine cal_init


  !> calender increment calculation
  !---------------------------------------------------------------------------
  subroutine cal_increment()
    implicit none

    real(DP) :: inc
    !---------------------------------------------------------------------------

    inc = HISTORY_DEFAULT_TINTERVAL * dble(INC_TSTEP)

    select case ( HISTORY_DEFAULT_TUNIT )
    case ( "SEC", "sec" )
       if ( inc < 60.0D0 ) then
          call cal_inc_sec( int(inc) )
       else
          inc = inc / 60.0D0
          if ( inc < 60.0D0 ) then
             call cal_inc_min( int(inc) )
          else
             inc = inc / 60.0D0
             call cal_inc_hour( int(inc) )
          endif
       endif
    case ( "MIN", "min" )
       if ( inc < 60.0D0 ) then
          call cal_inc_min( int(inc) )
       else
          inc = inc / 60.0D0
          call cal_inc_hour( int(inc) )
       endif
    case ( "HOUR", "hour" )
       call cal_inc_hour( int(inc) )
    case default
        call err_abort( 0, __LINE__ )
    end select

    return
  end subroutine cal_increment


  !> calender increment calculation: for [sec] increment
  !---------------------------------------------------------------------------
  subroutine cal_inc_sec( &
      inc  ) ! [in]
    implicit none

    integer, intent(in)  :: inc

    integer  :: eday
    !---------------------------------------------------------------------------

    call cal_date( yy, mm, eday )

    sc = sc + inc
    if ( sc >= 60 ) then
       mn = mn + 1
       sc = sc - 60
       if ( mn >= 60 ) then
          hh = hh + 1
          mn = mn - 60
          if ( hh >= 24 ) then
             dd = dd + 1
             hh = hh - 24
             if ( dd > eday ) then
                mm = mm + 1
                dd = dd - eday
                if ( mm > 12 ) then
                   yy = yy + 1
                   mm = mm - 12
                endif
             endif
          endif
       endif
    endif

    return
  end subroutine cal_inc_sec


  !> calender increment calculation: for [min] increment
  !---------------------------------------------------------------------------
  subroutine cal_inc_min( &
      inc  ) ! [in]
    implicit none

    integer, intent(in)  :: inc

    integer  :: eday
    !---------------------------------------------------------------------------

    call cal_date( yy, mm, eday )

    mn = mn + inc
    if ( mn >= 60 ) then
       hh = hh + 1
       mn = mn - 60
       if ( hh >= 24 ) then
          dd = dd + 1
          hh = hh - 24
          if ( dd > eday ) then
             mm = mm + 1
             dd = dd - eday
             if ( mm > 12 ) then
                yy = yy + 1
                mm = mm - 12
             endif
          endif
       endif
    endif

    return
  end subroutine cal_inc_min


  !> calender increment calculation: for [hour] increment
  !---------------------------------------------------------------------------
  subroutine cal_inc_hour( &
      inc  ) ! [in]
    implicit none

    integer, intent(in)  :: inc

    integer  :: eday
    !---------------------------------------------------------------------------

    call cal_date( yy, mm, eday )

    hh = hh + inc
    if ( hh >= 24 ) then
       dd = dd + 1
       hh = hh - 24
       if ( dd > eday ) then
          mm = mm + 1
          dd = dd - eday
          if ( mm > 12 ) then
             yy = yy + 1
             mm = mm - 12
          endif
       endif
    endif

    return
  end subroutine cal_inc_hour


  !> calender date calculation
  !---------------------------------------------------------------------------
  subroutine cal_date( &
      yy,  & ! [in ]
      mm,  & ! [in ]
      dd   ) ! [out]
    implicit none

    integer, intent(in)  :: yy
    integer, intent(in)  :: mm
    integer, intent(out) :: dd

    integer :: rem4, rem100, rem400
    !---------------------------------------------------------------------------

    select case ( mm )
    case ( 4, 6, 9, 11 )
       dd = 30
    case ( 1, 3, 5, 7, 8, 10, 12 )
       dd = 31
    case ( 2 )
       rem4   = int( mod(real(yy), 4.0  ) )
       rem100 = int( mod(real(yy), 100.0) )
       rem400 = int( mod(real(yy), 400.0) )
       dd = 28
       if ( rem4 == 0 ) then            ! T -> leap year
          if ( rem100 == 0 ) then       ! F -> leap year
             if ( rem400 == 0 ) dd = 29 ! T -> leap year
          else
             dd = 29
          endif
       endif
    case default
        call err_abort( 0, __LINE__ )
    end select

    return
  end subroutine cal_date


  !> allocation of data arrays
  !---------------------------------------------------------------------------
  subroutine allocation( &
      irank   ) ! [in]
    implicit none

    integer, intent(in) :: irank
    !---------------------------------------------------------------------------

    allocate( p_2dt       (mnxp,       mnyp              ) )
    allocate( p_2d        (mnxp,       mnyp,       1     ) )
    allocate( p_3d        (mnxp,       mnyp,       nz, 1 ) )
    allocate( p_3d_urban  (mnxp,       mnyp,       uz, 1 ) )
    allocate( p_3d_land   (mnxp,       mnyp,       lz, 1 ) )
    allocate( p_var       (mnxp,       mnyp*nmnge        ) )
    allocate( p_cdx       (nxgp*nmnge                    ) )
    allocate( p_cdy       (nygp*nmnge                    ) )
    allocate( p_cx        (nxgp*nmnge                    ) )
    allocate( p_cy        (nygp*nmnge                    ) )
    allocate( cdz         (nz+2*nzh                      ) )
    allocate( cz          (nz+2*nzh                      ) )

    allocate( sendbuf     (mnxp,       mnyp*nmnge        ) )
    allocate( sendbuf_gx  (nxgp*nmnge                    ) )
    allocate( sendbuf_gy  (nygp*nmnge                    ) )

    if ( irank == master ) then
       allocate( var_2d    (nx,               ny               ) )
       allocate( cx        (nx                                 ) )
       allocate( cy        (ny                                 ) )
       allocate( cdx       (nx                                 ) )
       allocate( cdy       (ny                                 ) )
       allocate( recvbuf   (mnxp,             mnyp*nmnge*tproc ) )
       allocate( cx_gather (nxgp*nmnge*tproc                   ) )
       allocate( cy_gather (nygp*nmnge*tproc                   ) )
       allocate( cdx_gather(nxgp*nmnge*tproc                   ) )
       allocate( cdy_gather(nygp*nmnge*tproc                   ) )
    else
       allocate( recvbuf   (1,                1                ) )
       allocate( cx_gather (1                                  ) )
       allocate( cy_gather (1                                  ) )
       allocate( cdx_gather(1                                  ) )
       allocate( cdy_gather(1                                  ) )
    endif

    return
  end subroutine allocation


  !> error handler for netcdf90 system
  !---------------------------------------------------------------------------
  subroutine handle_err( &
      istat,  & ! [in]
      nline   ) ! [in]
    implicit none

    integer, intent(in) :: istat
    integer, intent(in) :: nline
    !---------------------------------------------------------------------------

    if ( LOUT ) write (*, *) nf90_strerror(istat)
    call err_abort( -1, nline )

    stop
  end subroutine handle_err


  !> process abort: error handing
  !---------------------------------------------------------------------------
  subroutine err_abort( &
      ecode,  & ! [in]
      nline   ) ! [in]
    implicit none

    integer, intent(in) :: ecode
    integer, intent(in) :: nline
    !---------------------------------------------------------------------------

    if ( ecode == err_internal ) then
       write (*, *) "##### ERROR: internal error"
    elseif ( ecode == err_netcdf ) then
       write (*, *) "##### ERROR: netcdf error"
    endif

    write (*, *) "***** Abort: by rank =", irank
    write (*, *) "*****        at Line =", nline

    if ( LOUT ) close ( FID_LOG )
    call MPI_ABORT( MPI_COMM_WORLD, ecode, ierr )

    stop
  end subroutine err_abort

end program netcdf2grads_h
