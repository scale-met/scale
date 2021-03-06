!-------------------------------------------------------------------------------------------
!> module NET2G main routine
!!
!! @par Description
!!          main routine of post-process for scale high-performance
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_netcdf2grads_h
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_net2g_vars
  use mod_net2g_error
  use mod_net2g_calender
  use mod_net2g_comm
  use mod_net2g_io
  use mod_net2g_netcdf
  use mod_net2g_setup
  use mod_net2g_anal

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: popsca

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: combine_grid
  private :: combine_vars_2d
  private :: read_conf_logout
  private :: read_conf
  private :: check_targ_zlev
  private :: make_vgrid
  private :: allocation

  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(SP),allocatable, private :: zlev(:), cz(:), cdz(:)
  real(SP),allocatable, private :: cx(:), cdx(:)
  real(SP),allocatable, private :: cy(:), cdy(:)
  real(SP),allocatable, private :: var_2d(:,:)
  real(SP),allocatable, private :: p_var(:,:)
  real(SP),allocatable, private :: p_ref(:,:,:)

  real(SP),allocatable, private :: recvbuf(:,:)
  real(SP),allocatable, private :: cx_gather(:)
  real(SP),allocatable, private :: cy_gather(:)
  real(SP),allocatable, private :: cdx_gather(:)
  real(SP),allocatable, private :: cdy_gather(:)

  real(SP),allocatable, private :: p_cx(:), p_cdx(:)     ! partial grid data
  real(SP),allocatable, private :: p_cy(:), p_cdy(:)     ! partial grid data

  real(SP),allocatable, private :: vgrid(:)

  integer, private :: vtype = -1
  integer, private :: atype = -1
  integer, private :: ctype = -1
  integer, private :: nt
  integer, private :: nx   ! num of x-dimension in the combined file
  integer, private :: ny   ! num of y-dimension in the combined file
  integer, private :: nz   ! num of z-dimension in the combined file
  integer, private :: uz   ! num of z-dimension in the combined file for urban
  integer, private :: lz   ! num of z-dimension in the combined file for land
  integer, private :: nzg  ! num of z-dimension in the combined file for grid
  integer, private :: nxh  ! num of halo grids for x-dimension
  integer, private :: nyh  ! num of halo grids for y-dimension
  integer, private :: nzh  ! num of halo grids for z-dimension
  integer, private :: nz_all(3)
  integer, private :: ks, ke
  integer, private :: nst, nen

  integer, private  :: xproc, yproc, tproc
  integer, private  :: it, ix, jy, iz, iv
  integer, private  :: zz
  integer, private  :: ndim
  integer, private  :: idom
  integer, private  :: ierr

  ! these are defined in "mod_net2g_comm2"
  !integer :: irank = 0
  !integer :: isize = 1

  integer, private :: nxp   ! num of x-dimension in partial file
  integer, private :: nyp   ! num of y-dimension in partial file
  integer, private :: nxgp  ! num of x-dimension in partial file for grid
  integer, private :: nygp  ! num of y-dimension in partial file for grid
  integer, private :: mnxp  ! maximum num of x-grid in partial file among all files
  integer, private :: mnyp  ! maximum num of y-grid in partial file among all files
  integer, private :: nxg_tproc
  integer, private :: nyg_tproc
  integer, private :: mnx   ! maximum num of x-dimension in partial file among all files
  integer, private :: mny   ! maximum num of y-dimension in partial file among all files

  integer, private :: nm
  integer, private :: nmnge  ! number of managing files
  integer, private :: work
  integer, private :: yy, mm, dd, hh, mn, sc
  integer,allocatable, private  :: rk_mnge(:)

  integer :: irec_time
  integer :: irec_timelev

  character(6), private     :: num
  character(CLNG), private  :: ncfile
  character(CMID), private  :: varname
  character(CMID), private  :: fconf
  character(CLNG), private  :: fname_save
  !-----------------------------------------------------------------------------------------
  namelist /LOGOUT/            &
    LOG_BASENAME,              &
    LOG_ALL_OUTPUT,            &
    LOG_LEVEL
  namelist /INFO/              &
    START_TSTEP,               &
    END_TSTEP,                 &
    INC_TSTEP,                 &
    DOMAIN_NUM,                &
    ZCOUNT,                    &
    ZSTART,                    &
    CONFFILE,                  &
    IDIR,                      &
    ODIR,                      &
    Z_LEV_TYPE,                &
    Z_MERGE_OUT,               &
    T_MERGE_OUT
  namelist /EXTRA/             &
    EXTRA_TINTERVAL,           &
    EXTRA_TUNIT
  namelist /ANAL/              &
    ANALYSIS
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
    HISTORY_DEFAULT_ZDIM   
  namelist  /PARAM_HIST/       &
    HIST_BND

  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------

  !> main routine of POPSCA
  !-----------------------------------------------------------------------------------------
  subroutine popsca( &
      LOCAL_COMM_WORLD,   &
      CONFFILE_IN         ) ! [in]
    implicit none

    integer,        intent(in) :: LOCAL_COMM_WORLD
    character(256), intent(in) :: CONFFILE_IN
    !---------------------------------------------------------------------------


  !-----------------------------------------------------------------------------------------

  !### initialization
  !-----------------------------------------------------------------------------------------
  call comm_initialize( LOCAL_COMM_WORLD )

  fconf = trim( CONFFILE_IN )

  call read_conf_logout
  call io_log_init( irank, LOUT )

  if ( LOUT ) write( FID_LOG, '(1X,A)') "-----------------------------------"
  if ( LOUT ) write( FID_LOG, '(1X,A)') "     START : netcdf to grads"
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
     write (*, *) "ERROR: specified num of mpi processes is not adequate."
     write (*, *) "*** specify the num; PRC_X*PRC_Y shuold be divisable by the num."
     call err_abort( 1, __LINE__, loc_main )
  elseif ( isize > tproc ) then
     write (*, *) "ERROR: num of mpi processes is larger than that of the scale-rm run."
     call err_abort( 1, __LINE__, loc_main )
  endif

  allocate ( rk_mnge (nmnge) )
  call set_rank_manage( irank, nmnge, rk_mnge )

  !### read and combine
  !-----------------------------------------------------------------------------------------
  write ( num,'(I6.6)' ) rk_mnge(1)
  ncfile = trim(IDIR)//"/"//trim(HISTORY_DEFAULT_BASENAME)//".pe"//num//".nc"
  call netcdf_retrieve_dims( ncfile,                          &
                             nxp, nxgp, nxh, nyp, nygp, nyh,  &
                             nz,  nzg,  nzh, uz, lz, ks, ke   )
  nz_all(1) = nz; nz_all(2) = uz; nz_all(3) = lz

  call set_array_size( nxp, nyp, nxgp, nygp, nx, ny, mnx, mny,  &
                       mnxp, mnyp, nxg_tproc, nyg_tproc         )
  call allocation( irank ) ! "irank" is right, not "rk_mnge"
  call comm_setup( mnxp, mnyp, nxgp, nygp, nmnge )
  call netcdf_setup( mnxp, mnyp, nz_all )

  call set_atype( atype )
  call set_ctype( ctype )
  if ( atype == a_conv ) then ! allocation
     call anal_setup( mnxp, mnyp, nz_all, nmnge )
     if ( ctype == c_height ) then !--- get reference height
        do nm = 1, nmnge
           it = START_TSTEP
           call netcdf_read_ref( rk_mnge(nm), nxp, nyp, mnxp, mnyp, &
                                 it, nz_all, ctype, p_ref )
           call anal_input_ref( p_ref, nm )
        enddo
     endif
  endif
  do nm = 1, nmnge
     call netcdf_read_grid( rk_mnge(nm), nm, nxgp, nygp,            &
                            zlev, cz, cdz, p_cx, p_cdx, p_cy, p_cdy )
  enddo

  call comm_gather_grid( nxgp, nygp, nmnge, p_cx, p_cdx, p_cy, p_cdy, &
                         cx_gather, cdx_gather, cy_gather, cdy_gather )

  if ( irank == master ) call combine_grid
  if ( Z_MERGE_OUT ) call make_vgrid

  !--- read data and combine
  nst = START_TSTEP
  nen = END_TSTEP
  nt  = (nen - nst + 1) / INC_TSTEP
  !if ( nt > max_tcount ) then
  !   write (*, *) "ERROR: overflow maximum of tcount"
  !   call err_abort( 1, __LINE__, loc_main )
  !endif

  call cal_init( yy, mm, dd, hh, mn, sc, DELT )

  irec_time = 1
  do it = nst, nen, INC_TSTEP !--- time loop

     call set_calender( yy, mm, dd, hh, mn, sc, STIME, FTIME )
     if ( LOUT ) write( FID_LOG, '(1X,A)' ) ""
     if ( LOUT ) write( FID_LOG, '(1X,A,I3)' ) "+++ TIME STEP:   ", it
     if ( LOUT ) write( FID_LOG, * ) "+++ Date: ", STIME

     if ( ctype == c_pres ) then !--- update reference pressure
        do nm = 1, nmnge
           call netcdf_read_ref( rk_mnge(nm), nxp, nyp, mnxp, mnyp, &
                                 it, nz_all, ctype, p_ref )
           call anal_input_ref( p_ref, nm )
        enddo
     endif

     do iv = 1, vcount !--- var loop
        varname = trim( VNAME(iv) )
        call netcdf_var_dim( ncfile, varname, ndim )
        call set_vtype( ndim, varname, vtype, atype )
        if ( LOUT ) write( FID_LOG, '(1X,A)' ) "|"
        if ( LOUT ) write( FID_LOG, '(1X,A,A,A,I1,A)' ) &
        "+++ VARIABLE: ", trim(varname), " (vtype = ", vtype, ")"

        select case( vtype )
        case ( vt_urban, vt_land, vt_height, vt_3d )
        !-----------------------------------------------------------------------------------

           irec_timelev = irec_time
           do iz = 1, ZCOUNT        !--- level loop
              if ( atype == a_slice ) then
                 zz = TARGET_ZLEV(iz)
                 if ( LOUT ) write( FID_LOG, '(1X,A,I6)' ) "+++ Z LEVEL: ", zz
                 call check_targ_zlev( zz )
              elseif ( atype == a_conv ) then
                 zz = TARGET_ZLEV(iz)
                 if ( ctype == c_pres ) then
                    if ( LOUT ) write( FID_LOG, '(1X,A,I6,A)' ) "+++ Z LEVEL: ", zz, " [hPa]"
                 elseif ( ctype == c_height ) then
                    if ( LOUT ) write( FID_LOG, '(1X,A,I6,A)' ) "+++ Z LEVEL: ", zz, " [m]"
                 endif
              endif

              do nm = 1, nmnge
                 call netcdf_read_var( rk_mnge(nm), nm, nxp, nyp, mnxp, mnyp,  &
                                       it, zz, nz_all, varname, atype, ctype, vtype, p_var )
              enddo

              call comm_gather_vars( mnxp, mnyp, nmnge, p_var, recvbuf )

              if ( irank == master ) then
                 call combine_vars_2d( recvbuf, var_2d )

                 if ( iz /= 1 .and. Z_MERGE_OUT ) then
                    call io_write_vars( var_2d, varname, atype, vtype, &
                                        idom, nx, ny, zz, irec_timelev )
                 else
                    if ( T_MERGE_OUT ) then
                       if (irec_time==1) call io_create_ctl( varname, atype, ctype, vtype, idom, &
                                                             nx, ny, zz, nt, cx, cy, vgrid, zlev )
                    else
                       call io_create_ctl( varname, atype, ctype, vtype, idom, &
                                           nx, ny, zz, 1, cx, cy, vgrid, zlev )
                    endif
                    call io_write_vars( var_2d, varname, atype, vtype, &
                                        idom, nx, ny, zz, irec_timelev )
                 endif
              endif

              irec_timelev = irec_timelev + 1
           enddo !--- level loop

        case ( vt_2d, vt_tpmsk )
        !-----------------------------------------------------------------------------------

           zz = 1
           do nm = 1, nmnge
              call netcdf_read_var( rk_mnge(nm), nm, nxp, nyp, mnxp, mnyp,   &
                                    it, zz, nz_all, varname, atype, ctype, vtype, p_var )
           enddo

           call comm_gather_vars( mnxp, mnyp, nmnge, p_var, recvbuf )

           !zz = 0 ! used as file name
           if ( irank == master ) then
              call combine_vars_2d( recvbuf, var_2d )

              if ( T_MERGE_OUT ) then
                 if (irec_time == 1) call io_create_ctl( varname, atype, ctype, vtype, idom, &
                                                         nx, ny, zz, nt, cx, cy, vgrid, zlev )
              else
                 call io_create_ctl( varname, atype, ctype, vtype, idom, &
                                     nx, ny, zz, 1, cx, cy, vgrid, zlev )
              endif
              call io_write_vars( var_2d, varname, atype, vtype, &
                                  idom, nx, ny, zz, irec_time )
           endif

        case default
        !-----------------------------------------------------------------------------------
           call err_abort( 0, __LINE__, loc_main )
        end select

     enddo !--- var loop

     call cal_increment( yy, mm, dd, hh, mn, sc )
    if ( ZCOUNT == 0 ) then
       write (*, *) "ERROR: at record increment"
       call err_abort( 0, __LINE__, loc_main )
    endif
     irec_time = irec_time + ZCOUNT
  enddo !--- time loop

  ! finalization
  if ( LOUT ) write( FID_LOG, '(1X,A)') ""
  if ( LOUT ) write( FID_LOG, '(1X,A)') "     FINISH : netcdf to grads"
  if ( LOUT ) write( FID_LOG, '(1X,A)') "-----------------------------------"
  call io_close_all

    return
  end subroutine popsca


  !> combine region divided data: grid data
  !-----------------------------------------------------------------------------------------
  subroutine combine_grid()
    implicit none

    integer :: ix, jy, iix, jjy
    integer :: is, ie, js, je
    integer :: gix, gjy
    integer :: gis, gie, gjs, gje
    !---------------------------------------------------------------------------

    jy = 1
    do ix = 1, xproc
       call set_index( ix, jy, nxp, nyp, is, ie, js, je )
       call set_index_grid( ix, jy, nxgp, nygp, gis, gie, gjs, gje )

       gix = gis
       do iix = is, ie
          cdx(iix) = cdx_gather(gix+nxh)
          cx(iix)  = cx_gather( gix+nxh)
          gix = gix + 1
       enddo
    enddo

    ix = 1
    do jy = 1, yproc
       call set_index( ix, jy, nxp, nyp, is, ie, js, je )
       call set_index_grid( ix, jy, nxgp, nygp, gis, gie, gjs, gje )

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
  !-----------------------------------------------------------------------------------------
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
       call set_index( ix, jy, nxp, nyp, is, ie, js, je )
       call set_index_gathered( ix, jy, mnxp, mnyp, gis, gie, gjs, gje )

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


  !> read configulation namelists
  !-----------------------------------------------------------------------------------------
  subroutine read_conf_logout()
    implicit none
    !---------------------------------------------------------------------------

    open ( FID_CONF, file=trim(fconf), status='old', &
           form="formatted", delim='apostrophe', iostat=ierr )
    if ( ierr /= 0 ) then
       write (*, *) "ERROR: fail to open configuration file"
       call err_abort( 1, __LINE__, loc_main )
    endif

    read  ( FID_CONF, nml=LOGOUT, iostat=ierr )

    close ( FID_CONF )
    return
  end subroutine read_conf_logout


  !> read configulation namelists
  !-----------------------------------------------------------------------------------------
  subroutine read_conf()
    implicit none

    integer :: n, m
    !---------------------------------------------------------------------------

    !--- read namelist file
    open ( FID_CONF, file=trim(fconf), status='old', &
           form="formatted", delim='apostrophe', iostat=ierr )
    if ( ierr /= 0 ) then
       write (*, *) "ERROR: fail to open net2g.conf file", trim(fconf)
       call err_abort( 1, __LINE__, loc_main )
    endif

    read  ( FID_CONF, nml=INFO, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=INFO )

    if ( ZCOUNT > max_zcount ) then
       write (*, *) "ERROR: overflow maximum of zcount"
       call err_abort( 1, __LINE__, loc_main )
    endif

    rewind( FID_CONF )
    read  ( FID_CONF, nml=EXTRA, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=EXTRA )

    rewind( FID_CONF )
    read  ( FID_CONF, nml=ANAL, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=ANAL )

    rewind( FID_CONF )
    read  ( FID_CONF, nml=VARI, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=VARI )
    if ( LOUT .and. LOG_DBUG ) write( FID_LOG,* ) ""

    VCOUNT = 0
    do n=1, max_vcount
       if ( trim(VNAME(n)).eq."" ) exit
       VCOUNT = VCOUNT + 1
    enddo
    if ( VCOUNT >= max_vcount ) then
       write (*, *) "ERROR: overflow maximum of vcount"
       call err_abort( 1, __LINE__, loc_main )
    elseif( VCOUNT < 1 ) then
       if ( LOUT ) write( FID_LOG, '(1X,A)' ) "+++ Use Default set of VNAME"
       VCOUNT = num_std_vname
       do n=1, VCOUNT
          VNAME(n) = trim(std_vname(n))
       enddo
    endif

    do n=1, VCOUNT
       VNAME(n) = trim(VNAME(n))
       if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,A)' ) "+++ Listing Vars: (", n, ") ", VNAME(n)
    enddo
    if ( LOUT ) write( FID_LOG,* ) ""

    if ( TARGET_ZLEV(1) >= 0 ) then
       if ( ZCOUNT == 0 ) then
          ZCOUNT = 0
          do n=1, max_zcount
             if ( TARGET_ZLEV(n) < 0 ) exit
             ZCOUNT = ZCOUNT + 1
          enddo
       endif

       Z_LEV_LIST = .true.
       select case( trim(Z_LEV_TYPE) )
       case ( "ZLEV", "zlev" )
          do n=1, ZCOUNT
             if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5,A)' ) "+++ Listing Levs: (", n, ") ",TARGET_ZLEV(n)," [m]"
          enddo
       case ( "PLEV", "plev" )
          do n=1, ZCOUNT
             if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5,A)' ) "+++ Listing Levs: (", n, ") ",TARGET_ZLEV(n)," [hPa]"
          enddo
       case default
          do n=1, ZCOUNT
             if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5,A)' ) "+++ Listing Levs: (", n, ") ",TARGET_ZLEV(n)," [grid]"
          enddo
       end select
    else
       select case( trim(Z_LEV_TYPE) )
       case ( "ZLEV", "zlev" )
          Z_LEV_LIST = .true.
          if ( LOUT ) write( FID_LOG, '(1X,A)' ) "+++ Use Default set of HEGIHT LEVELs"
          ZCOUNT = num_std_zlev
          do n=1, ZCOUNT
             TARGET_ZLEV(n) = std_zlev(n)
             if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5,A)' ) "+++ Listing Levs: (", n, ") ",TARGET_ZLEV(n)," [m]"
          enddo
       case ( "PLEV", "plev" )
          Z_LEV_LIST = .true.
          if ( LOUT ) write( FID_LOG, '(1X,A)' ) "+++ Use Default set of PRESSURE LEVELs"
          ZCOUNT = num_std_plev
          do n=1, ZCOUNT
             TARGET_ZLEV(n) = std_plev(n)
             if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5,A)' ) "+++ Listing Levs: (", n, ") ",TARGET_ZLEV(n)," [hPa]"
          enddo
       case default
          if ( ZCOUNT == 0 ) then
             write (*, *) "ERROR: plz, specify ZCOUNT or TARGET_ZLEV."
             call err_abort( 1, __LINE__, loc_main )
          endif
          Z_LEV_LIST = .false.
          m = ZSTART
          do n=1, ZCOUNT
             TARGET_ZLEV(n) = m
             if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,I5,A)' ) "+++ Listing Levs: (", n, ") ",TARGET_ZLEV(n)," [grid]"
             m = m + 1
          enddo
       end select
    endif
    if ( LOUT ) write( FID_LOG,* ) ""

    rewind( FID_CONF )
    read  ( FID_CONF, nml=GRADS, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=GRADS )

    close ( FID_CONF )


    !--- read run.conf file
    open ( FID_RCNF, file=trim(CONFFILE), status='old', &
           form="formatted", delim='apostrophe', iostat=ierr )
    if ( ierr /= 0 ) then
       write (*, *) "ERROR: fail to open running *.conf file", trim(CONFFILE)
       call err_abort( 1, __LINE__, loc_main )
    endif

    read  ( FID_RCNF, nml=PARAM_TIME, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=PARAM_TIME )

    rewind( FID_RCNF )
    read  ( FID_RCNF, nml=PARAM_PRC, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=PARAM_PRC )

    rewind( FID_RCNF )
    read  ( FID_RCNF, nml=PARAM_HIST, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=PARAM_HIST )

    rewind( FID_RCNF )
    read  ( FID_RCNF, nml=PARAM_HISTORY, iostat=ierr )
    if ( LOUT .and. LOG_DBUG ) write ( FID_LOG, nml=PARAM_HISTORY )

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
       call err_abort( 1, __LINE__, loc_main )
    endif

    return
  end subroutine read_conf


  !> check of zlev compatibility
  !-----------------------------------------------------------------------------------------
  subroutine check_targ_zlev( &
      zz   ) ! [in]
    implicit none

    integer, intent(in)  :: zz
    !---------------------------------------------------------------------------

    select case( vtype )
    case ( vt_urban )
       if ( zz < 1 .or. zz > uz ) then
          write (*, *) "ERROR: requested level is in K-HALO [urban]"
          call err_abort( 1, __LINE__, loc_main )
       endif
    case ( vt_land )
       if ( zz < 1 .or. zz > lz ) then
          write (*, *) "ERROR: requested level is in K-HALO [land]"
          call err_abort( 1, __LINE__, loc_main )
       endif
    case ( vt_height )
       if ( zz < ks .or. zz > ke ) then
           write (*, *) "ERROR: requested level is in K-HALO [height]"
          call err_abort( 1, __LINE__, loc_main )
       endif
    case ( vt_3d )
       if ( zz < ks .or. zz > ke ) then
          write (*, *) "ERROR: requested level is in K-HALO [3D data]"
          call err_abort( 1, __LINE__, loc_main )
       endif
    case ( vt_2d, vt_tpmsk )
       if ( zz /= 1 ) then
          write (*, *) "ERROR: requested level is not exist [2D data]"
          call err_abort( 1, __LINE__, loc_main )
       endif
    case default
        call err_abort( 0, __LINE__, loc_main )
    end select

    return
  end subroutine check_targ_zlev

  !> make vgrid for merged-z output
  !-----------------------------------------------------------------------------------------
  subroutine make_vgrid()
    implicit none

    integer :: iz, ik, k
    real(DP) :: diff, mini
    !---------------------------------------------------------------------------

    allocate( vgrid( ZCOUNT ) )

    if ( atype == a_conv ) then
       do iz = 1, ZCOUNT
          vgrid(iz) = TARGET_ZLEV(iz)
          if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,F8.2)' ) "+++ Target Height: (", iz, ") ", vgrid(iz)
       enddo

    else
       do iz = 1, ZCOUNT
          vgrid(iz) = zlev( TARGET_ZLEV(iz) )
          if ( LOUT ) write( FID_LOG, '(1X,A,I3,A,F8.2,A,I3)' ) &
                      "+++ Target Height: (", iz, ") ", vgrid(iz), " - req. lev = ", TARGET_ZLEV(iz)
       enddo

    endif
    if ( LOUT ) write( FID_LOG, '(1X,A)' ) ""

    return
  end subroutine make_vgrid


  !> allocation of data arrays
  !---------------------------------------------------------------------------
  subroutine allocation( &
      irank   ) ! [in]
    implicit none

    integer, intent(in) :: irank
    !---------------------------------------------------------------------------

    allocate( p_var       (mnxp,       mnyp*nmnge     ) )
    allocate( p_ref       (mnxp,       mnyp,       nz ) )
    allocate( p_cdx       (nxgp*nmnge                 ) )
    allocate( p_cdy       (nygp*nmnge                 ) )
    allocate( p_cx        (nxgp*nmnge                 ) )
    allocate( p_cy        (nygp*nmnge                 ) )
    allocate( cdz         (nz+2*nzh                   ) )
    allocate( cz          (nz+2*nzh                   ) )
    allocate( zlev        (nz                         ) )

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

end module mod_netcdf2grads_h
