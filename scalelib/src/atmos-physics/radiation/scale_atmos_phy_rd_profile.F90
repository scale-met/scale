!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation / Vertical profile
!!
!! @par Description
!!          Climatological vertical profile (up to mesopause)
!!          for Atmospheric radiation transfer process
!!          this module use following dataset
!!           pressure, temperature         : CIRA86
!!           H2O,CO2,O3,N2O,CO,CH4,O2,CFCs : MIPAS2001
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-06 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_rd_profile
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_PROFILE_setup
  public :: ATMOS_PHY_RD_PROFILE_setup_zgrid
  public :: ATMOS_PHY_RD_PROFILE_read

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: PROFILE_setup_CIRA86
  private :: PROFILE_setup_MIPAS2001
  private :: readfile_MIPAS2001
  private :: PROFILE_read_climatology
  private :: PROFILE_read_CIRA86
  private :: PROFILE_read_MIPAS2001
  private :: PROFILE_read_user
  private :: PROFILE_interp

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP),              private, save :: PROFILE_TOA = 100.0_RP              !< top of atmosphere [km]
  logical,               private, save :: PROFILE_use_climatology = .true.    !< use climatology?
  character(len=H_LONG), private, save :: PROFILE_CIRA86_fname    = "cira.nc" !< file (CIRA86,netCDF format)
  character(len=H_LONG), private, save :: PROFILE_MIPAS2001_dir   = "."       !< dir  (MIPAS2001,ASCII format)
  character(len=H_LONG), private, save :: PROFILE_USER_fname      = ""        !< file (user,ASCII format)
  logical,               private, save :: debug                   = .false.   !< debug mode?

  integer,  private,              save :: CIRA_ntime
  integer,  private,              save :: CIRA_nplev
  integer,  private,              save :: CIRA_nlat
  real(RP), private, allocatable, save :: CIRA_nd  (:)     ! [day]
  real(RP), private, allocatable, save :: CIRA_plog(:)     ! log([hPa])
  real(RP), private, allocatable, save :: CIRA_lat (:)     ! [rad]
  real(RP), private, allocatable, save :: CIRA_temp(:,:,:) ! [K]
  real(RP), private, allocatable, save :: CIRA_z   (:,:,:) ! [km]

  real(RP), private, allocatable, save :: interp_temp(:)   ! [K]
  real(RP), private, allocatable, save :: interp_z   (:)   ! [km]

  integer,  private, parameter :: MIPAS_kmax  = 121
  integer,  private, parameter :: MIPAS_ntime = 2
  real(RP), private, save      :: MIPAS_nd  (0:MIPAS_ntime+1) ! [day]
  real(RP), private, save      :: MIPAS_lat (5)               ! [rad]
  real(RP), private, save      :: MIPAS_z   (MIPAS_kmax,4)    ! [km]
  real(RP), private, save      :: MIPAS_pres(MIPAS_kmax,4)    ! (not used) [hPa]
  real(RP), private, save      :: MIPAS_temp(MIPAS_kmax,4)    ! (not used) [K]
  real(RP), private, save      :: MIPAS_gas (MIPAS_kmax,30,4) ! [ppmv]

  integer,  private, parameter :: I_tropic   =  1
  integer,  private, parameter :: I_midlat   =  2
  integer,  private, parameter :: I_polarsum =  3
  integer,  private, parameter :: I_polarwin =  4

  integer,  private, parameter :: I_N2     =  1
  integer,  private, parameter :: I_O2     =  2
  integer,  private, parameter :: I_CO2    =  3
  integer,  private, parameter :: I_O3     =  4
  integer,  private, parameter :: I_H2O    =  5
  integer,  private, parameter :: I_CH4    =  6
  integer,  private, parameter :: I_N2O    =  7
  integer,  private, parameter :: I_HNO3   =  8
  integer,  private, parameter :: I_CO     =  9
  integer,  private, parameter :: I_NO2    = 10
  integer,  private, parameter :: I_N2O5   = 11
  integer,  private, parameter :: I_ClO    = 12
  integer,  private, parameter :: I_HOCl   = 13
  integer,  private, parameter :: I_ClONO2 = 14
  integer,  private, parameter :: I_NO     = 15
  integer,  private, parameter :: I_HNO4   = 16
  integer,  private, parameter :: I_HCN    = 17
  integer,  private, parameter :: I_NH3    = 18
  integer,  private, parameter :: I_F11    = 19
  integer,  private, parameter :: I_F12    = 20
  integer,  private, parameter :: I_F14    = 21
  integer,  private, parameter :: I_F22    = 22
  integer,  private, parameter :: I_CCl4   = 23
  integer,  private, parameter :: I_COF2   = 24
  integer,  private, parameter :: I_H2O2   = 25
  integer,  private, parameter :: I_C2H2   = 26
  integer,  private, parameter :: I_C2H6   = 27
  integer,  private, parameter :: I_OCS    = 28
  integer,  private, parameter :: I_SO2    = 29
  integer,  private, parameter :: I_SF6    = 30

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_PROFILE_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP)              :: ATMOS_PHY_RD_PROFILE_TOA
    logical               :: ATMOS_PHY_RD_PROFILE_use_climatology
    character(len=H_LONG) :: ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_PROFILE_USER_IN_FILENAME

    namelist / PARAM_ATMOS_PHY_RD_PROFILE / &
       ATMOS_PHY_RD_PROFILE_TOA,                   &
       ATMOS_PHY_RD_PROFILE_use_climatology,       &
       ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME,    &
       ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME, &
       ATMOS_PHY_RD_PROFILE_USER_IN_FILENAME,      &
       debug

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD PROFILE]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ climatological profile'

    ATMOS_PHY_RD_PROFILE_TOA                   = PROFILE_TOA
    ATMOS_PHY_RD_PROFILE_use_climatology       = PROFILE_use_climatology
    ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME    = PROFILE_CIRA86_fname
    ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = PROFILE_MIPAS2001_dir
    ATMOS_PHY_RD_PROFILE_USER_IN_FILENAME      = PROFILE_USER_fname

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_PROFILE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD_PROFILE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_RD_PROFILE)

    PROFILE_TOA             = ATMOS_PHY_RD_PROFILE_TOA
    PROFILE_use_climatology = ATMOS_PHY_RD_PROFILE_use_climatology
    PROFILE_CIRA86_fname    = ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME
    PROFILE_MIPAS2001_dir   = ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME
    PROFILE_USER_fname      = ATMOS_PHY_RD_PROFILE_USER_IN_FILENAME

    if ( PROFILE_use_climatology ) then

       call PROFILE_setup_CIRA86

       call PROFILE_setup_MIPAS2001

    endif

    return
  end subroutine ATMOS_PHY_RD_PROFILE_setup

  !-----------------------------------------------------------------------------
  !> Setup CIRA86 climatological data (temperature, pressure)
  subroutine PROFILE_setup_CIRA86
    use netcdf
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
         CONST_D2R
    use scale_calendar, only: &
         CALENDAR_date2daysec
    implicit none

    integer  :: status, ncid, varid, dimid ! for netCDF

    integer, allocatable :: CIRA_date(:,:)
    integer  :: nday
    real(RP) :: nsec, subsec = 0.0_RP

    real(4), allocatable :: tmp1d(:)
    real(4), allocatable :: tmp3d(:,:,:)

    logical  :: exist
    integer  :: n, m, t
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** read CIRA86 climatology'
    if( IO_L ) write(IO_FID_LOG,*) '*** FILENAME:', trim(PROFILE_CIRA86_fname)

    inquire( file=trim(PROFILE_CIRA86_fname), exist=exist )
    if ( .NOT. exist ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** File not found. check!'
       call PRC_MPIstop
    endif

    ! open CIRA86 datafile (netCDF)
    status = nf90_open( PROFILE_CIRA86_fname, nf90_NoWrite, ncid )
    !if( status /= nf90_noerr ) call handle_err(status)

    status = nf90_inq_dimid( ncid, "time", dimid )
    !if( status /= nf90_NoErr ) call handle_err(status)
    status = nf90_inquire_dimension( ncid, dimid, len=CIRA_ntime )
    !if( status /= nf90_NoErr ) call handle_err(status)

    status = nf90_inq_dimid( ncid, "plev", dimid )
    !if( status /= nf90_NoErr ) call handle_err(status)
    status = nf90_inquire_dimension( ncid, dimid, len=CIRA_nplev )
    !if( status /= nf90_NoErr ) call handle_err(status)

    status = nf90_inq_dimid( ncid, "latitude", dimid )
    !if( status /= nf90_NoErr ) call handle_err(status)
    status = nf90_inquire_dimension( ncid, dimid, len=CIRA_nlat )
    !if( status /= nf90_NoErr ) call handle_err(status)

    !print *, "CIRA_ntime", CIRA_ntime
    !print *, "CIRA_nplev", CIRA_nplev
    !print *, "CIRA_nlat",  CIRA_nlat

    allocate( CIRA_nd(  0:CIRA_ntime+1) )

    allocate( CIRA_plog(CIRA_nplev) )
    allocate( CIRA_lat (CIRA_nlat ) )

    allocate( CIRA_temp(CIRA_nplev,CIRA_nlat,0:CIRA_ntime+1) )
    allocate( CIRA_z   (CIRA_nplev,CIRA_nlat,0:CIRA_ntime+1) )

    ! read time list
    allocate( CIRA_date(6,0:CIRA_ntime+1) )

    CIRA_date(:, 0) = (/ 1985, 12, 15, 12, 0, 0 /)
    CIRA_date(:, 1) = (/ 1986,  1, 15, 12, 0, 0 /)
    CIRA_date(:, 2) = (/ 1986,  2, 15, 12, 0, 0 /)
    CIRA_date(:, 3) = (/ 1986,  3, 15, 12, 0, 0 /)
    CIRA_date(:, 4) = (/ 1986,  4, 15, 12, 0, 0 /)
    CIRA_date(:, 5) = (/ 1986,  5, 15, 12, 0, 0 /)
    CIRA_date(:, 6) = (/ 1986,  6, 15, 12, 0, 0 /)
    CIRA_date(:, 7) = (/ 1986,  7, 15, 12, 0, 0 /)
    CIRA_date(:, 8) = (/ 1986,  8, 15, 12, 0, 0 /)
    CIRA_date(:, 9) = (/ 1986,  9, 15, 12, 0, 0 /)
    CIRA_date(:,10) = (/ 1986, 10, 15, 12, 0, 0 /)
    CIRA_date(:,11) = (/ 1986, 11, 15, 12, 0, 0 /)
    CIRA_date(:,12) = (/ 1986, 12, 15, 12, 0, 0 /)
    CIRA_date(:,13) = (/ 1987,  1, 15, 12, 0, 0 /)

    do t = 0, CIRA_ntime+1
       call CALENDAR_date2daysec( nday, nsec,            & ! [OUT]
                                  CIRA_date(:,t), subsec ) ! [IN]

       CIRA_nd(t) = real(nday,kind=RP) + nsec / 86400.0_RP
    enddo
    deallocate( CIRA_date )

    ! read pressure level [hPa]
    allocate( tmp1d(CIRA_nplev) )

    status = nf90_inq_varid( ncid, "plev", varid )
    !if( status /= nf90_NoErr ) call handle_err(status)
    status = nf90_get_var( ncid, varid, tmp1d(:) )
    !if( status /= nf90_NoErr ) call handle_err(status)

    do n = 1, CIRA_nplev
       CIRA_plog(n) = log( real(tmp1d(n),kind=RP) )
    enddo
    deallocate( tmp1d )

    ! read latitude bin
    allocate( tmp1d(CIRA_nlat) )

    status = nf90_inq_varid( ncid, "latitude", varid )
    !if( status /= nf90_NoErr ) call handle_err(status)
    status = nf90_get_var( ncid, varid, tmp1d(:) )
    !if( status /= nf90_NoErr ) call handle_err(status)

    do n = 1, CIRA_nlat
       CIRA_lat(n) = real(tmp1d(n),kind=RP) * CONST_D2R ! [deg]->[rad]
    enddo
    deallocate( tmp1d )

    !print *, "CIRA_plog", CIRA_plog
    !print *, "CIRA_pres", exp(CIRA_plog)
    !print *, "CIRA_lat", CIRA_lat

    ! read temperature [K]
    allocate( tmp3d(CIRA_nlat,CIRA_nplev,CIRA_ntime) )

    status = nf90_inq_varid( ncid, "ta", varid )
    !if( status /= nf90_NoErr ) call handle_err(status)
    status = nf90_get_var( ncid, varid, tmp3d(:,:,:) )
    !if( status /= nf90_NoErr ) call handle_err(status)

    do m = 1, CIRA_nlat
    do n = 1, CIRA_nplev
    do t = 1, CIRA_ntime
       CIRA_temp(n,m,t) = real(tmp3d(m,n,t),kind=RP)
    enddo
    enddo
    enddo
    ! cyclic condition
    CIRA_temp(:,:,0           ) = CIRA_temp(:,:,CIRA_ntime)
    CIRA_temp(:,:,CIRA_ntime+1) = CIRA_temp(:,:,1         )

    ! avoid missing value
    do m = 1, CIRA_nlat
    do t = 1, CIRA_ntime
    do n = 2, CIRA_nplev
       if( CIRA_temp(n,m,t) >= 999.9_RP ) CIRA_temp(n,m,t) = CIRA_temp(n-1,m,t)
    enddo
    enddo
    enddo

    !print *, "CIRA_temp", CIRA_temp

    ! read geopotencial height [m]
    status = nf90_inq_varid( ncid, "zg", varid )
    !if( status /= nf90_NoErr ) call handle_err(status)
    status = nf90_get_var( ncid, varid, tmp3d(:,:,:) )
    !if( status /= nf90_NoErr ) call handle_err(status)

    do m = 1, CIRA_nlat
    do n = 1, CIRA_nplev
    do t = 1, CIRA_ntime
       CIRA_z(n,m,t) = real(tmp3d(m,n,t),kind=RP) * 1.E-3_RP ! [m]->[km]
    enddo
    enddo
    enddo
    ! cyclic condition
    CIRA_z(:,:,0           ) = CIRA_z(:,:,CIRA_ntime)
    CIRA_z(:,:,CIRA_ntime+1) = CIRA_z(:,:,1         )

    ! avoid missing value
    do m = 1, CIRA_nlat
    do t = 1, CIRA_ntime
    do n = 2, CIRA_nplev
       if( CIRA_z(n,m,t) == 0.999_RP ) CIRA_z(n,m,t) = CIRA_z(n-1,m,t)
    enddo
    enddo
    enddo

    deallocate( tmp3d )

    ! close CIRA86 datafile (netCDF)
    status = nf90_close(ncid)
    !if( status /= nf90_NoErr ) call handle_err(status)

    allocate( interp_temp(CIRA_nplev) )
    allocate( interp_z   (CIRA_nplev) )

    return
  end subroutine PROFILE_setup_CIRA86

  !-----------------------------------------------------------------------------
  !> Setup MIPAS2001 climatological data (gas)
  subroutine PROFILE_setup_MIPAS2001
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_D2R
    use scale_calendar, only: &
       CALENDAR_date2daysec
    implicit none

    character(len=H_LONG) :: fname

    character(len=7), parameter :: MIPAS_fname(4) = (/"equ.atm","day.atm","sum.atm","win.atm"/)

    integer  :: MIPAS_date(6,0:MIPAS_ntime+1)
    integer  :: nday
    real(RP) :: nsec, subsec = 0.0_RP

    character(len=H_LONG) :: dummy
    integer  :: fid, ierr
    integer  :: t, l, rgn
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** read MIPAS2001 gas profile'

    MIPAS_date(:, 0) = (/ 2000, 12, 22, 12, 0, 0 /)
    MIPAS_date(:, 1) = (/ 2001,  6, 21, 12, 0, 0 /)
    MIPAS_date(:, 2) = (/ 2001, 12, 22, 12, 0, 0 /)
    MIPAS_date(:, 3) = (/ 2002,  6, 21, 12, 0, 0 /)

    do t = 0, MIPAS_ntime+1
       call CALENDAR_date2daysec( nday, nsec,             & ! [OUT]
                                  MIPAS_date(:,t), subsec ) ! [IN]

       MIPAS_nd(t) = real(nday,kind=RP) + nsec / 86400.0_RP
    enddo

    MIPAS_lat(1) = -90.0_RP * CONST_D2R
    MIPAS_lat(2) = -45.0_RP * CONST_D2R
    MIPAS_lat(3) =   0.0_RP * CONST_D2R
    MIPAS_lat(4) =  45.0_RP * CONST_D2R
    MIPAS_lat(5) =  90.0_RP * CONST_D2R

    do rgn = I_tropic, I_polarwin
       fname = trim(PROFILE_MIPAS2001_dir)//'/'//MIPAS_fname(rgn)
       if( IO_L ) write(IO_FID_LOG,*) '*** FILENAME:', trim(fname)

       fid = IO_get_available_fid()
       open( unit   = fid,         &
             file   = trim(fname), &
             form   = 'formatted', &
             status = 'old',       &
             iostat = ierr         )

          if ( ierr /= 0 ) then !--- missing
             if( IO_L ) write(IO_FID_LOG,*) '*** File not found. check!'
             call PRC_MPIstop
          endif

          do l = 1, 24
             read(fid,*) dummy
          enddo

          call readfile_MIPAS2001( fid, MIPAS_z   (:,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_pres(:,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_temp(:,rgn) )

          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_N2    ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_O2    ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_CO2   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_O3    ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_H2O   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_CH4   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_N2O   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_HNO3  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_CO    ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_NO2   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_N2O5  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_ClO   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_HOCl  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_ClONO2,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_NO    ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_HNO4  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_HCN   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_NH3   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_F11   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_F12   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_F14   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_F22   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_CCl4  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_COF2  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_H2O2  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_C2H2  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_C2H6  ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_OCS   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_SO2   ,rgn) )
          call readfile_MIPAS2001( fid, MIPAS_gas(:,I_SF6   ,rgn) )
       close(fid)
    enddo

    return
  end subroutine PROFILE_setup_MIPAS2001

  !-----------------------------------------------------------------------------
  !> Setup MIPAS2001 climatological data (gas)
  subroutine readfile_MIPAS2001( &
       fid, &
       var  )
    implicit none

    integer,  intent(in)  :: fid
    real(RP), intent(out) :: var(121)

    character(len=H_LONG) :: dummy
    real(RP)              :: tmp5(5), tmp1

    integer  :: nstr, l
    !---------------------------------------------------------------------------

    read(fid,*) dummy
    !if( IO_L ) write(IO_FID_LOG,*) dummy

    nstr = MIPAS_kmax
    do l = 1, 24
       read(fid,*) tmp5(:)
       !if( IO_L ) write(IO_FID_LOG,*) l, tmp5(:)
       var(nstr  ) = tmp5(1)
       var(nstr-1) = tmp5(2)
       var(nstr-2) = tmp5(3)
       var(nstr-3) = tmp5(4)
       var(nstr-4) = tmp5(5)

       nstr = nstr - 5
    enddo

    read(fid,*) tmp1
    var(nstr) = tmp1

    return
  end subroutine readfile_MIPAS2001

  !-----------------------------------------------------------------------------
  !> Read profile for radiation
  subroutine ATMOS_PHY_RD_PROFILE_read( &
       kmax,         &
       ngas,         &
       ncfc,         &
       naero,        &
       lat,          &
       now_date,     &
       zh,           &
       z,            &
       rhodz,        &
       pres,         &
       presh,        &
       temp,         &
       temph,        &
       gas,          &
       cfc,          &
       aerosol_conc, &
       aerosol_radi, &
       cldfrac       )
    use scale_const, only: &
         GRAV  => CONST_GRAV
    implicit none

    integer,  intent(in)  :: kmax                     !< Number of layer
    integer,  intent(in)  :: ngas                     !< Number of gas species
    integer,  intent(in)  :: ncfc                     !< Number of CFCs
    integer,  intent(in)  :: naero                    !< Number of aerosol(particle) categories
    real(RP), intent(in)  :: lat                      !< latitude [rad]
    integer,  intent(in)  :: now_date(6)              !< date
    real(RP), intent(in)  :: zh(kmax+1)               !< altitude    at the interface [km]
    real(RP), intent(in)  :: z (kmax)                 !< altitude    at the center    [km]
    real(RP), intent(out) :: rhodz       (kmax)       !< density * delta z            [kg/m2]
    real(RP), intent(out) :: pres        (kmax)       !< pressure    at the center    [hPa]
    real(RP), intent(out) :: presh       (kmax+1)     !< pressure    at the interface [hPa]
    real(RP), intent(out) :: temp        (kmax)       !< temperature at the center    [K]
    real(RP), intent(out) :: temph       (kmax+1)     !< temperature at the interface [K]
    real(RP), intent(out) :: gas         (kmax,ngas)  !< gas species   volume mixing ratio [ppmv]
    real(RP), intent(out) :: cfc         (kmax,ncfc)  !< CFCs          volume mixing ratio [ppmv]
    real(RP), intent(out) :: aerosol_conc(kmax,naero) !< cloud/aerosol volume mixing ratio [ppmv]
    real(RP), intent(out) :: aerosol_radi(kmax,naero) !< cloud/aerosol effective radius    [cm]
    real(RP), intent(out) :: cldfrac     (kmax)       !< cloud fraction    [0-1]

    integer  :: k
    !---------------------------------------------------------------------------

    if ( PROFILE_use_climatology ) then

       call PROFILE_read_climatology( kmax,     & ! [IN]
                                      ngas,     & ! [IN]
                                      ncfc,     & ! [IN]
                                      naero,    & ! [IN]
                                      lat,      & ! [IN], tentative treatment
                                      now_date, & ! [IN]
                                      zh,       & ! [IN]
                                      z,        & ! [IN]
                                      pres,     & ! [OUT]
                                      presh,    & ! [OUT]
                                      temp,     & ! [OUT]
                                      temph,    & ! [OUT]
                                      gas,      & ! [OUT]
                                      cfc       ) ! [OUT]

    else

       call PROFILE_read_user( kmax,  & ! [IN]
                               ngas,  & ! [IN]
                               ncfc,  & ! [IN]
                               naero, & ! [IN]
                               zh,    & ! [IN]
                               z,     & ! [IN]
                               pres,  & ! [OUT]
                               presh, & ! [OUT]
                               temp,  & ! [OUT]
                               temph, & ! [OUT]
                               gas,   & ! [OUT]
                               cfc    ) ! [OUT]

    endif


    do k = 1, kmax
       rhodz(k) = ( presh(k+1) - presh(k) ) * 100.0_RP / GRAV
    enddo

    ! no cloud/aerosol
    aerosol_conc(:,:) = 0.0_RP
    aerosol_radi(:,:) = 0.0_RP
    cldfrac     (:)   = 0.0_RP

    !----- report data -----
    if ( debug ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|===============      Vertical  Coordinate      =================|'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|            -GRID CENTER-             -GRID INTERFACE-          |'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|    k       z      pres    temp       zh      pres    temp    k |'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|         [km]     [hPa]     [K]     [km]     [hPa]     [K]      |'
       k = 1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,F10.4,F8.2,I5,A)') &
       '|                                ',zh(k),presh(k),temph(k),k, ' | TOA'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,F10.4,F8.2,A)') &
       '|',k,z(k),pres(k),temp(k),    '                                 | '
       do k = 2, kmax-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,F10.4,F8.2,I5,A)') &
       '|                                ',zh(k),presh(k),temph(k),k, ' | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,F10.4,F8.2,A)') &
       '|',k,z(k),pres(k),temp(k),    '                                 | '
       enddo
       k = kmax
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,F10.4,F8.2,I5,A)') &
       '|                                ',zh(k),presh(k),temph(k),k, ' | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,F10.4,F8.2,A)') &
       '|',k,z(k),pres(k),temp(k),    '                                 | '
       k = kmax+1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,F10.4,F8.2,I5,A)') &
       '|                                ',zh(k),presh(k),temph(k),k, ' | Ground'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|================================================================|'

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|=====================================================================================|'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|                                         -Gas concetrations [ppmv]-                  |'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|    k        z       H2O       CO2        O3       N2O        CO       CH4        O2 |'
       do k = 1, kmax
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,1F9.3,7E10.3,A)') '|',k,z(k),gas(k,:),' | '
       enddo
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|=====================================================================================|'

    endif

    return
  end subroutine ATMOS_PHY_RD_PROFILE_read

  !-----------------------------------------------------------------------------
  !> Read climatological profile from file
  subroutine PROFILE_read_climatology( &
       kmax,         &
       ngas,         &
       ncfc,         &
       naero,        &
       lat,          &
       now_date,     &
       zh,           &
       z,            &
       pres,         &
       presh,        &
       temp,         &
       temph,        &
       gas,          &
       cfc           )
    implicit none

    integer,  intent(in)  :: kmax             !< Number of layer
    integer,  intent(in)  :: ngas             !< Number of gas species
    integer,  intent(in)  :: ncfc             !< Number of CFCs
    integer,  intent(in)  :: naero            !< Number of aerosol(particle) categories
    real(RP), intent(in)  :: lat              !< latitude [rad]
    integer,  intent(in)  :: now_date(6)      !< date
    real(RP), intent(in)  :: zh   (kmax+1)    !< altitude    at the interface [km]
    real(RP), intent(in)  :: z    (kmax)      !< altitude    at the center    [km]
    real(RP), intent(out) :: pres (kmax)      !< pressure    at the center    [hPa]
    real(RP), intent(out) :: presh(kmax+1)    !< pressure    at the interface [hPa]
    real(RP), intent(out) :: temp (kmax)      !< temperature at the center    [K]
    real(RP), intent(out) :: temph(kmax+1)    !< temperature at the interface [K]
    real(RP), intent(out) :: gas  (kmax,ngas) !< gas species   volume mixing ratio [ppmv]
    real(RP), intent(out) :: cfc  (kmax,ncfc) !< CFCs          volume mixing ratio [ppmv]
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** [RD_PROFILE] generate climatological profile'

    call PROFILE_read_CIRA86( kmax,        & ! [IN]
                              lat,         & ! [IN]
                              now_date(:), & ! [IN]
                              zh      (:), & ! [IN]
                              z       (:), & ! [IN]
                              presh   (:), & ! [OUT]
                              temph   (:), & ! [OUT]
                              pres    (:), & ! [OUT]
                              temp    (:)  ) ! [OUT]

    call PROFILE_read_MIPAS2001( kmax,          & ! [IN]
                                 ngas,          & ! [IN]
                                 ncfc,          & ! [IN]
                                 lat,           & ! [IN]
                                 now_date(:),   & ! [IN]
                                 z       (:),   & ! [IN]
                                 gas     (:,:), & ! [OUT]
                                 cfc     (:,:)  ) ! [OUT]

    ! no vapor
    gas(:,1) = 0.0_RP

    return
  end subroutine PROFILE_read_climatology

  !-----------------------------------------------------------------------------
  !> interpolate CIRA86 climatological data (temperature, pressure)
  subroutine PROFILE_read_CIRA86( &
       kmax,     &
       lat,      &
       now_date, &
       zh,       &
       z,        &
       presh,    &
       temph,    &
       pres,     &
       temp      )
    use scale_calendar, only: &
         CALENDAR_date2daysec
    implicit none

    integer,  intent(in)  :: kmax          !< Number of layer
    real(RP), intent(in)  :: lat           !< latitude [rad]
    integer,  intent(in)  :: now_date(6)   !< date
    real(RP), intent(in)  :: zh   (kmax+1) !< altitude    at the interface [km]
    real(RP), intent(in)  :: z    (kmax)   !< altitude    at the center    [km]
    real(RP), intent(out) :: presh(kmax+1) !< pressure    at the interface [hPa]
    real(RP), intent(out) :: temph(kmax+1) !< temperature at the interface [K]
    real(RP), intent(out) :: pres (kmax)   !< pressure    at the center    [hPa]
    real(RP), intent(out) :: temp (kmax)   !< temperature at the center    [K]

    real(RP) :: plogh(kmax+1)
    real(RP) :: plog (kmax)

    integer  :: now_date_mod(6), nday
    real(RP) :: nd, nsec, subsec = 0.0_RP

    integer  :: nplev_mod
    integer  :: indexLAT, indexD
    real(RP) :: factLAT, factD

    integer  :: n, t
    !---------------------------------------------------------------------------

    ! latitude interpolation
    if    ( lat <  CIRA_lat(1) ) then ! extrapolation
       indexLAT = 1
       factLAT  = 0.0_RP
    elseif( lat >= CIRA_lat(CIRA_nlat) ) then ! extrapolation
       indexLAT = CIRA_nlat - 1
       factLAT  = 1.0_RP
    else
       do n = 1, CIRA_nlat-1
          if (       lat >= CIRA_lat(n  ) &
               .AND. lat <  CIRA_lat(n+1) ) then ! interpolation
             indexLAT = n
             factLAT  = ( lat-CIRA_lat(n) ) / ( CIRA_lat(n+1)-CIRA_lat(n) )
          endif
       enddo
    endif

    ! time interpolation
    now_date_mod(2:6) = now_date(2:6)
    now_date_mod(1)   = 1986

    call CALENDAR_date2daysec( nday, nsec,             & ! [OUT]
                               now_date_mod(:), subsec ) ! [IN]

    nd = real(nday,kind=RP) + nsec / 86400.0_RP

    do t = 0, CIRA_ntime
       if (       nd >= CIRA_nd(t  ) &
            .AND. nd <  CIRA_nd(t+1) ) then ! interpolation
             indexD = t
             factD  = ( nd-CIRA_nd(t) ) / ( CIRA_nd(t+1)-CIRA_nd(t) )
       endif
    enddo

    interp_z(:) = CIRA_z(:,indexLAT  ,indexD  ) * ( 1.0_RP-factLAT ) * ( 1.0_RP-factD ) &
                + CIRA_z(:,indexLAT+1,indexD  ) * (        factLAT ) * ( 1.0_RP-factD ) &
                + CIRA_z(:,indexLAT  ,indexD+1) * ( 1.0_RP-factLAT ) * (        factD ) &
                + CIRA_z(:,indexLAT+1,indexD+1) * (        factLAT ) * (        factD )

    interp_temp(:) = CIRA_temp(:,indexLAT  ,indexD  ) * ( 1.0_RP-factLAT ) * ( 1.0_RP-factD ) &
                   + CIRA_temp(:,indexLAT+1,indexD  ) * (        factLAT ) * ( 1.0_RP-factD ) &
                   + CIRA_temp(:,indexLAT  ,indexD+1) * ( 1.0_RP-factLAT ) * (        factD ) &
                   + CIRA_temp(:,indexLAT+1,indexD+1) * (        factLAT ) * (        factD )

    ! avoid missing value
    nplev_mod = CIRA_nplev
    do n = CIRA_nplev, 1, -1
       if ( interp_temp(n) == interp_temp(n-1) ) then
          nplev_mod = nplev_mod-1
       else
          exit
       endif
    enddo

    !do n = 1, nplev_mod
    !   print *, n, interp_z(n), interp_temp(n), indexLAT, indexD, factLAT, factD
    !enddo

    ! pressure interpolation
    call PROFILE_interp( nplev_mod,              & ! [IN]
                         interp_z (1:nplev_mod), & ! [IN]
                         CIRA_plog(1:nplev_mod), & ! [IN]
                         kmax+1,                 & ! [IN]
                         zh(:),                  & ! [IN]
                         plogh(:)                ) ! [OUT]

    presh(:) = exp( plogh(:) )

    call PROFILE_interp( kmax+1, zh(:), plogh(:), kmax, z(:), plog(:) )
    pres (:) = exp( plog (:) )

    call PROFILE_interp( nplev_mod,                & ! [IN]
                         interp_z   (1:nplev_mod), & ! [IN]
                         interp_temp(1:nplev_mod), & ! [IN]
                         kmax,                     & ! [IN]
                         z(:),                     & ! [IN]
                         temp(:)                   ) ! [OUT]

    call PROFILE_interp( nplev_mod,                & ! [IN]
                         interp_z   (1:nplev_mod), & ! [IN]
                         interp_temp(1:nplev_mod), & ! [IN]
                         kmax+1,                   & ! [IN]
                         zh(:),                    & ! [IN]
                         temph(:)                  ) ! [OUT]

    return
  end subroutine PROFILE_read_CIRA86

  !-----------------------------------------------------------------------------
  !> interpolate MIPAS2001 climatological data (gas)
  subroutine PROFILE_read_MIPAS2001( &
       kmax,     &
       ngas,     &
       ncfc,     &
       lat,      &
       now_date, &
       z,        &
       gas,      &
       cfc       )
    use scale_calendar, only: &
         CALENDAR_date2daysec
    implicit none

    integer,  intent(in)    :: kmax           !< Number of layer
    integer,  intent(in)    :: ngas           !< Number of gas species
    integer,  intent(in)    :: ncfc           !< Number of CFCs
    real(RP), intent(in)    :: lat            !< latitude [rad]
    integer,  intent(in)    :: now_date(6)    !< date
    real(RP), intent(in)    :: z  (kmax)      !< altitude at the center [km]
    real(RP), intent(inout) :: gas(kmax,ngas) !< gas species   volume mixing ratio [ppmv]
    real(RP), intent(inout) :: cfc(kmax,ncfc) !< CFCs          volume mixing ratio [ppmv]

    real(RP) :: interp_gas(MIPAS_kmax,30)
    real(RP) :: interp_z  (MIPAS_kmax)

    integer  :: now_date_mod(6), nday
    real(RP) :: nd, nsec, subsec = 0.0_RP

    integer  :: indexD1, indexD2
    real(RP) :: factLAT, factD
    !---------------------------------------------------------------------------

    ! time interpolation
    now_date_mod(2:6) = now_date(2:6)
    now_date_mod(1)   = 2001

    call CALENDAR_date2daysec( nday, nsec,             & ! [OUT]
                               now_date_mod(:), subsec ) ! [IN]

    nd = real(nday,kind=RP) + nsec / 86400.0_RP

    if    ( nd >= MIPAS_nd(0) .AND. nd < MIPAS_nd(1) ) then ! winter-summer

       indexD1 = I_polarwin
       indexD2 = I_polarsum
       factD   = ( nd-MIPAS_nd(0) ) / ( MIPAS_nd(1)-MIPAS_nd(0) )

    elseif( nd >= MIPAS_nd(1) .AND. nd < MIPAS_nd(2) ) then ! summer-winter

       indexD1 = I_polarsum
       indexD2 = I_polarwin
       factD   = ( nd-MIPAS_nd(1) ) / ( MIPAS_nd(2)-MIPAS_nd(1) )

    elseif( nd >= MIPAS_nd(2) .AND. nd < MIPAS_nd(3) ) then ! winter-summer

       indexD1 = I_polarwin
       indexD2 = I_polarsum
       factD   = ( nd-MIPAS_nd(2) ) / ( MIPAS_nd(3)-MIPAS_nd(2) )

    endif

    ! latitude interpolation
    if    (                           lat < MIPAS_lat(1) ) then ! south pole

       interp_gas(:,:) = MIPAS_gas(:,:,indexD1 ) * ( 1.0_RP-factD ) &
                       + MIPAS_gas(:,:,indexD2 ) * (        factD )

       interp_z(:) = MIPAS_z(:,indexD1 ) * ( 1.0_RP-factD ) &
                   + MIPAS_z(:,indexD2 ) * (        factD )

    elseif( lat >= MIPAS_lat(1) .AND. lat < MIPAS_lat(2) ) then ! south pole-SH mid

       factLAT = ( lat-MIPAS_lat(1) ) / ( MIPAS_lat(2)-MIPAS_lat(1) )

       interp_gas(:,:) = MIPAS_gas(:,:,indexD1) * ( 1.0_RP-factD ) * ( 1.0_RP-factLAT ) &
                       + MIPAS_gas(:,:,indexD2) * (        factD ) * ( 1.0_RP-factLAT ) &
                       + MIPAS_gas(:,:,I_midlat)                   * (        factLAT )

       interp_z(:) = MIPAS_z(:,indexD1) * ( 1.0_RP-factD ) * ( 1.0_RP-factLAT ) &
                   + MIPAS_z(:,indexD2) * (        factD ) * ( 1.0_RP-factLAT ) &
                   + MIPAS_z(:,I_midlat)                   * (        factLAT )

    elseif( lat >= MIPAS_lat(2) .AND. lat < MIPAS_lat(3) ) then ! SH mid-EQ

       factLAT = ( lat-MIPAS_lat(2) ) / ( MIPAS_lat(3)-MIPAS_lat(2) )

       interp_gas(:,:) = MIPAS_gas(:,:,I_midlat) * ( 1.0_RP-factLAT ) &
                       + MIPAS_gas(:,:,I_tropic) * (        factLAT )

       interp_z(:) = MIPAS_z(:,I_midlat) * ( 1.0_RP-factLAT ) &
                   + MIPAS_z(:,I_tropic) * (        factLAT )

    elseif( lat >= MIPAS_lat(3) .AND. lat < MIPAS_lat(4) ) then ! EQ-NH mid

       factLAT = ( lat-MIPAS_lat(3) ) / ( MIPAS_lat(4)-MIPAS_lat(3) )

       interp_gas(:,:) = MIPAS_gas(:,:,I_tropic) * ( 1.0_RP-factLAT ) &
                       + MIPAS_gas(:,:,I_midlat) * (        factLAT )

       interp_z(:) = MIPAS_z(:,I_tropic) * ( 1.0_RP-factLAT ) &
                   + MIPAS_z(:,I_midlat) * (        factLAT )

    elseif( lat >= MIPAS_lat(4) .AND. lat < MIPAS_lat(5) ) then ! NH mid-north pole

       factLAT = ( lat-MIPAS_lat(4) ) / ( MIPAS_lat(5)-MIPAS_lat(4) )

       interp_gas(:,:) = MIPAS_gas(:,:,I_midlat)                   * ( 1.0_RP-factLAT ) &
                       + MIPAS_gas(:,:,indexD2) * ( 1.0_RP-factD ) * (        factLAT ) &
                       + MIPAS_gas(:,:,indexD1) * (        factD ) * (        factLAT )

       interp_z(:) = MIPAS_z(:,I_midlat)                   * ( 1.0_RP-factLAT ) &
                   + MIPAS_z(:,indexD2) * ( 1.0_RP-factD ) * (        factLAT ) &
                   + MIPAS_z(:,indexD1) * (        factD ) * (        factLAT )

    elseif( lat >= MIPAS_lat(5)                          ) then ! north pole

       interp_gas(:,:) = MIPAS_gas(:,:,indexD2) * ( 1.0_RP-factD ) &
                       + MIPAS_gas(:,:,indexD1) * (        factD )

       interp_z(:) = MIPAS_z(:,indexD2) * ( 1.0_RP-factD ) &
                   + MIPAS_z(:,indexD1) * (        factD )

    endif

    gas(:,:) = 0.0_RP
    cfc(:,:) = 0.0_RP

    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_H2O   ), kmax, z(:), gas(:,1) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_CO2   ), kmax, z(:), gas(:,2) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_O3    ), kmax, z(:), gas(:,3) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_N2O   ), kmax, z(:), gas(:,4) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_CO    ), kmax, z(:), gas(:,5) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_CH4   ), kmax, z(:), gas(:,6) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_O2    ), kmax, z(:), gas(:,7) )

    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_F11   ), kmax, z(:), cfc(:, 1) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_F12   ), kmax, z(:), cfc(:, 2) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_F14   ), kmax, z(:), cfc(:, 4) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_F22   ), kmax, z(:), cfc(:, 9) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_SF6   ), kmax, z(:), cfc(:,22) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_ClONO2), kmax, z(:), cfc(:,23) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_CCl4  ), kmax, z(:), cfc(:,24) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_N2O5  ), kmax, z(:), cfc(:,25) )
    call PROFILE_interp( MIPAS_kmax, interp_z(:), interp_gas(:,I_HNO4  ), kmax, z(:), cfc(:,27) )

    return
  end subroutine PROFILE_read_MIPAS2001

  !-----------------------------------------------------------------------------
  !> Linear interpolation & extrapolation
  subroutine PROFILE_interp( imax1, x1, y1, imax2, x2, y2 )
    implicit none

    integer,  intent(in)  :: imax1
    real(RP), intent(in)  :: x1(imax1)
    real(RP), intent(in)  :: y1(imax1)
    integer,  intent(in)  :: imax2
    real(RP), intent(in)  :: x2(imax2)
    real(RP), intent(out) :: y2(imax2)

    real(RP) :: fact

    integer  :: i1, i2
    !---------------------------------------------------------------------------

    do i2 = 1, imax2

       if ( x2(i2) > x1(1) ) then ! extrapolation

             fact = ( x1(1) - x2(i2) ) / ( x1(2) - x1(1) )

             y2(i2) = y1(1) * ( 1.0_RP-fact ) &
                    + y1(2) * (        fact )

       elseif( x2(i2) <= x1(imax1) ) then ! extrapolation

             fact = ( x1(imax1) - x2(i2) ) / ( x1(imax1) - x1(imax1-1) )

             y2(i2) = y1(imax1-1) * (        fact ) &
                    + y1(imax1  ) * ( 1.0_RP-fact )

       else
          do i1 = 1, imax1-1
             if (       x2(i2) <= x1(i1  ) &
                  .AND. x2(i2) >  x1(i1+1) ) then ! interpolation

                fact = ( x2(i2) - x1(i1) ) / ( x1(i1+1) - x1(i1) )

                y2(i2) = y1(i1  ) * ( 1.0_RP-fact ) &
                       + y1(i1+1) * (        fact )

                exit
             endif
          enddo
       endif

    enddo

    return
  end subroutine PROFILE_interp

  !-----------------------------------------------------------------------------
  !> Setup vertical grid for radiation
  subroutine ATMOS_PHY_RD_PROFILE_setup_zgrid( &
       kmax, &
       kadd, &
       zh,   &
       z     )
    use scale_grid, only: &
       CZ  => GRID_CZ,  &
       FZ  => GRID_FZ
    implicit none

    integer,  intent(in)    :: kmax       !< number of vertical grid
    integer,  intent(in)    :: kadd       !< number of additional vertical grid
    real(RP), intent(inout) :: zh(kmax+1) !< altitude at the interface [km]
    real(RP), intent(inout) :: z (kmax)   !< altitude at the center    [km]

    real(RP) :: dz
    integer  :: k, RD_k
    !---------------------------------------------------------------------------

    if ( kadd > 0 ) then
       !--- additional layer over the computational domain
       dz = ( PROFILE_TOA - FZ(KE)*1.E-3_RP ) / real( kadd, kind=RP )

       zh(1) = PROFILE_TOA
       do k = 2, kadd
          zh(k) = zh(k-1) - dz
       enddo
       zh(kadd+1) = FZ(KE)*1.E-3_RP

       ! linear interpolation for center of layers
       do k = 1, kadd
          z(k) = 0.5_RP * ( zh(k+1) + zh(k) )
       enddo
    endif

    !--- in computational domain [NOTE] layer index is reversed
    do k = KS-1, KE
       RD_k = kmax - ( k - KS )
       zh(RD_k) = FZ(k)*1.E-3_RP
    enddo
    do k = KS, KE
       RD_k = kmax - ( k - KS )
       z(RD_k) = CZ(k)*1.E-3_RP
    enddo

    !----- report data -----
    if ( debug ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|===============  Vertical Coordinate  ===============|'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|      -GRID CENTER-            -GRID INTERFACE-      |'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '| RD_k       z    k GRID_CZ GRID_FZ    k      zh RD_k |'
       if ( kadd > 0 ) then
       RD_k = 1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,A)') &
       '|                                       ',zh(RD_k),RD_k,    ' | TOA'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,A)') &
       '|',RD_k,z(RD_k),     '                                        | '
       do RD_k = 2, kadd-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,A)') &
       '|                                       ',zh(RD_k),RD_k,    ' | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,A)') &
       '|',RD_k,z(RD_k),     '                                        | '
       enddo
       RD_k = kadd
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,A)') &
       '|                                       ',zh(RD_k),RD_k,    ' | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,A)') &
       '|',RD_k,z(RD_k),     '                                        | KADD'
       RD_k = kadd+1
       k    = kmax - RD_k + KS
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,F8.3,I5,A)') &
       '|                          ',FZ(k)*1.E-3_RP,k,zh(RD_k),RD_k,' | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,I5,F8.3,A)') &
       '|',RD_k,z(RD_k),k,CZ(k)*1.E-3_RP, '                           | KADD+1=KE'
       else
       RD_k = 1
       k    = kmax - RD_k + KS
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,F8.3,I5,A)') &
       '|                          ',FZ(k)*1.E-3_RP,k,zh(RD_k),RD_k,' | TOA=KE'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,I5,F8.3,A)') &
       '|',RD_k,z(RD_k),k,CZ(k)*1.E-3_RP, '                           | '
       endif
       do RD_k = kadd+2, kmax-1
       k    = kmax - RD_k + KS
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,F8.3,I5,A)') &
       '|                          ',FZ(k)*1.E-3_RP,k,zh(RD_k),RD_k,' | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,I5,F8.3,A)') &
       '|',RD_k,z(RD_k),k,CZ(k)*1.E-3_RP, '                           | '
       enddo
       RD_k = kmax
       k    = kmax - RD_k + KS
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,F8.3,I5,A)') &
       '|                          ',FZ(k),k,zh(RD_k),RD_k, ' | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I5,F8.3,I5,F8.3,A)') &
       '|',RD_k,z(RD_k),k,CZ(k)*1.E-3_RP, '                           | RD_KMAX=KS'
       RD_k = kmax+1
       k    = kmax - RD_k + KS
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,I5,F8.3,I5,A)') &
       '|                          ',FZ(k)*1.E-3_RP,k,zh(RD_k),RD_k,' | Ground'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)')                   &
       '|=====================================================|'

    endif

    return
  end subroutine ATMOS_PHY_RD_PROFILE_setup_zgrid

  !-----------------------------------------------------------------------------
  !> Read user-defined profile from file
  subroutine PROFILE_read_user( &
       kmax,         &
       ngas,         &
       ncfc,         &
       naero,        &
       zh,           &
       z,            &
       pres,         &
       presh,        &
       temp,         &
       temph,        &
       gas,          &
       cfc           )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       Mdry => CONST_Mdry, &
       Mvap => CONST_Mvap, &
       PPM  => CONST_PPM
    implicit none

    integer,  intent(in)  :: kmax             !< Number of layer
    integer,  intent(in)  :: ngas             !< Number of gas species
    integer,  intent(in)  :: ncfc             !< Number of CFCs
    integer,  intent(in)  :: naero            !< Number of aerosol(particle) categories
    real(RP), intent(in)  :: zh   (kmax+1)    !< altitude    at the interface [km]
    real(RP), intent(in)  :: z    (kmax)      !< altitude    at the center    [km]
    real(RP), intent(out) :: pres (kmax)      !< pressure    at the center    [hPa]
    real(RP), intent(out) :: presh(kmax+1)    !< pressure    at the interface [hPa]
    real(RP), intent(out) :: temp (kmax)      !< temperature at the center    [K]
    real(RP), intent(out) :: temph(kmax+1)    !< temperature at the interface [K]
    real(RP), intent(out) :: gas  (kmax,ngas) !< gas species   volume mixing ratio [ppmv]
    real(RP), intent(out) :: cfc  (kmax,ncfc) !< CFCs          volume mixing ratio [ppmv]

    integer, parameter :: USER_klim = 500
    integer  :: USER_kmax
    real(RP) :: USER_z   (USER_klim)
    real(RP) :: USER_pres(USER_klim)
    real(RP) :: USER_temp(USER_klim)
    real(RP) :: USER_qv  (USER_klim)
    real(RP) :: USER_o3  (USER_klim)

    real(RP), allocatable :: work_z(:)
    real(RP), allocatable :: work  (:)

    real(RP) :: plog (kmax)
    real(RP) :: plogh(kmax+1)

    character(len=H_LONG) :: dummy

    integer  :: fid, ierr
    integer  :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** [RD_PROFILE] user-defined profile'

    gas(:,:) = 0.0_RP
    cfc(:,:) = 0.0_RP

    if( IO_L ) write(IO_FID_LOG,*) '*** FILENAME:', trim(PROFILE_USER_fname)

    fid = IO_get_available_fid()
    open( unit   = fid,                      &
          file   = trim(PROFILE_USER_fname), &
          form   = 'formatted',              &
          status = 'old',                    &
          iostat = ierr                      )

       if ( ierr /= 0 ) then !--- missing
          if( IO_L ) write(IO_FID_LOG,*) '*** File not found. check!'
          call PRC_MPIstop
       endif

       read(fid,*) dummy

       do k = 1, USER_klim
          read(fid,*,iostat=ierr) USER_z(k), USER_pres(k), USER_temp(k), USER_qv(k), USER_o3(k)
          if ( ierr /= 0 ) exit
       enddo
       USER_kmax = k - 1
    close(fid)

    allocate( work_z(USER_kmax) )
    allocate( work  (USER_kmax) )

    do k = 1, USER_kmax
       work_z(k) = USER_z(k) / 1000.0_RP ! [m->km]
       work  (k) = log( USER_pres(k)/100.0_RP ) ! log[Pa->hPA]
    enddo

    ! pressure interpolation
    call PROFILE_interp( USER_kmax, & ! [IN]
                         work_z(:), & ! [IN]
                         work  (:), & ! [IN]
                         kmax+1,    & ! [IN]
                         zh(:),     & ! [IN]
                         plogh(:)   ) ! [OUT]

    presh(:) = exp( plogh(:) )

    call PROFILE_interp( kmax+1, zh(:), plogh(:), kmax, z(:), plog(:) )
    pres (:) = exp( plog (:) )

    do k = 1, USER_kmax
       work(k) = USER_temp(k)
    enddo

    call PROFILE_interp( USER_kmax, & ! [IN]
                         work_z(:), & ! [IN]
                         work  (:), & ! [IN]
                         kmax,      & ! [IN]
                         z(:),      & ! [IN]
                         temp(:)    ) ! [OUT]

    call PROFILE_interp( USER_kmax, & ! [IN]
                         work_z(:), & ! [IN]
                         work  (:), & ! [IN]
                         kmax+1,    & ! [IN]
                         zh(:),     & ! [IN]
                         temph(:)   ) ! [OUT]

    do k = 1, USER_kmax
       work(k) = USER_qv(k) / Mvap * Mdry / PPM ! [kg/kg->PPM]
    enddo

    call PROFILE_interp( USER_kmax, & ! [IN]
                         work_z(:), & ! [IN]
                         work  (:), & ! [IN]
                         kmax,      & ! [IN]
                         z(:),      & ! [IN]
                         gas(:,1)   ) ! [OUT]

    do k = 1, USER_kmax
       work(k) = USER_o3(k) / 48.0_RP * Mdry / PPM ! [kg/kg->PPM]
    enddo

    call PROFILE_interp( USER_kmax, & ! [IN]
                         work_z(:), & ! [IN]
                         work  (:), & ! [IN]
                         kmax,      & ! [IN]
                         z(:),      & ! [IN]
                         gas(:,3)   ) ! [OUT]

    return
  end subroutine PROFILE_read_user

end module scale_atmos_phy_rd_profile
!-------------------------------------------------------------------------------
