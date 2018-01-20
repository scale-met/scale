!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!      Put atmospheric data for offline (between snow/land and atmos) test
!!      Input data is provided from observation at JMA Sapporo
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private            :: USER_do   = .false. !< do user step?
  real(DP), private            :: INPUT_UPDATE_DT = 0.0_DP
  integer,  private            :: INPUT_UPDATE_NSTEP
  integer,  private            :: nowstep
  integer,  private            :: totalstep
  integer,  private            :: stepnum, stepnum1,stepnum2

  integer,  private            :: LAND_NSTEP

  integer,            private :: data_length
  integer, parameter, private :: data_length_max = 10000    ! maximum data length
  real(RP)          , private :: SNOWIN (data_length_max)
  real(RP)          , private :: TAIN   (data_length_max)
  real(RP)          , private :: RHIN   (data_length_max)
  real(RP)          , private :: WINDIN (data_length_max)
  real(RP)          , private :: SHORTIN(data_length_max)
  real(RP)          , private :: LONGIN (data_length_max)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_atmos_admin, only: &
       ATMOS_do,        &
       ATMOS_sw_dyn,    &
       ATMOS_sw_phy_mp, &
       ATMOS_sw_phy_ae, &
       ATMOS_sw_phy_ch, &
       ATMOS_sw_phy_rd, &
       ATMOS_sw_phy_sf, &
       ATMOS_sw_phy_tb, &
       ATMOS_sw_phy_cp
    use scale_time, only:   &
       dt_LND => TIME_DTSEC_LAND,  & !< time interval of land step  [sec]
       TIME_NSTEP,                 & !< total steps
       TIME_DSTEP_LAND               !< step interval of land step
    implicit none

    namelist / PARAM_USER / &
       USER_do,             &
       INPUT_UPDATE_DT

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    if( IO_L ) write(IO_FID_LOG,*)

    ! atmosphric model set to off
    ATMOS_do        = .false.
    ATMOS_sw_dyn    = .false.
    ATMOS_sw_phy_mp = .false.
    ATMOS_sw_phy_ae = .false.
    ATMOS_sw_phy_ch = .false.
    ATMOS_sw_phy_rd = .false.
    ATMOS_sw_phy_sf = .false.
    ATMOS_sw_phy_tb = .false.
    ATMOS_sw_phy_cp = .false.

    call read_input_atm_data(data_length)
    nowstep   = -1
    totalstep = 0
    stepnum   = 1

    LAND_NSTEP = TIME_NSTEP / TIME_DSTEP_LAND + 1 ! +1 is initial

    if ( INPUT_UPDATE_DT <= 0.0_DP ) then
       write(*,*) 'xxx You need specify value larger than 0.0 to INPUT_UPDATE_DT'
       call PRC_MPIstop
    endif
    INPUT_UPDATE_NSTEP = nint( INPUT_UPDATE_DT / dt_LND )
    if ( abs(INPUT_UPDATE_NSTEP * dt_LND - INPUT_UPDATE_DT) > 1E-10_DP ) then
       write(*,*) 'xxx INPUT_UPDATE_DT is not multiple of LAND DT'
       call PRC_MPIstop
    end if

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    call USER_step

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       I_LW  => CONST_I_LW,  &
       I_SW  => CONST_I_SW,  &
       D2R   => CONST_D2R,   &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &    ! specific gas constant (dry air)
       CPdry => CONST_CPdry       ! specific heat (dry air,constant pressure) [J/kg/K]
    use scale_atmos_saturation, only:  &
       qsatf => ATMOS_SATURATION_pres2qsat_all  ! better to  change name from qsatf to qsat
    use scale_grid_real, only: &
       REAL_lon
    use scale_time, only:   &
       NOWSEC => TIME_NOWSEC,      & ! subday part  of current time [sec]
       dt_LND => TIME_DTSEC_LAND     ! time interval of land step  [sec]
    use scale_history, only: &
       HIST_in
    use mod_cpl_vars, only: &
       TMPA  => LND_ATM_TEMP,        &
       PRSA  => LND_ATM_PRES,        &
       WA    => LND_ATM_W,           &
       UA    => LND_ATM_U,           &
       VA    => LND_ATM_V,           &
       RHOA  => LND_ATM_DENS,        &
       QVA   => LND_ATM_QV,          &
       PBL   => LND_ATM_PBL,         &
       RHOS  => LND_ATM_SFC_DENS,    &
       PRSS  => LND_ATM_SFC_PRES,    &
       RWD   => LND_ATM_SFLX_rad_dn, &
       RAIN  => LND_ATM_SFLX_rain,   &
       SNOW  => LND_ATM_SFLX_snow
    implicit none

    real(RP) :: LWD (IA,JA)
    real(RP) :: SWD (IA,JA)
    real(RP) :: WORK(IA,JA)
    real(RP) :: RH
    real(RP) :: QAsat

    !real(RP) :: LON
    !real(RP) :: dsec
    real(RP) :: fact

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( USER_do ) then

       ! counter of total step of land model
       totalstep = totalstep + 1

       ! counter for reading input data
       nowstep = nowstep + 1
       if( nowstep >= INPUT_UPDATE_NSTEP )then
          nowstep = 0
          stepnum = stepnum + 1
       endif
       fact = real(nowstep, kind=RP) / real(INPUT_UPDATE_NSTEP, kind=RP)

       if( totalstep == LAND_NSTEP )then ! last
          stepnum1 = stepnum
          stepnum2 = stepnum
       else
          stepnum1 = stepnum
          stepnum2 = stepnum + 1
       endif
       if( max(stepnum1,stepnum2) > data_length )then
          write(*,*) "xxx Number of input data is not enough: ",data_length,"<",max(stepnum1,stepnum2)
          call PRC_MPIstop
       endif

       TMPA(:,:)         = ( (1.0_RP-fact) * TAIN(stepnum1)  &
                           +         fact  * TAIN(stepnum2) )
       PRSA(:,:)         = 100000.0_RP
       PRSS(:,:)         = 100120.0_RP
       UA  (:,:)         = ( (1.0_RP-fact) * WINDIN(stepnum1) &
                           +         fact  * WINDIN(stepnum2) )
       VA  (:,:)         = 0.0_RP
       WA  (:,:)         = 0.0_RP
       RHOA(:,:)         = 1.13_RP
       RH                = ( (1.0_RP-fact) * RHIN(stepnum1) &
                           +         fact  * RHIN(stepnum2) )
       PBL (:,:)         = 100.0_RP
       RWD (:,:,I_SW,1)  = ( (1.0_RP-fact) * SHORTIN(stepnum1) &
                           +         fact  * SHORTIN(stepnum2) ) ! direct
       RWD (:,:,I_SW,2)  = 0.0_RP                                 ! duffusion
       RWD (:,:,I_LW,1)  = 0.0_RP
       RWD (:,:,I_LW,2)  = ( (1.0_RP-fact) * LONGIN(stepnum1) &
                           +         fact  * LONGIN(stepnum2) )
       SNOW(:,:)         = ( (1.0_RP-fact) * SNOWIN(stepnum1) &
                           +         fact  * SNOWIN(stepnum2) )
       RAIN(:,:)         = 0.0_RP


       do j = 1, JA
       do i = 1, IA

          !LON = REAL_lon(i,j) / D2R
          !tloc = mod(int(NOWSEC/3600.0_RP)+int(LON/15.0_RP),24)
          !dsec = mod(NOWSEC,3600.0_RP) / 3600.0_RP

          RHOS(i,j) = PRSS(i,j) / ( Rdry * TMPA(i,j) )
          LWD (i,j) = RWD(i,j,I_LW,1) + RWD(i,j,I_LW,2)
          SWD (i,j) = RWD(i,j,I_SW,1) + RWD(i,j,I_SW,2)

          call qsatf( QAsat, TMPA(i,j), PRSA(i,j) )
          QVA (i,j) = QAsat * (RH * 0.01)
          SNOW(i,j) = SNOW(i,j) / real(INPUT_UPDATE_DT, kind=RP) ! [mm/h->kg/m2/s]

       enddo
       enddo

       call HIST_in( TMPA (:,:), 'TA_snow',   'Potential air temperature',    'K'     )
       call HIST_in( QVA  (:,:), 'QA_snow',   'Specific humidity',            'kg/kg' )
       call HIST_in( UA   (:,:), 'UA_snow',   'Wind speed',                   'm/s'   )
       call HIST_in( SWD  (:,:), 'SWD_snow',  'Downward shortwave radiation', 'W/m2'  )
       call HIST_in( LWD  (:,:), 'LWD_snow',  'Downward longwave  radiation', 'W/m2'  )
       !WORK(:,:) = ( RAIN (:,:) + SNOW(:,:) ) * dt_URB
       !call HIST_in( WORK (:,:), 'RAIN_snow', 'Precipitation',                'kg/m2' )
       call HIST_in( RAIN (:,:), 'RAIN_snow', 'Rainfall',                     'kg/m2/s' )
       call HIST_in( SNOW (:,:), 'SNOW_snow', 'Snowfall',                     'kg/m2/s' )

    endif

    return
  end subroutine USER_step

!----------------------------------------------------------
  subroutine read_input_atm_data(timestep)
    implicit none

    integer, intent(out) :: timestep
    integer              :: PREProws, TArows, RHrows, WINDrows, SHORTrows, LONGrows

    call inputtext('./input/SNOW_input.txt',    SNOWIN,  PREProws )
    call inputtext('./input/TEMPAIR_input.txt', TAIN,    TArows   )
    call inputtext('./input/RH_input.txt',      RHIN,    RHrows   )
    call inputtext('./input/WIND_input.txt',    WINDIN,  WINDrows )
    call inputtext('./input/SW_input.txt',      SHORTIN, SHORTrows)
    call inputtext('./input/LW_input.txt',      LONGIN,  LONGrows )

    write(*,*) "SNOWIN  =", SNOWIN(1), SNOWIN(PREProws),  PREProws
    write(*,*) "TAIN    =", TAIN(1),   TAIN(TArows),      TArows
    write(*,*) "RHIN    =", RHIN(1),   RHIN(RHrows),      RHrows
    write(*,*) "WINDIN  =", WINDIN(1), WINDIN(WINDrows),  WINDrows
    write(*,*) "SHORTIN =", SHORTIN(1),SHORTIN(SHORTrows),SHORTrows
    write(*,*) "LONGIN  =", LONGIN(1), LONGIN(LONGrows),  LONGrows

    timestep = min(PREProws, TArows, RHrows, WINDrows, SHORTrows, LONGrows)
    write(*,*) "Timestep = ",timestep

    return
  end subroutine read_input_atm_data

  subroutine inputtext (filename, data, rows)
    implicit none
    character(len=*), intent(in)  :: filename
    real(RP),         intent(out) :: data(data_length_max)
    integer,          intent(out) :: rows

    integer :: IO_FID_SNOW_TEST
    integer :: iymd,ihh

    integer :: i, io

    IO_FID_SNOW_TEST = IO_get_available_fid()
    open(IO_FID_SNOW_TEST, file = trim(filename), status = 'OLD')
    do i=1,data_length_max
       read(IO_FID_SNOW_TEST,*,end=100) iymd,ihh,data(i)
       rows = i
    end do
100 close(IO_FID_SNOW_TEST)

    !if ( rows > data_length_max ) then
    !   write(*,*) "data array is not enough. Please check. ", trim(filename)
    !   write(*,*) "data_length_max = ",data_length_max, "row number= ",rows
    !endif

    return
  end subroutine inputtext

 !----------------------------------------------

end module mod_user

