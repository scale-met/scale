!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!      Put atmospheric data for urban test
!!      Test is based on Aoyagi et al. (2011,JAMC)
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
  logical, private :: USER_do   = .false. !< do user step?

  integer, parameter, private :: data_length_max = 10000    ! maximum data length
  real(RP)          , private :: SNOWIN (data_length_max)
  real(RP)          , private :: TAIN   (data_length_max)
  real(RP)          , private :: RHIN   (data_length_max)
  real(RP)          , private :: WINDIN (data_length_max)
  real(RP)          , private :: SHORTIN(data_length_max)
  real(RP)          , private :: LONGIN (data_length_max)
  integer           , private :: stepnum
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
    implicit none

    namelist / PARAM_USER / &
       USER_do

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


    call read_input_atm_data
    stepnum = 0

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
       NOWSEC => TIME_NOWSEC,      & !< subday part  of current time [sec]
       dt_URB => TIME_DTSEC_URBAN    !< time interval of urban step  [sec]
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

    !real(RP) :: PTA (IA,JA)
    real(RP) :: LWD (IA,JA)
    real(RP) :: SWD (IA,JA)
    real(RP) :: WORK(IA,JA)
    !real(RP) :: RovCP
    real(RP) :: RH
    real(RP) :: QAsat

    !real(RP) :: SWtot
    real(RP) :: LON
    real(RP) :: dsec
    integer  :: tloc

    integer  :: i, j
    !---------------------------------------------------------------------------

    !RovCP = Rdry / CPdry

    if ( USER_do ) then

       stepnum = stepnum + 1     ! step number of

       TMPA(:,:)         = TAIN(stepnum)
       PRSA(:,:)         = 100000.0_RP
       PRSS(:,:)         = 100120.0_RP
       UA  (:,:)         = WINDIN(stepnum)
       VA  (:,:)         = 0.0_RP
       WA  (:,:)         = 0.0_RP
       RHOA(:,:)         = 1.13_RP
       RH                = RHIN(stepnum)*0.01
       PBL (:,:)         = 100.0_RP
       RWD (:,:,I_SW,1)  = SHORTIN(stepnum) ! direct
       RWD (:,:,I_SW,2)  = 0.0_RP ! duffusion
       RWD (:,:,I_LW,1)  = 0.0_RP ! direct
       RWD (:,:,I_LW,2)  = LONGIN(stepnum)
       RAIN(:,:)         = 0.0_RP
       SNOW(:,:)         = SNOWIN(stepnum)


       do j = 1, JA
       do i = 1, IA

          LON = REAL_lon(i,j) / D2R

          tloc = mod(int(NOWSEC/3600.0_RP)+int(LON/15.0_RP),24)

          dsec = mod(NOWSEC,3600.0_RP) / 3600.0_RP

          !SWtot = ( ( 1.0_RP-dsec ) * SW(tloc  ) &
          !        + (        dsec ) * SW(tloc+1) )

          !RWD (i,j,I_SW,1) = (        SRATIO ) * SWtot ! direct
          !RWD (i,j,I_SW,2) = ( 1.0_RP-SRATIO ) * SWtot ! diffuse

          !PTA (i,j) = ( ( 1.0_RP-dsec ) * PT(tloc  ) &
          !            + (        dsec ) * PT(tloc+1) )
          !TMPA(i,j) = PTA(i,j) * ( PRSA(i,j) / PRE00 )**RovCP   ! air temp, but now PRSA = 100000Pa

          RHOS(i,j) = PRSS(i,j) / ( Rdry * TMPA(i,j) )

          !UA  (i,j) = ( ( 1.0_RP-dsec ) * Wind(tloc  ) &
          !            + (        dsec ) * Wind(tloc+1) )

          LWD (i,j) = RWD(i,j,I_LW,1) + RWD(i,j,I_LW,2)
          SWD (i,j) = RWD(i,j,I_SW,1) + RWD(i,j,I_SW,2)

          !QVA (i,j) = ( ( 1.0_RP-dsec ) * Qvapor(tloc  ) &
          !            + (        dsec ) * Qvapor(tloc+1) )
          !QVA (i,j) = QVA(i,j) / (1.0_RP + QVA(i,j))         ! [mixing ratio->specific humidity]

          call qsatf( QAsat, TMPA(i,j), PRSA(i,j) )
          QVA  (i,j) =  QAsat* RH

          !RAIN(i,j) = ( ( 1.0_RP-dsec ) * Prcp(tloc  ) &
          !            + (        dsec ) * Prcp(tloc+1) ) / 3600.0_RP ! [mm/h->kg/m2/s]

          SNOW(i,j) = SNOW(i,j) / 3600.0_RP ! [mm/h->kg/m2/s]

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
  subroutine read_input_atm_data
    implicit none

    integer :: timestep
    integer :: PREProws, TArows, RHrows, WINDrows, SHORTrows, LONGrows

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

