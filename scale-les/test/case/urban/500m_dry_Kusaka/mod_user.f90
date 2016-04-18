!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!      Put atmospheric data for urban test
!!      Test is based on Kusaka et al. (2000,BLM)
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

  !-----------------------------------------------------------------------------
contains
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
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

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

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    ! replace urban variables and input to history buffer
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
       TEM00 => CONST_TEM00, &
       PRE00 => CONST_PRE00, &
       LHV   => CONST_LHV,   &    ! ELL : latent heat of vaporization [J/kg]
       Rdry  => CONST_Rdry,  &    ! specific gas constant (dry air)
       CPdry => CONST_CPdry       ! specific heat (dry air,constant pressure) [J/kg/K]
    use scale_grid_real, only: &
       REAL_lon
    use mod_cpl_vars, only: &
       TMPA  => URB_ATM_TEMP,        &  ! air temperature
       PRSA  => URB_ATM_PRES,        &
       WA    => URB_ATM_W,           &
       UA    => URB_ATM_U,           &
       VA    => URB_ATM_V,           &
       RHOA  => URB_ATM_DENS,        &
       QVA   => URB_ATM_QV,          &
       PBL   => URB_ATM_PBL,         &
       PRSS  => URB_ATM_SFC_PRES,    &
       RWD   => URB_ATM_SFLX_rad_dn, &
       RAIN  => URB_ATM_SFLX_rain,   &
       SNOW  => URB_ATM_SFLX_snow
    use scale_time, only:   &
       NOWSEC => TIME_NOWSEC,      & !< subday part  of current time [sec]
       dt_URB => TIME_DTSEC_URBAN    !< time interval of urban step  [sec]
    use mod_admin_time, only: &
       TIME_DOURBAN_step             !< execute urban component in this step?
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP) :: WORK(IA,JA)
    real(RP) :: LON, LAT
    real(RP) :: dsec
    integer  :: tloc

    real(RP) :: LWD(IA,JA)
    real(RP) :: SWD(IA,JA)

    real(RP), parameter :: SRATIO = 0.75_RP ! ratio between direct/total solar [-]
    real(RP)            :: SWtot

    real(RP) :: SW  (0:24)
    real(RP) :: PT  (0:24)
    real(RP) :: Wind(0:24)

    data SW / 0.0,0.0,0.0,0.0,0.0,0.0,50.0,240.0,420.0,600.0,690.0,765.0,800.0, &
              765.0,690.0,600.0,420.0,240.0,50.0,0.0,0.0,0.0,0.0,0.0,0.0/

    data PT / 28.0,27.8,27.65,27.5,27.35,27.2,27.1,27.5,27.85,28.25,28.8,29.4,30.0, &
              30.2,30.4,30.6,30.35,30.1,29.85,29.55,29.15,28.75,28.5,28.25,28.0/

    data Wind / 2.75,2.75,2.75,2.75,2.75,2.75,2.8,3.0,3.25,3.5,3.65,3.65,3.5, &
                3.4,3.27,3.15,3.05,2.95,2.85,2.8,2.75,2.7,2.72,2.75,2.75/

    real(RP) :: THETA
    real(RP) :: RovCP
    integer :: k, i, j
    !---------------------------------------------------------------------------

    RovCP = Rdry / CPdry

    if ( USER_do ) then

       VA  (:,:)        =      0.0_RP
       WA  (:,:)        =      0.0_RP
       RHOA(:,:)        =     1.13_RP
       PBL (:,:)        =    100.0_RP
       RWD (:,:,I_LW,1) =      0.0_RP ! direct
       RWD (:,:,I_LW,2) =    400.0_RP ! diffuse
       PRSA(:,:)        = 100000.0_RP
       PRSS(:,:)        = 100120.0_RP
       QVA (:,:)        =    0.015_RP
       RAIN(:,:)        =      0.0_RP
       SNOW(:,:)        =      0.0_RP

       do j = 1, JA
       do i = 1, IA

          LON = REAL_lon(i,j) / D2R

          tloc = mod(int(NOWSEC/3600.0_RP)+int(LON/15.0_RP),24)

          dsec = mod(NOWSEC,3600.0_RP) / 3600.0_RP

          SWtot = ( ( 1.0_RP-dsec ) * SW(tloc  ) &
                  + (        dsec ) * SW(tloc+1) )

          RWD (i,j,I_SW,1) = (        SRATIO ) * SWtot ! direct
          RWD (i,j,I_SW,2) = ( 1.0_RP-SRATIO ) * SWtot ! diffuse

          THETA     = ( ( 1.0_RP-dsec ) * PT(tloc  ) &
                      + (        dsec ) * PT(tloc+1) ) + TEM00
          WORK(i,j) = THETA                                  ! potential temp
          TMPA(i,j) = THETA * ( PRSA(i,j) / PRE00 )**RovCP   ! air temp, but now PRSA = 100000Pa

          UA  (i,j) = ( ( 1.0_RP-dsec ) * Wind(tloc  ) &
                      + (        dsec ) * Wind(tloc+1) )

          LWD(i,j) = RWD(i,j,I_LW,1) + RWD(i,j,I_LW,2)
          SWD(i,j) = RWD(i,j,I_SW,1) + RWD(i,j,I_SW,2)

       enddo
       enddo

       call HIST_in( WORK (:,:), 'PT_urb',   'Potential air temperature',    'K'     )
       call HIST_in( QVA  (:,:), 'QA_urb',   'Specific humidity',            'kg/kg' )
       call HIST_in( UA   (:,:), 'UA_urb',   'Wind speed',                   'm/s'   )
       call HIST_in( SWD  (:,:), 'SWD_urb',  'Downward shortwave radiation', 'W/m2'  )
       call HIST_in( LWD  (:,:), 'LWD_urb',  'Downward longwave radiation',  'W/m2'  )

       WORK(:,:) = ( RAIN(:,:) + SNOW(:,:) ) * dt_URB
       call HIST_in( WORK (:,:), 'RAIN_urb', 'Precipitation',                'kg/m2' )

    endif

    return
  end subroutine USER_step

end module mod_user
