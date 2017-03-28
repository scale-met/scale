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
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &    ! specific gas constant (dry air)
       CPdry => CONST_CPdry       ! specific heat (dry air,constant pressure) [J/kg/K]
    use scale_time, only:   &
       dt_URB => TIME_DTSEC_URBAN    !< time interval of urban step  [sec]
    use scale_history, only: &
       HIST_in
    use mod_cpl_vars, only: &
       TMPA  => URB_ATM_TEMP,        &
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
    use mod_urban_vars, only: &
       URBAN_TR,         &
       URBAN_TB,         &
       URBAN_TG,         &
       URBAN_TC,         &
       URBAN_QC,         &
       URBAN_UC,         &
       URBAN_TRL,        &
       URBAN_TBL,        &
       URBAN_TGL,        &
       URBAN_ROFF
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP), parameter :: SRATIO = 0.75_RP ! ratio between direct/total solar [-]

    real(RP) :: SX, SB, SG
    real(RP) :: VFGS, VFWS, VFGW, VFWG

    real(RP) :: PTA (IA,JA)
    real(RP) :: LWD (IA,JA)
    real(RP) :: SWD (IA,JA)
    real(RP) :: WORK(IA,JA)
    real(RP) :: RovCP

    integer  :: i, j
    !---------------------------------------------------------------------------

    RovCP = Rdry / CPdry

    if ( USER_do ) then

       UA  (:,:)        = 18.95165494678601_RP
       VA  (:,:)        =               0.0_RP
       WA  (:,:)        =               0.0_RP
       RHOA(:,:)        = 1.193221659609323_RP
       PBL (:,:)        =             100.0_RP
       RWD (:,:,I_LW,1) =               0.0_RP ! direct
       RWD (:,:,I_LW,2) = 434.6034964144717_RP ! diffuse
       PRSA(:,:)        =          100000.0_RP
       PRSS(:,:)        = 102400.6750905938_RP
       QVA (:,:)        = 1.612903266525567E-02_RP
       RAIN(:,:)        =   5.0_RP / 3600.0_RP
       SNOW(:,:)        =               0.0_RP

       URBAN_ROFF(:,:)   = 0.0_RP

       URBAN_TB  (:,:)   = 298.4421947975652_RP
       URBAN_TG  (:,:)   = 298.6575894549702_RP

       URBAN_TC  (:,:)   = 297.7168673290382_RP
       URBAN_QC  (:,:)   = 1.830501916072291E-02_RP
       URBAN_UC  (:,:)   = 7.978626675770874_RP

       URBAN_TBL (1,:,:) = 298.4421947975652_RP
       URBAN_TBL (2,:,:) = 298.9831365435758_RP
       URBAN_TBL (3,:,:) = 299.6502356623731_RP
       URBAN_TBL (4,:,:) = 299.6819173427097_RP
       URBAN_TBL (5,:,:) = 299.0036771970256_RP

       URBAN_TGL (1,:,:) = 298.6575894549702_RP
       URBAN_TGL (2,:,:) = 299.2942430001926_RP
       URBAN_TGL (3,:,:) = 300.0760515131021_RP
       URBAN_TGL (4,:,:) = 300.0731793271447_RP
       URBAN_TGL (5,:,:) = 299.1892611738443_RP

       SB   = 1.429919681745362_RP
       SG   = 2.090724511380252_RP
       VFGS = 0.5335010498145294_RP
       VFWS = 0.3498742126391030_RP
       VFGW = 0.4664989501854706_RP
       VFWG = 0.3498742126391030_RP

       SX = SG / ( VFGS * 0.8_RP +  VFWS * 0.8_RP * 0.2_RP / 0.8_RP * VFGW * 0.8_RP  )
       SX = SB / ( VFWS * 0.8_RP +  VFGS * 0.8_RP * 0.2_RP / 0.8_RP * VFWG * 0.8_RP  )

       RWD(:,:,I_SW,1) = (        SRATIO ) * SX ! direct
       RWD(:,:,I_SW,2) = ( 1.0_RP-SRATIO ) * SX ! diffuse

       PTA (:,:) = 293.7453140572144_RP
       TMPA(:,:) = PTA(:,:) * ( PRSA(:,:) / PRE00 )**RovCP   ! air temp, but now PRSA = 100000Pa

       LWD (:,:) = RWD(:,:,I_LW,1) + RWD(:,:,I_LW,2)
       SWD (:,:) = RWD(:,:,I_SW,1) + RWD(:,:,I_SW,2)

       call HIST_in( PTA (:,:), 'PT_urb',   'Potential air temperature',    'K'     )
       call HIST_in( QVA (:,:), 'QA_urb',   'Specific humidity',            'kg/kg' )
       call HIST_in( UA  (:,:), 'UA_urb',   'Wind speed',                   'm/s'   )
       call HIST_in( SWD (:,:), 'SWD_urb',  'Downward shortwave radiation', 'W/m2'  )
       call HIST_in( LWD (:,:), 'LWD_urb',  'Downward longwave  radiation', 'W/m2'  )
       WORK(:,:) = ( RAIN(:,:) + SNOW(:,:) ) * dt_URB
       call HIST_in( WORK(:,:), 'RAIN_urb', 'Precipitation',                'kg/m2' )

    endif

    return
  end subroutine USER_step

end module mod_user
