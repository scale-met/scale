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
#include "scalelib.h"
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_finalize
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  !> Tracer setup
  subroutine USER_tracer_setup
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
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

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'
    LOG_INFO("USER_setup",*) 'User procedure in test/case/urban/offline'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    LOG_NEWLINE

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

    call USER_update

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Finalization
  subroutine USER_finalize
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_finalize

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine USER_calc_tendency
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_update
    use scale_const, only: &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &    ! specific gas constant (dry air)
       CPdry => CONST_CPdry       ! specific heat (dry air,constant pressure) [J/kg/K]
    use scale_atmos_hydrometeor, only: &
       CV_WATER
    use scale_time, only:   &
       dt_URB => TIME_DTSEC_URBAN    !< time interval of urban step  [sec]
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_cpl_vars, only: &
       TMPA  => URB_ATM_TEMP,        &
       PRSA  => URB_ATM_PRES,        &
       WA    => URB_ATM_W,           &
       UA    => URB_ATM_U,           &
       VA    => URB_ATM_V,           &
       RHOA  => URB_ATM_DENS,        &
       QVA   => URB_ATM_QV,          &
       PBL   => URB_ATM_PBL,         &
       RHOS  => URB_ATM_SFC_DENS,    &
       PRSS  => URB_ATM_SFC_PRES,    &
       RWD   => URB_ATM_SFLX_rad_dn, &
       PREC  => URB_ATM_SFLX_water,  &
       ENGI  => URB_ATM_SFLX_ENGI
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
    use scale_file_history, only: &
       FILE_HISTORY_in
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
       RWD (:,:,I_R_direct ,I_R_IR) =               0.0_RP ! direct
       RWD (:,:,I_R_diffuse,I_R_IR) = 434.6034964144717_RP ! diffuse
       PRSA(:,:)        =          100000.0_RP
       PRSS(:,:)        = 102400.6750905938_RP
       QVA (:,:)        = 1.612903266525567E-02_RP
       PREC(:,:)        =   5.0_RP / 3600.0_RP

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

       RWD(:,:,I_R_direct ,I_R_NIR) = 0.0_RP
       RWD(:,:,I_R_diffuse,I_R_NIR) = 0.0_RP
       RWD(:,:,I_R_direct ,I_R_VIS) = (        SRATIO ) * SX ! direct
       RWD(:,:,I_R_diffuse,I_R_VIS) = ( 1.0_RP-SRATIO ) * SX ! diffuse

       PTA (:,:) = 293.7453140572144_RP
       TMPA(:,:) = PTA(:,:) * ( PRSA(:,:) / PRE00 )**RovCP   ! air temp, but now PRSA = 100000Pa

       RHOS(:,:) = PRSS(:,:) / ( Rdry * TMPA(:,:) )

       ENGI(:,:) = PREC(:,:) * CV_WATER * TMPA(:,:)

       LWD (:,:) = RWD(:,:,I_R_direct ,I_R_IR) + RWD(:,:,I_R_direct ,I_R_VIS)
       SWD (:,:) = RWD(:,:,I_R_diffuse,I_R_IR) + RWD(:,:,I_R_diffuse,I_R_VIS)

       call FILE_HISTORY_in( PTA (:,:), 'PT_urb',   'Potential air temperature',    'K'     )
       call FILE_HISTORY_in( QVA (:,:), 'QA_urb',   'Specific humidity',            'kg/kg' )
       call FILE_HISTORY_in( UA  (:,:), 'UA_urb',   'Wind speed',                   'm/s'   )
       call FILE_HISTORY_in( SWD (:,:), 'SWD_urb',  'Downward shortwave radiation', 'W/m2'  )
       call FILE_HISTORY_in( LWD (:,:), 'LWD_urb',  'Downward longwave  radiation', 'W/m2'  )
       WORK(:,:) = PREC(:,:) * dt_URB
       call FILE_HISTORY_in( WORK(:,:), 'RAIN_urb', 'Precipitation',                'kg/m2' )

    endif

    return
  end subroutine USER_update

end module mod_user
