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
    LOG_INFO("USER_setup",*) 'User procedure in test/case/urban/500m_wet_Aoyagi'

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
       D2R   => CONST_D2R,   &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &    ! specific gas constant (dry air)
       CPdry => CONST_CPdry       ! specific heat (dry air,constant pressure) [J/kg/K]
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_lon
    use scale_time, only:   &
       NOWSEC => TIME_NOWSEC,      & !< subday part  of current time [sec]
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
    use scale_atmos_hydrometeor, only: &
       CV_WATER
    implicit none

    real(RP) :: SW    (0:24)
    real(RP) :: PT    (0:24)
    real(RP) :: Wind  (0:24)
    real(RP) :: Prcp  (0:24)
    real(RP) :: Qvapor(0:24) ! mixing ratio [kg/kg]

    data SW /       0.0_RP, & ! 00LST
                    0.0_RP, & ! 01LST
                    0.0_RP, & ! 02LST
                    0.0_RP, & ! 03LST
                    0.0_RP, & ! 04LST
                   11.0_RP, & ! 05LST
                   26.5_RP, & ! 06LST
                  161.5_RP, & ! 07LST
                  434.2_RP, & ! 08LST
                  694.6_RP, & ! 09LST
                  877.5_RP, & ! 10LST
                  807.4_RP, & ! 11LST
                  907.4_RP, & ! 12LST
                  847.2_RP, & ! 13LST
                  725.2_RP, & ! 14LST
                  559.0_RP, & ! 15LST
                  357.9_RP, & ! 16LST
                  158.6_RP, & ! 17LST
                   14.7_RP, & ! 18LST
                    0.0_RP, & ! 19LST
                    0.0_RP, & ! 20LST
                    0.0_RP, & ! 21LST
                    0.0_RP, & ! 22LST
                    0.0_RP, & ! 23LST
                    0.0_RP  / ! 24LST

    data PT /     296.6_RP, & ! 00LST
                  296.6_RP, & ! 01LST
                  296.7_RP, & ! 02LST
                  296.5_RP, & ! 03LST
                  296.4_RP, & ! 04LST
                  296.3_RP, & ! 05LST
                  295.9_RP, & ! 06LST
                  296.1_RP, & ! 07LST
                  296.7_RP, & ! 08LST
                  298.0_RP, & ! 09LST
                  299.2_RP, & ! 10LST
                  300.0_RP, & ! 11LST
                  300.6_RP, & ! 12LST
                  301.7_RP, & ! 13LST
                  302.2_RP, & ! 14LST
                  302.4_RP, & ! 15LST
                  302.0_RP, & ! 16LST
                  301.1_RP, & ! 17LST
                  300.5_RP, & ! 18LST
                  299.8_RP, & ! 19LST
                  299.1_RP, & ! 20LST
                  298.9_RP, & ! 21LST
                  298.6_RP, & ! 22LST
                  298.3_RP, & ! 23LST
                  296.6_RP  / ! 24LST

    data Wind /     3.4_RP, & ! 00LST
                    3.2_RP, & ! 01LST
                    3.2_RP, & ! 02LST
                    3.2_RP, & ! 03LST
                    3.2_RP, & ! 04LST
                    3.0_RP, & ! 05LST
                    3.8_RP, & ! 06LST
                    3.6_RP, & ! 07LST
                    3.5_RP, & ! 08LST
                    3.6_RP, & ! 09LST
                    4.0_RP, & ! 10LST
                    3.8_RP, & ! 11LST
                    4.0_RP, & ! 12LST
                    3.9_RP, & ! 13LST
                    4.0_RP, & ! 14LST
                    4.1_RP, & ! 15LST
                    3.7_RP, & ! 16LST
                    3.6_RP, & ! 17LST
                    4.2_RP, & ! 18LST
                    3.9_RP, & ! 19LST
                    3.7_RP, & ! 20LST
                    3.4_RP, & ! 21LST
                    3.8_RP, & ! 22LST
                    4.0_RP, & ! 23LST
                    3.4_RP  / ! 24LST

    data Qvapor / 0.018_RP, & ! 00LST
                  0.018_RP, & ! 01LST
                  0.018_RP, & ! 02LST
                  0.018_RP, & ! 03LST
                  0.018_RP, & ! 04LST
                  0.018_RP, & ! 05LST
                  0.017_RP, & ! 06LST
                  0.017_RP, & ! 07LST
                  0.017_RP, & ! 08LST
                  0.018_RP, & ! 09LST
                  0.018_RP, & ! 10LST
                  0.018_RP, & ! 11LST
                  0.018_RP, & ! 12LST
                  0.017_RP, & ! 13LST
                  0.017_RP, & ! 14LST
                  0.016_RP, & ! 15LST
                  0.017_RP, & ! 16LST
                  0.018_RP, & ! 17LST
                  0.018_RP, & ! 18LST
                  0.017_RP, & ! 19LST
                  0.018_RP, & ! 20LST
                  0.018_RP, & ! 21LST
                  0.018_RP, & ! 22LST
                  0.018_RP, & ! 23LST
                  0.018_RP  / ! 24LST

    data Prcp /     0.0_RP, & ! 00LST
                    0.0_RP, & ! 01LST
                    0.0_RP, & ! 02LST
                    0.0_RP, & ! 03LST
                    0.1_RP, & ! 04LST
                    0.2_RP, & ! 05LST
                    3.0_RP, & ! 06LST
                    0.4_RP, & ! 07LST
                    0.0_RP, & ! 08LST
                    0.0_RP, & ! 09LST
                    0.0_RP, & ! 10LST
                    0.0_RP, & ! 11LST
                    0.0_RP, & ! 12LST
                    0.0_RP, & ! 13LST
                    0.0_RP, & ! 14LST
                    0.0_RP, & ! 15LST
                    0.0_RP, & ! 16LST
                    0.0_RP, & ! 17LST
                    0.0_RP, & ! 18LST
                    0.0_RP, & ! 19LST
                    0.0_RP, & ! 20LST
                    0.0_RP, & ! 21LST
                    0.0_RP, & ! 22LST
                    0.0_RP, & ! 23LST
                    0.0_RP  / ! 24LST

    real(RP), parameter :: SRATIO = 0.75_RP ! ratio between direct/total solar [-]

    real(RP) :: PTA (IA,JA)
    real(RP) :: LWD (IA,JA)
    real(RP) :: SWD (IA,JA)
    real(RP) :: WORK(IA,JA)
    real(RP) :: RovCP

    real(RP) :: SWtot
    real(RP) :: LON
    real(RP) :: dsec
    integer  :: tloc

    integer  :: i, j
    !---------------------------------------------------------------------------

    RovCP = Rdry / CPdry

    if ( USER_do ) then

       VA  (:,:)        =      0.0_RP
       WA  (:,:)        =      0.0_RP
       RHOA(:,:)        =     1.13_RP
       PBL (:,:)        =    100.0_RP
       RWD (:,:,I_R_direct ,I_R_IR) =   0.0_RP ! direct
       RWD (:,:,I_R_diffuse,I_R_IR) = 400.0_RP ! diffuse
       PRSA(:,:)        = 100000.0_RP
       PRSS(:,:)        = 100120.0_RP

       do j = 1, JA
       do i = 1, IA

          LON = ATMOS_GRID_CARTESC_REAL_lon(i,j) / D2R

          tloc = mod(int(NOWSEC/3600.0_RP)+int(LON/15.0_RP),24)

          dsec = mod(NOWSEC,3600.0_RP) / 3600.0_RP

          SWtot = ( ( 1.0_RP-dsec ) * SW(tloc  ) &
                  + (        dsec ) * SW(tloc+1) )

          RWD (i,j,I_R_direct ,I_R_NIR) = 0.0_RP
          RWD (i,j,I_R_diffuse,I_R_NIR) = 0.0_RP
          RWD (i,j,I_R_direct ,I_R_VIS) = (        SRATIO ) * SWtot ! direct
          RWD (i,j,I_R_diffuse,I_R_VIS) = ( 1.0_RP-SRATIO ) * SWtot ! diffuse

          PTA (i,j) = ( ( 1.0_RP-dsec ) * PT(tloc  ) &
                      + (        dsec ) * PT(tloc+1) )
          TMPA(i,j) = PTA(i,j) * ( PRSA(i,j) / PRE00 )**RovCP   ! air temp, but now PRSA = 100000Pa

          RHOS(i,j) = PRSS(i,j) / ( Rdry * TMPA(i,j) )

          UA  (i,j) = ( ( 1.0_RP-dsec ) * Wind(tloc  ) &
                      + (        dsec ) * Wind(tloc+1) )

          LWD (i,j) = RWD(i,j,I_R_direct ,I_R_IR) + RWD(i,j,I_R_direct ,I_R_VIS)
          SWD (i,j) = RWD(i,j,I_R_diffuse,I_R_IR) + RWD(i,j,I_R_diffuse,I_R_VIS)

          QVA (i,j) = ( ( 1.0_RP-dsec ) * Qvapor(tloc  ) &
                      + (        dsec ) * Qvapor(tloc+1) )
          QVA (i,j) = QVA(i,j) / (1.0_RP + QVA(i,j))         ! [mixing ratio->specific humidity]

          PREC(i,j) = ( ( 1.0_RP-dsec ) * Prcp(tloc  ) &
                      + (        dsec ) * Prcp(tloc+1) ) / 3600.0_RP ! [mm/h->kg/m2/s]

          ENGI(i,j) = PREC(i,j) * CV_WATER * TMPA(i,j)

!          if ( with_rain ) then
!             if (( tloc >= 1 ).and.( tloc < 10 )) then
!                RAIN(i,j) = 5.0_RP / 3600.0_RP
!             endif
!                RAIN(i,j) = 1.0_RP / 3600.0_RP
!          endif

       enddo
       enddo

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
