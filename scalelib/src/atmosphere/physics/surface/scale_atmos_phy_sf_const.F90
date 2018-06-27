!-------------------------------------------------------------------------------
!> module atmosphere / physics / surface / const
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Constant flux, domain-uniform
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_sf_const
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_const_setup
  public :: ATMOS_PHY_SF_const_flux

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
  integer,  private            :: ATMOS_PHY_SF_FLG_MOM_FLUX = 0 ! application type for momentum flux
                                                                ! 0: Bulk coefficient  is constant
                                                                ! 1: Friction velocity is constant

  real(RP), private, parameter :: ATMOS_PHY_SF_U_maxM      =  100.0_RP ! maximum limit of absolute velocity for momentum [m/s]
  real(RP), private            :: ATMOS_PHY_SF_U_minM      =    0.0_RP ! minimum limit of absolute velocity for momentum [m/s]
  real(RP), private, parameter :: ATMOS_PHY_SF_Cm_max      = 2.5E-3_RP ! maximum limit of bulk coefficient for momentum [NIL]
  real(RP), private            :: ATMOS_PHY_SF_Cm_min      = 1.0E-5_RP ! minimum limit of bulk coefficient for momentum [NIL]

  real(RP), private            :: ATMOS_PHY_SF_Const_Ustar =   0.25_RP ! constant friction velocity [m/s]
  real(RP), private            :: ATMOS_PHY_SF_Const_Cm    = 0.0011_RP ! constant bulk coefficient for momentum [NIL]
  real(RP), private            :: ATMOS_PHY_SF_Const_SH    =   15.0_RP ! constant surface sensible heat flux [W/m2]
  real(RP), private            :: ATMOS_PHY_SF_Const_LH    =  115.0_RP ! constant surface latent   heat flux [W/m2]

  logical,  private            :: ATMOS_PHY_SF_FLG_SH_DIURNAL = .false. ! diurnal modulation for sensible heat flux?
  real(RP), private            :: ATMOS_PHY_SF_Const_FREQ     = 24.0_RP ! frequency of sensible heat flux modulation [hour]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_const_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_ATMOS_PHY_SF_CONST / &
       ATMOS_PHY_SF_FLG_MOM_FLUX,   &
       ATMOS_PHY_SF_U_minM,         &
       ATMOS_PHY_SF_CM_min,         &
       ATMOS_PHY_SF_Const_Ustar,    &
       ATMOS_PHY_SF_Const_Cm,       &
       ATMOS_PHY_SF_Const_SH,       &
       ATMOS_PHY_SF_Const_LH,       &
       ATMOS_PHY_SF_FLG_SH_DIURNAL, &
       ATMOS_PHY_SF_Const_FREQ

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_const_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_SF_const_setup",*) 'Constant flux'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_CONST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_SF_const_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_SF_const_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SF_CONST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SF_CONST)

    return
  end subroutine ATMOS_PHY_SF_const_setup

  !-----------------------------------------------------------------------------
  !> Constant flux
  subroutine ATMOS_PHY_SF_const_flux( &
       IA, IS, IE, JA, JS, JE, &
       ATM_W, ATM_U, ATM_V, ATM_TEMP,               &
       ATM_Z1, SFC_DENS,                            &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH, &
       SFLX_QV,                                     &
       U10, V10                                     )
    use scale_const, only: &
       PI    => CONST_PI
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_time, only: &
       TIME_NOWSEC
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: ATM_W   (IA,JA) ! velocity w  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_U   (IA,JA) ! velocity u  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_V   (IA,JA) ! velocity v  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_TEMP(IA,JA) ! temperature at the lowermost layer (cell center) [K]
    real(RP), intent(in) :: ATM_Z1  (IA,JA) ! height of the lowermost grid from surface (cell center) [m]
    real(RP), intent(in) :: SFC_DENS(IA,JA) ! density     at the surface atmosphere [kg/m3]

    real(RP), intent(out) :: SFLX_MW(IA,JA) ! surface flux for z-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_MU(IA,JA) ! surface flux for x-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_MV(IA,JA) ! surface flux for y-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_SH(IA,JA) ! surface flux for sensible heat (area center)   [J/m2/s]
    real(RP), intent(out) :: SFLX_LH(IA,JA) ! surface flux for latent   heat (area center)   [J/m2/s]
    real(RP), intent(out) :: SFLX_QV(IA,JA) ! surface flux for qv            (area center)   [kg/m2/s]
    real(RP), intent(out) :: U10    (IA,JA) ! velocity u        at 10m height
    real(RP), intent(out) :: V10    (IA,JA) ! velocity v        at 10m height

    real(RP) :: ATM_Uabs(IA,JA) ! absolute velocity at z1 [m/s]

    real(RP) :: Cm(IA,JA)       ! bulk coefficient for momentum
    real(RP) :: R10

    real(RP) :: modulation
    real(RP) :: LHV(IA,JA)

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / surface flux / const'

    !omp parallel do
    do j = JS, JE
    do i = IS, IE
       ATM_Uabs(i,j) = min( ATMOS_PHY_SF_U_maxM, max( ATMOS_PHY_SF_U_minM, &
            sqrt( ATM_W(i,j)**2 + ATM_U(i,j)**2 + ATM_V(i,j)**2 ) ) ) ! at cell center
    enddo
    enddo

    if   ( ATMOS_PHY_SF_FLG_MOM_FLUX == 0 ) then ! Bulk coefficient is constant
       !omp parallel do
       do j = JS, JE
       do i = IS, IE
          Cm(i,j) = ATMOS_PHY_SF_Const_Cm
       enddo
       enddo
    elseif( ATMOS_PHY_SF_FLG_MOM_FLUX == 1 ) then ! Friction velocity is constant
       !omp parallel do
       do j = JS, JE
       do i = IS, IE
          Cm(i,j) = ( ATMOS_PHY_SF_Const_Ustar / ATM_Uabs(i,j) )**2
          Cm(i,j) = min( max( Cm(i,j), ATMOS_PHY_SF_Cm_min ), ATMOS_PHY_SF_Cm_max )
       enddo
       enddo
    endif

    !-----< momentum >-----

    !omp parallel do
    do j = JS, JE
    do i = IS, IE
       SFLX_MW(i,j) = -Cm(i,j) * ATM_Uabs(i,j) * SFC_DENS(i,j) * ATM_W(i,j)
       SFLX_MU(i,j) = -Cm(i,j) * ATM_Uabs(i,j) * SFC_DENS(i,j) * ATM_U(i,j)
       SFLX_MV(i,j) = -Cm(i,j) * ATM_Uabs(i,j) * SFC_DENS(i,j) * ATM_V(i,j)
    enddo
    enddo

    !-----< heat flux >-----

    if ( ATMOS_PHY_SF_FLG_SH_DIURNAL ) then
       modulation = sin( 2.0_RP * PI * TIME_NOWSEC / 3600.0_RP / ATMOS_PHY_SF_Const_FREQ )
    else
       modulation = 1.0_RP
    endif

    !omp parallel do
    do j = JS, JE
    do i = IS, IE
       SFLX_SH(i,j) = ATMOS_PHY_SF_Const_SH * modulation
       SFLX_LH(i,j) = ATMOS_PHY_SF_Const_LH
    enddo
    enddo

    !-----< mass flux >-----
    call HYDROMETEOR_LHV( &
         IA, IS, IE, JA, JS, JE, &
         ATM_TEMP(:,:), & ! [IN]
         LHV(:,:)       ) ! [OUT]

    !omp parallel do
    do j = JS, JE
    do i = IS, IE
       SFLX_QV(i,j) = SFLX_LH(i,j) / LHV(i,j)
    enddo
    enddo

    !-----< U10, V10 >-----

    !omp parallel do
    do j = JS, JE
    do i = IS, IE
       R10 = 10.0_RP / ATM_Z1(i,j)

       U10   (i,j) = R10 * ATM_U(i,j)
       V10   (i,j) = R10 * ATM_V(i,j)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_const_flux

end module scale_atmos_phy_sf_const
