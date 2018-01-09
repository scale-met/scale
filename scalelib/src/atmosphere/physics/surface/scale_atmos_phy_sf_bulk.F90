!-------------------------------------------------------------------------------
!> module atmosphere / physics / surface / bulk
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Bulk Method
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_sf_bulk
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_bulk_setup
  public :: ATMOS_PHY_SF_bulk_flux

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
  real(RP), private            :: ATMOS_PHY_SF_beta   =   1.0_RP ! evaporation efficiency (0-1)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_bulk_setup
    use scale_process, only: &
       PRC_abort
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_SF_BULK / &
       ATMOS_PHY_SF_beta

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[surface bulk] / Categ[atmosphere physics] / Origin[SCALE lib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Bulk scheme'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_BULK,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_BULK. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_SF_BULK)

    return
  end subroutine ATMOS_PHY_SF_bulk_setup

  !-----------------------------------------------------------------------------
  !> Calculate surface flux
  subroutine ATMOS_PHY_SF_bulk_flux( &
       IA, IS, IE, JA, JS, JE, &
       ATM_W, ATM_U, ATM_V,          &
       ATM_TEMP, ATM_PRES, ATM_QV,   &
       SFC_DENS, SFC_TEMP, SFC_PRES, &
       SFC_Z0M, SFC_Z0H, SFC_Z0E,    &
       PBL, ATM_Z1,                  &
       SFLX_MW, SFLX_MU, SFLX_MV,    &
       SFLX_SH, SFLX_LH, SFLX_QV,    &
       U10, V10, T2, Q2              )
    use scale_const, only: &
       PRE00 => CONST_PRE00, &
       CPdry => CONST_CPdry, &
       Rdry  => CONST_Rdry
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       I_QV
    use scale_atmos_saturation, only: &
       SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
       BULKFLUX
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: ATM_W   (IA,JA) ! velocity w  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_U   (IA,JA) ! velocity u  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_V   (IA,JA) ! velocity v  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in) :: ATM_TEMP(IA,JA) ! temperature at the lowermost layer (cell center) [K]
    real(RP), intent(in) :: ATM_PRES(IA,JA) ! pressure    at the lowermost layer (cell center) [Pa]
    real(RP), intent(in) :: ATM_QV  (IA,JA) ! qv          at the lowermost layer (cell center) [kg/kg]
    real(RP), intent(in) :: SFC_DENS(IA,JA) ! density     at the surface atmosphere [kg/m3]
    real(RP), intent(in) :: SFC_TEMP(IA,JA) ! temperature at the surface skin [K]
    real(RP), intent(in) :: SFC_PRES(IA,JA) ! pressure    at the surface atmosphere [Pa]
    real(RP), intent(in) :: SFC_Z0M (IA,JA) ! surface roughness length (momentum) [m]
    real(RP), intent(in) :: SFC_Z0H (IA,JA) ! surface roughness length (heat) [m]
    real(RP), intent(in) :: SFC_Z0E (IA,JA) ! surface roughness length (vapor) [m]
    real(RP), intent(in) :: PBL     (IA,JA) ! depth of the PBL [m]
    real(RP), intent(in) :: ATM_Z1  (IA,JA) ! height of the lowermost grid from surface (cell center) [m]

    real(RP), intent(out) :: SFLX_MW(IA,JA) ! surface flux for z-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_MU(IA,JA) ! surface flux for x-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_MV(IA,JA) ! surface flux for y-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out) :: SFLX_SH(IA,JA) ! surface flux for sensible heat (area center)   [J/m2/s]
    real(RP), intent(out) :: SFLX_LH(IA,JA) ! surface flux for latent   heat (area center)   [J/m2/s]
    real(RP), intent(out) :: SFLX_QV(IA,JA) ! surface flux for qv            (area center)   [kg/m2/s]
    real(RP), intent(out) :: U10    (IA,JA) ! velocity u        at 10m height
    real(RP), intent(out) :: V10    (IA,JA) ! velocity v        at 10m height
    real(RP), intent(out) :: T2     (IA,JA) ! temperature t     at  2m height
    real(RP), intent(out) :: Q2     (IA,JA) ! water vapor q     at  2m height

    real(RP) :: SFC_QSAT (IA,JA) ! saturatad water vapor mixing ratio [kg/kg]
    real(RP) :: LHV      (IA,JA)

    real(RP) :: Ustar   ! friction velocity [m]
    real(RP) :: Tstar   ! friction temperature [K]
    real(RP) :: Qstar   ! friction mixing rate [kg/kg]
    real(RP) :: Uabs    ! modified absolute velocity [m/s]
    real(RP) :: Ra      ! Aerodynamic resistance (=1/Ce) [1/s]
    real(RP) :: SFC_QV  ! water vapor mixing ratio [kg/kg]

    real(RP) :: FracU10 ! calculation parameter for U10 [-]
    real(RP) :: FracT2  ! calculation parameter for T2 [-]
    real(RP) :: FracQ2  ! calculation parameter for Q2 [-]

    integer  :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Surface flux(bulk)'

    ! ToDo consider ATM_TEMP is appropriate
    call HYDROMETEOR_LHV( &
         IA, IS, IE, JA, JS, JE, &
         ATM_TEMP(:,:), & ! [IN]
         LHV(:,:)       ) ! [OUT]

    call SATURATION_pres2qsat_all( &
         IA, IS, IE, JA, JS, JE, &
         SFC_TEMP(:,:), & ! [IN]
         SFC_PRES(:,:), & ! [IN]
         SFC_QSAT(:,:)  ) ! [OUT]

    !omp parallel do default(none) &
    !omp private(SFC_QV,Ustar,Tstar,Qstar,Uabs,Ra,FracU10,FracT2,FracQ2) &
    !omp shared (IS,IE,JS,JE, ATMOS_PHY_SF_beta, &
    !omp         ATM_TEMP,ATM_PRES,ATM_QV,ATM_U,ATM_V,ATM_Z1, &
    !omp         SFC_TEMP,SFC_PRES,SFC_Z0M,SFC_Z0H,SFC_Z0E,PBL, &
    !omp         SFLX_MW,SFLX_MU,SFLX_MW,SFLX_QV,U10,V10,T2,Q2)
    do j = JS, JE
    do i = IS, IE
       SFC_QV = ( 1.0_RP - ATMOS_PHY_SF_beta ) * ATM_QV(i,j) + ATMOS_PHY_SF_beta * SFC_QSAT(i,j)

       call BULKFLUX( Ustar,         & ! [OUT]
                      Tstar,         & ! [OUT]
                      Qstar,         & ! [OUT]
                      Uabs,          & ! [OUT]
                      Ra,            & ! [OUT]
                      FracU10,       & ! [OUT]
                      FracT2,        & ! [OUT]
                      FracQ2,        & ! [OUT]
                      ATM_TEMP(i,j), & ! [IN]
                      SFC_TEMP(i,j), & ! [IN]
                      ATM_PRES(i,j), & ! [IN]
                      SFC_PRES(i,j), & ! [IN]
                      ATM_QV  (i,j), & ! [IN]
                      SFC_QV       , & ! [IN]
                      ATM_U   (i,j), & ! [IN]
                      ATM_V   (i,j), & ! [IN]
                      ATM_Z1  (i,j), & ! [IN]
                      PBL     (i,j), & ! [IN]
                      SFC_Z0M (i,j), & ! [IN]
                      SFC_Z0H (i,j), & ! [IN]
                      SFC_Z0E (i,j)  ) ! [IN]

       !-----< momentum >-----
       SFLX_MW(i,j) = -SFC_DENS(i,j) * Ustar * Ustar / Uabs * ATM_W(i,j)
       SFLX_MU(i,j) = -SFC_DENS(i,j) * Ustar * Ustar / Uabs * ATM_U(i,j)
       SFLX_MV(i,j) = -SFC_DENS(i,j) * Ustar * Ustar / Uabs * ATM_V(i,j)

       !-----< heat flux >-----
       SFLX_SH(i,j) = -SFC_DENS(i,j) * Ustar * Tstar * CPdry
       SFLX_LH(i,j) = -SFC_DENS(i,j) * Ustar * Qstar * LHV(i,j)

       !-----< mass flux >-----
       SFLX_QV(i,j) = SFLX_LH(i,j) / LHV(i,j)

       !-----< U10, T2, q2 >-----
       !U10(i,j) = FracU10 * ATM_U(i,j)
       !V10(i,j) = FracU10 * ATM_V(i,j)
       !T2 (i,j) = ( 1.0_RP - FracT2 ) * SFC_TEMP(i,j) + FracT2 * ATM_TEMP(i,j)
       !Q2 (i,j) = ( 1.0_RP - FracQ2 ) * SFC_QV        + FracQ2 * ATM_QV  (i,j)

       U10(i,j) = ATM_U(i,j) * log( 10.0_RP / SFC_Z0M(i,j) ) / log( ATM_Z1(i,j) / SFC_Z0M(i,j) )
       V10(i,j) = ATM_V(i,j) * log( 10.0_RP / SFC_Z0M(i,j) ) / log( ATM_Z1(i,j) / SFC_Z0M(i,j) )
       T2 (i,j) = SFC_TEMP(i,j) + ( ATM_TEMP(i,j) - SFC_TEMP(i,j) ) &
                                * ( log(      2.0_RP / SFC_Z0M(i,j) ) * log(      2.0_RP / SFC_Z0H(i,j) ) ) &
                                / ( log( ATM_Z1(i,j) / SFC_Z0M(i,j) ) * log( ATM_Z1(i,j) / SFC_Z0H(i,j) ) )
       Q2 (i,j) = SFC_QV        + ( ATM_QV  (i,j) - SFC_QV        ) &
                                * ( log(      2.0_RP / SFC_Z0M(i,j) ) * log(      2.0_RP / SFC_Z0E(i,j) ) ) &
                                / ( log( ATM_Z1(i,j) / SFC_Z0M(i,j) ) * log( ATM_Z1(i,j) / SFC_Z0E(i,j) ) )
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_bulk_flux

end module scale_atmos_phy_sf_bulk
