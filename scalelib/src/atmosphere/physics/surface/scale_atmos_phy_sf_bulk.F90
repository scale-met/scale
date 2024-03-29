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
#include "scalelib.h"
module scale_atmos_phy_sf_bulk
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
  real(RP), private :: ATMOS_PHY_SF_BULK_beta = 1.0_RP ! evaporation efficiency (0-1)
  !$acc declare create(ATMOS_PHY_SF_BULK_beta)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_bulk_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_ATMOS_PHY_SF_BULK / &
       ATMOS_PHY_SF_BULK_beta

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_bulk_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_SF_bulk_setup",*) 'Bulk scheme'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_BULK,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_SF_bulk_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_SF_bulk_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SF_BULK. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SF_BULK)

    !$acc update device(ATMOS_PHY_SF_BULK_beta)

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
       Ustar, Tstar, Qstar, Wstar,   &
       RLmo,                         &
       U10, V10, T2, Q2              )
    use scale_const, only: &
       EPS    => CONST_EPS, &
       CPdry  => CONST_CPdry, &
       EPSvap => CONST_EPSvap
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_atmos_saturation, only: &
       SATURATION_psat_all => ATMOS_SATURATION_psat_all
    use scale_bulkflux, only: &
       BULKFLUX, &
       BULKFLUX_diagnose_surface
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
    real(RP), intent(out) :: Ustar  (IA,JA) ! friction velocity
    real(RP), intent(out) :: Tstar  (IA,JA) ! temperatuer scale
    real(RP), intent(out) :: Qstar  (IA,JA) ! moisture scale
    real(RP), intent(out) :: Wstar  (IA,JA) ! convective veolocity scale
    real(RP), intent(out) :: RLmo   (IA,JA) ! inverse of Obukhov length
    real(RP), intent(out) :: U10    (IA,JA) ! velocity u        at 10m height
    real(RP), intent(out) :: V10    (IA,JA) ! velocity v        at 10m height
    real(RP), intent(out) :: T2     (IA,JA) ! temperature t     at  2m height
    real(RP), intent(out) :: Q2     (IA,JA) ! water vapor q     at  2m height

    real(RP) :: SFC_PSAT (IA,JA) ! saturatad water vapor pressure [Pa]
    real(RP) :: LHV      (IA,JA)

    real(RP) :: Uabs          ! modified absolute velocity [m/s]
    real(RP) :: Ra            ! Aerodynamic resistance (=1/Ce) [1/s]
    real(RP) :: SFC_QSAT      ! saturatad water vapor mixing ratio [kg/kg]
    real(RP) :: SFC_QV(IA,JA) ! water vapor mixing ratio [kg/kg]

    real(RP) :: FracU10(IA,JA) ! calculation parameter for U10 [-]
    real(RP) :: FracT2 (IA,JA) ! calculation parameter for T2 [-]
    real(RP) :: FracQ2 (IA,JA) ! calculation parameter for Q2 [-]

    real(RP) :: MFLUX

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / surface flux / bulk'

    !$acc data copyin(ATM_W,ATM_U,ATM_V,ATM_TEMP,ATM_PRES,ATM_QV,SFC_DENS,SFC_TEMP,SFC_PRES,SFC_Z0M,SFC_Z0H,SFC_Z0E,PBL,ATM_Z1) &
    !$acc      copyout(SFLX_MW,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_LH,SFLX_QV,Ustar,Tstar,Qstar,Wstar,RLmo,U10,V10,T2,Q2) &
    !$acc      create(SFC_PSAT,LHV,SFC_QV,FracU10,FracT2,FracQ2)


    call HYDROMETEOR_LHV( IA, IS, IE, JA, JS, JE, &
                          SFC_TEMP(:,:), & ! [IN]
                          LHV(:,:)       ) ! [OUT]

    call SATURATION_psat_all( IA, IS, IE, JA, JS, JE, &
                              SFC_TEMP(:,:), & ! [IN]
                              SFC_PSAT(:,:)  ) ! [OUT]

    !$omp parallel do &
#ifndef __GFORTRAN__
    !$omp default(none) &
    !$omp shared (IS,IE,JS,JE,EPSvap,ATMOS_PHY_SF_BULK_beta,EPS,CPdry,LHV,bulkflux,&
    !$omp         ATM_TEMP,ATM_PRES,ATM_QV,ATM_W,ATM_U,ATM_V,ATM_Z1,PBL, &
    !$omp         SFC_DENS,SFC_TEMP,SFC_PRES,SFC_QV,SFC_PSAT,SFC_Z0M,SFC_Z0H,SFC_Z0E, &
    !$omp         SFLX_MW,SFLX_MU,SFLX_MV,SFLX_SH,SFLX_LH,SFLX_QV, &
    !$omp         FracU10,FracT2,FracQ2, &
    !$omp         Ustar,Tstar,Qstar,Wstar,RLmo,U10,V10,T2,Q2) &
#else
    !$omp default(shared) &
#endif
    !$omp private(SFC_QSAT,Uabs,Ra,MFLUX)
    !$acc kernels
    !$acc loop independent collapse(2) &
    !$acc private(SFC_QSAT,Uabs,Ra,MFLUX)
    do j = JS, JE
    do i = IS, IE
       ! qdry = 1 - psat
       SFC_QSAT = EPSvap * SFC_PSAT(i,j) / ( SFC_PRES(i,j) - ( 1.0_RP-EPSvap ) * SFC_PSAT(i,j) )

       SFC_QV(i,j) = ( 1.0_RP - ATMOS_PHY_SF_BULK_beta ) * ATM_QV(i,j) + ATMOS_PHY_SF_BULK_beta * SFC_QSAT
       Uabs = sqrt( ATM_W(i,j)**2 + ATM_U(i,j)**2 + ATM_V(i,j)**2 )

       call BULKFLUX( ATM_TEMP(i,j), SFC_TEMP(i,j),             & ! [IN]
                      ATM_PRES(i,j), SFC_PRES(i,j),             & ! [IN]
                      ATM_QV  (i,j), SFC_QV  (i,j),             & ! [IN]
                      Uabs, ATM_Z1(i,j), PBL(i,j),              & ! [IN]
                      SFC_Z0M(i,j), SFC_Z0H(i,j), SFC_Z0E(i,j), & ! [IN]
                      Ustar(i,j), Tstar(i,j), Qstar(i,j),       & ! [OUT]
                      Wstar(i,j), RLmo(i,j), Ra,                & ! [OUT]
                      FracU10(i,j), FracT2(i,j), FracQ2(i,j)    ) ! [OUT]

       !-----< momentum >-----
       if ( Uabs < EPS ) then
          SFLX_MW(i,j) = 0.0_RP
          SFLX_MU(i,j) = 0.0_RP
          SFLX_MV(i,j) = 0.0_RP
       else
          MFLUX = -SFC_DENS(i,j) * Ustar(i,j)**2
          SFLX_MW(i,j) = MFLUX * ATM_W(i,j) / Uabs
          SFLX_MU(i,j) = MFLUX * ATM_U(i,j) / Uabs
          SFLX_MV(i,j) = MFLUX * ATM_V(i,j) / Uabs
       end if

       !-----< heat flux >-----
       SFLX_SH(i,j) = -SFC_DENS(i,j) * Ustar(i,j) * Tstar(i,j) * CPdry
       SFLX_LH(i,j) = -SFC_DENS(i,j) * Ustar(i,j) * Qstar(i,j) * LHV(i,j)

       !-----< mass flux >-----
       SFLX_QV(i,j) = SFLX_LH(i,j) / LHV(i,j)

    enddo
    enddo
    !$acc end kernels

    call BULKFLUX_diagnose_surface( IA, IS, IE, JA, JS, JE, &
                                    ATM_U(:,:), ATM_V(:,:),                    &
                                    ATM_TEMP(:,:), ATM_QV(:,:),                &
                                    SFC_TEMP(:,:), SFC_QV(:,:),                &
                                    ATM_Z1(:,:),                               &
                                    SFC_Z0M(:,:), SFC_Z0H(:,:), SFC_Z0E(:,:),  &
                                    U10(:,:), V10(:,:), T2(:,:), Q2(:,:),      &
                                    FracU10 = FracU10(:,:),                    &
                                    FracT2 = FracT2(:,:), FracQ2 = FracQ2(:,:) )

    !$acc end data

    return
  end subroutine ATMOS_PHY_SF_bulk_flux

end module scale_atmos_phy_sf_bulk
