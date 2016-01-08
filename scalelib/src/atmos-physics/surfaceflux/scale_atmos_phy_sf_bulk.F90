!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Bulk Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-03 (Y.Miyamoto)  [new]
!! @li      2011-12-11 (H.Yashiro)   [mod] integrate to SCALE-LES ver.3
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-04-10 (Y.Miyamoto)  [mod] introduce coefficients for interpolation
!! @li      2012-09-11 (S.Nishizawa) [mod] bugfix based on the scale document
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
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_bulk_setup
  public :: ATMOS_PHY_SF_bulk

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
  real(RP), private, parameter :: ATMOS_PHY_SF_U_maxM = 100.0_RP ! maximum limit of absolute velocity for momentum [m/s]
  real(RP), private, parameter :: ATMOS_PHY_SF_U_maxH = 100.0_RP ! maximum limit of absolute velocity for heat     [m/s]
  real(RP), private, parameter :: ATMOS_PHY_SF_U_maxE = 100.0_RP ! maximum limit of absolute velocity for vapor    [m/s]
  real(RP), private            :: ATMOS_PHY_SF_U_minM =   0.0_RP ! minimum limit of absolute velocity for momentum [m/s]
  real(RP), private            :: ATMOS_PHY_SF_U_minH =   0.0_RP ! minimum limit of absolute velocity for heat     [m/s]
  real(RP), private            :: ATMOS_PHY_SF_U_minE =   0.0_RP ! minimum limit of absolute velocity for vapor    [m/s]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_SF_bulk_setup( ATMOS_PHY_SF_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_phy_sf_bulkcoef, only: &
       SF_bulkcoef_setup => ATMOS_PHY_SF_bulkcoef_setup
    implicit none

    character(len=*), intent(in) :: ATMOS_PHY_SF_TYPE

    NAMELIST / PARAM_ATMOS_PHY_SF_BULK / &
       ATMOS_PHY_SF_U_minM, &
       ATMOS_PHY_SF_U_minH, &
       ATMOS_PHY_SF_U_minE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SURFACE FLUX] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Bulk scheme'

    if ( ATMOS_PHY_SF_TYPE /= 'BULK' ) then
       write(*,*) 'xxx ATMOS_PHY_SF_TYPE is not BULK. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_BULK,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_BULK. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_BULK)

    call SF_bulkcoef_setup

    return
  end subroutine ATMOS_PHY_SF_bulk_setup

  !-----------------------------------------------------------------------------
  !> Calculate surface flux
  subroutine ATMOS_PHY_SF_bulk( &
       ATM_TEMP, ATM_PRES, ATM_W, ATM_U, ATM_V,     &
       ATM_DENS,                                    &
       ATM_QTRC,                                    &
       ATM_Z1, dt,                                  &
       SFC_DENS, SFC_PRES,                          &
       SFLX_LW_dn, SFLX_SW_dn,                      &
       SFC_TEMP, SFC_albedo, SFC_beta,              &
       SFC_Z0M, SFC_Z0H, SFC_Z0E,                   &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH, &
       SFLX_QTRC,                                   &
       U10, V10, T2, Q2                             )
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       CPdry => CONST_CPdry, &
       Rdry  => CONST_Rdry
    use scale_atmos_phy_sf_bulkcoef, only: &
       SF_bulkcoef => ATMOS_PHY_SF_bulkcoef
    use scale_atmos_saturation, only: &
       SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_templhv
    use scale_roughness, only: &
       ROUGHNESS
    implicit none

    real(RP), intent(in)    :: ATM_TEMP  (IA,JA)    ! temperature at the lowermost layer (cell center) [K]
    real(RP), intent(in)    :: ATM_PRES  (IA,JA)    ! pressure    at the lowermost layer (cell center) [Pa]
    real(RP), intent(in)    :: ATM_W     (IA,JA)    ! velocity w  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in)    :: ATM_U     (IA,JA)    ! velocity u  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in)    :: ATM_V     (IA,JA)    ! velocity v  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in)    :: ATM_DENS  (IA,JA)    ! density     at the lowermost layer (cell center) [kg/m3]
    real(RP), intent(in)    :: ATM_QTRC  (IA,JA,QA) ! tracer      at the lowermost layer (cell center) [kg/kg]
    real(RP), intent(in)    :: ATM_Z1    (IA,JA)    ! height of the lowermost grid from surface (cell center) [m]
    real(DP), intent(in)    :: dt                   ! delta time
    real(RP), intent(in)    :: SFC_DENS  (IA,JA)    ! density     at the surface atmosphere [kg/m3]
    real(RP), intent(in)    :: SFC_PRES  (IA,JA)    ! pressure    at the surface atmosphere [Pa]
    real(RP), intent(in)    :: SFLX_LW_dn(IA,JA)    ! downward longwave  radiation flux at the surface [J/m2/s]
    real(RP), intent(in)    :: SFLX_SW_dn(IA,JA)    ! downward shortwave radiation flux at the surface [J/m2/s]
    real(RP), intent(in)    :: SFC_TEMP  (IA,JA)    ! temperature at the surface skin [K]
    real(RP), intent(in)    :: SFC_albedo(IA,JA,2)  ! surface albedo (LW/SW) [0-1]
    real(RP), intent(in)    :: SFC_beta  (IA,JA)    ! evaporation efficiency [0-1]
    real(RP), intent(inout) :: SFC_Z0M   (IA,JA)    ! surface roughness length (momentum) [m]
    real(RP), intent(inout) :: SFC_Z0H   (IA,JA)    ! surface roughness length (heat) [m]
    real(RP), intent(inout) :: SFC_Z0E   (IA,JA)    ! surface roughness length (vapor) [m]
    real(RP), intent(out)   :: SFLX_MW   (IA,JA)    ! surface flux for z-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out)   :: SFLX_MU   (IA,JA)    ! surface flux for x-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out)   :: SFLX_MV   (IA,JA)    ! surface flux for y-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out)   :: SFLX_SH   (IA,JA)    ! surface flux for sensible heat (area center)   [J/m2/s]
    real(RP), intent(out)   :: SFLX_LH   (IA,JA)    ! surface flux for latent   heat (area center)   [J/m2/s]
    real(RP), intent(out)   :: SFLX_QTRC (IA,JA,QA) ! surface flux for tracer mass   (area center)   [kg/m2/s]
    real(RP), intent(out)   :: U10       (IA,JA)    ! velocity u        at 10m height
    real(RP), intent(out)   :: V10       (IA,JA)    ! velocity v        at 10m height
    real(RP), intent(out)   :: T2        (IA,JA)    ! temperature t     at  2m height
    real(RP), intent(out)   :: Q2        (IA,JA)    ! water vapor q     at  2m height

    real(RP) :: SFC_Z0M_t(IA,JA)
    real(RP) :: SFC_Z0H_t(IA,JA)
    real(RP) :: SFC_Z0E_t(IA,JA)

    real(RP) :: ATM_Uabs(IA,JA) ! absolute velocity at z1 [m/s]
    real(RP) :: ATM_POTT(IA,JA) ! potential temperature at z1, based on the local surface pressure [K]
    real(RP) :: SFC_QSAT(IA,JA) ! saturatad water vapor mixing ratio [kg/kg]

    real(RP) :: Cm  (IA,JA)     ! bulk coefficient for momentum
    real(RP) :: Ch  (IA,JA)     ! bulk coefficient for heat
    real(RP) :: Ce  (IA,JA)     ! bulk coefficient for vapor
    real(RP) :: R10M(IA,JA)
    real(RP) :: R2H (IA,JA)
    real(RP) :: R2E (IA,JA)

    real(RP) :: Uabs_lim, RovCP
    real(RP) :: LHV (IA,JA)
    integer  :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface flux(bulk)'

    RovCP = Rdry / CPdry

    do j = JS, JE
    do i = IS, IE
       ATM_Uabs(i,j) = sqrt( ATM_W(i,j)*ATM_W(i,j) &
                           + ATM_U(i,j)*ATM_U(i,j) &
                           + ATM_V(i,j)*ATM_V(i,j) ) ! at cell center
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       ! potential temperature based on the local surface pressure (not 1000hPa)
       ATM_POTT(i,j) = ATM_TEMP(i,j) * ( SFC_PRES(i,j) / ATM_PRES(i,j) )**RovCP
    enddo
    enddo

    call ROUGHNESS( SFC_Z0M_t(:,:), & ! [OUT]
                    SFC_Z0H_t(:,:), & ! [OUT]
                    SFC_Z0E_t(:,:), & ! [OUT]
                    SFC_Z0M  (:,:), & ! [IN]
                    SFC_Z0H  (:,:), & ! [IN]
                    SFC_Z0E  (:,:), & ! [IN]
                    ATM_U    (:,:), & ! [IN]
                    ATM_V    (:,:), & ! [IN]
                    ATM_Z1   (:,:), & ! [IN]
                    dt              ) ! [IN]

    SFC_Z0M(:,:) = SFC_Z0M(:,:) + SFC_Z0M_t(:,:) * dt
    SFC_Z0H(:,:) = SFC_Z0H(:,:) + SFC_Z0H_t(:,:) * dt
    SFC_Z0E(:,:) = SFC_Z0E(:,:) + SFC_Z0E_t(:,:) * dt

    call SF_bulkcoef( ATM_Uabs(:,:), & ! [IN]
                      ATM_POTT(:,:), & ! [IN]
                      ATM_Z1  (:,:), & ! [IN]
                      SFC_TEMP(:,:), & ! [IN]
                      SFC_Z0M (:,:), & ! [IN]
                      SFC_Z0H (:,:), & ! [IN]
                      SFC_Z0E (:,:), & ! [IN]
                      Cm      (:,:), & ! [OUT]
                      Ch      (:,:), & ! [OUT]
                      Ce      (:,:), & ! [OUT]
                      R10M    (:,:), & ! [OUT]
                      R2H     (:,:), & ! [OUT]
                      R2E     (:,:)  ) ! [OUT]

    !-----< momentum >-----

    do j = JS, JE
    do i = IS, IE
       Uabs_lim = min( max( ATM_Uabs(i,j), ATMOS_PHY_SF_U_minM ), ATMOS_PHY_SF_U_maxM )

       SFLX_MW(i,j) = -Cm(i,j) * Uabs_lim * SFC_DENS(i,j) * ATM_W(i,j)
       SFLX_MU(i,j) = -Cm(i,j) * Uabs_lim * SFC_DENS(i,j) * ATM_U(i,j)
       SFLX_MV(i,j) = -Cm(i,j) * Uabs_lim * SFC_DENS(i,j) * ATM_V(i,j)
    enddo
    enddo

    !-----< heat flux >-----

    do j = JS, JE
    do i = IS, IE
       Uabs_lim = min( max( ATM_Uabs(i,j), ATMOS_PHY_SF_U_minH ), ATMOS_PHY_SF_U_maxH )

       SFLX_SH(i,j) = Ch(i,j) * Uabs_lim * SFC_DENS(i,j) * CPdry * ( SFC_TEMP(i,j) - ATM_POTT(i,j) )
    enddo
    enddo

    call SATURATION_pres2qsat_all( SFC_QSAT(:,:), & ! [OUT]
                                   SFC_TEMP(:,:), & ! [IN]
                                   SFC_PRES(:,:)  ) ! [IN]

    call ATMOS_THERMODYN_templhv( LHV, ATM_TEMP )

    do j = JS, JE
    do i = IS, IE
       Uabs_lim = min( max( ATM_Uabs(i,j), ATMOS_PHY_SF_U_minE ), ATMOS_PHY_SF_U_maxE )

       SFLX_LH(i,j) = Ce(i,j) * Uabs_lim * SFC_DENS(i,j) * LHV(i,j) * SFC_beta(i,j) * ( SFC_QSAT(i,j) - ATM_QTRC(i,j,I_QV) )
    enddo
    enddo

    !-----< mass flux >-----

    SFLX_QTRC(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       SFLX_QTRC(i,j,I_QV) = SFLX_LH(i,j) / LHV(i,j)
    enddo
    enddo

    !-----< U10, T2, q2 >-----

    do j = JS, JE
    do i = IS, IE
       U10   (i,j) = R10M(i,j) * ATM_U   (i,j)
       V10   (i,j) = R10M(i,j) * ATM_V   (i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       T2(i,j) = (          R2H(i,j) ) * ATM_TEMP(i,j) &
               + ( 1.0_RP - R2H(i,j) ) * SFC_TEMP(i,j)
       Q2(i,j) = (          R2E(i,j) ) * ATM_QTRC(i,j,I_QV) &
               + ( 1.0_RP - R2E(i,j) ) * SFC_beta(i,j) * SFC_QSAT(i,j)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_bulk

end module scale_atmos_phy_sf_bulk
