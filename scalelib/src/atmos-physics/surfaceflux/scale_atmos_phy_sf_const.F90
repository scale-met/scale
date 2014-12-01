!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Constant flux
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-03 (Y.Miyamoto)  [new]
!! @li      2011-12-11 (H.Yashiro)   [mod] integrate to SCALE-LES ver.3
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-04-10 (Y.Miyamoto)  [mod] introduce coefficients for interpolation
!! @li      2012-09-11 (S.Nishizawa) [mod] bugfix based on the scale document
!! @li      2012-09-12 (Y.Sato)    [renew] constant FLUX version
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_sf_const
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
  public :: ATMOS_PHY_SF_const_setup
  public :: ATMOS_PHY_SF_const

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
  subroutine ATMOS_PHY_SF_const_setup( ATMOS_PHY_SF_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: ATMOS_PHY_SF_TYPE

    NAMELIST / PARAM_ATMOS_PHY_SF_CONST / &
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

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SURFACE FLUX] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Constant flux'

    if ( ATMOS_PHY_SF_TYPE /= 'CONST' ) then
       write(*,*) 'xxx ATMOS_PHY_SF_TYPE is not CONST. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_CONST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_CONST. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_CONST)

    return
  end subroutine ATMOS_PHY_SF_const_setup

  !-----------------------------------------------------------------------------
  !> Constant flux
  subroutine ATMOS_PHY_SF_const( &
       ATM_TEMP, ATM_PRES, ATM_W, ATM_U, ATM_V,     &
       ATM_DENS,                                    &
       ATM_QTRC,                                    &
       ATM_Z1,                                      &
       SFC_DENS, SFC_PRES,                          &
       SFLX_LW_dn, SFLX_SW_dn,                      &
       SFC_TEMP, SFC_albedo, SFC_beta,              &
       SFC_Z0,                                      &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH, &
       SFLX_QTRC,                                   &
       Uabs10, U10, V10, T2, Q2                     )
    use scale_const, only: &
       PI   => CONST_PI,   &
       LHV  => CONST_LHV,  &
       STB  => CONST_STB
    use scale_time, only: &
       TIME_NOWSEC
    implicit none

    real(RP), intent(in)    :: ATM_TEMP  (IA,JA)    ! temperature at the lowermost layer (cell center) [K]
    real(RP), intent(in)    :: ATM_PRES  (IA,JA)    ! pressure    at the lowermost layer (cell center) [Pa]
    real(RP), intent(in)    :: ATM_W     (IA,JA)    ! velocity w  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in)    :: ATM_U     (IA,JA)    ! velocity u  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in)    :: ATM_V     (IA,JA)    ! velocity v  at the lowermost layer (cell center) [m/s]
    real(RP), intent(in)    :: ATM_DENS  (IA,JA)    ! density     at the lowermost layer (cell center) [kg/m3]
    real(RP), intent(in)    :: ATM_QTRC  (IA,JA,QA) ! tracer      at the lowermost layer (cell center) [kg/kg]
    real(RP), intent(in)    :: ATM_Z1    (IA,JA)    ! height of the lowermost grid from surface (cell center) [m]
    real(RP), intent(in)    :: SFC_DENS  (IA,JA)    ! density     at the surface atmosphere [kg/m3]
    real(RP), intent(in)    :: SFC_PRES  (IA,JA)    ! pressure    at the surface atmosphere [Pa]
    real(RP), intent(in)    :: SFLX_LW_dn(IA,JA)    ! downward longwave  radiation flux at the surface [J/m2/s]
    real(RP), intent(in)    :: SFLX_SW_dn(IA,JA)    ! downward shortwave radiation flux at the surface [J/m2/s]
    real(RP), intent(in)    :: SFC_TEMP  (IA,JA)    ! temperature at the surface skin [K]
    real(RP), intent(in)    :: SFC_albedo(IA,JA,2)  ! surface albedo (LW/SW) [0-1]
    real(RP), intent(in)    :: SFC_beta  (IA,JA)    ! evaporation efficiency [0-1]
    real(RP), intent(inout) :: SFC_Z0    (IA,JA)    ! surface roughness length (momentum) [m]
    real(RP), intent(out)   :: SFLX_MW   (IA,JA)    ! surface flux for z-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out)   :: SFLX_MU   (IA,JA)    ! surface flux for x-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out)   :: SFLX_MV   (IA,JA)    ! surface flux for y-momentum    (area center)   [m/s*kg/m2/s]
    real(RP), intent(out)   :: SFLX_SH   (IA,JA)    ! surface flux for sensible heat (area center)   [J/m2/s]
    real(RP), intent(out)   :: SFLX_LH   (IA,JA)    ! surface flux for latent   heat (area center)   [J/m2/s]
    real(RP), intent(out)   :: SFLX_QTRC (IA,JA,QA) ! surface flux for tracer mass   (area center)   [kg/m2/s]
    real(RP), intent(out)   :: Uabs10    (IA,JA)    ! absolute velocity at 10m height
    real(RP), intent(out)   :: U10       (IA,JA)    ! velocity u        at 10m height
    real(RP), intent(out)   :: V10       (IA,JA)    ! velocity v        at 10m height
    real(RP), intent(out)   :: T2        (IA,JA)    ! temperature t     at  2m height
    real(RP), intent(out)   :: Q2        (IA,JA)    ! water vapor q     at  2m height

    real(RP) :: ATM_Uabs(IA,JA) ! absolute velocity at z1 [m/s]
    real(RP) :: SFC_QSAT(IA,JA) ! saturatad water vapor mixing ratio [kg/kg]

    real(RP) :: Cm(IA,JA)       ! bulk coefficient for momentum
    real(RP) :: R10, R2

    real(RP) :: modulation
    real(RP) :: Uabs_lim
    integer  :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface flux(const)'

    do j = JS, JE
    do i = IS, IE
       ATM_Uabs(i,j) = sqrt( ATM_W(i,j)*ATM_W(i,j) &
                           + ATM_U(i,j)*ATM_U(i,j) &
                           + ATM_V(i,j)*ATM_V(i,j) ) ! at cell center
    enddo
    enddo

    if   ( ATMOS_PHY_SF_FLG_MOM_FLUX == 0 ) then ! Bulk coefficient is constant
       do j = JS, JE
       do i = IS, IE
          Cm(i,j) = ATMOS_PHY_SF_Const_Cm
       enddo
       enddo
    elseif( ATMOS_PHY_SF_FLG_MOM_FLUX == 1 ) then ! Friction velocity is constant
       do j = JS, JE
       do i = IS, IE
          Cm(i,j) = ( ATMOS_PHY_SF_Const_Ustar / ATM_Uabs(i,j) )**2
          Cm(i,j) = min( max( Cm(i,j), ATMOS_PHY_SF_Cm_min ), ATMOS_PHY_SF_Cm_max )
       enddo
       enddo
    endif

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

    if ( ATMOS_PHY_SF_FLG_SH_DIURNAL ) then
       modulation = sin( 2.0_RP * PI * TIME_NOWSEC / 3600.0_RP / ATMOS_PHY_SF_Const_FREQ )
    else
       modulation = 1.0_RP
    endif

    do j = JS, JE
    do i = IS, IE
       SFLX_SH(i,j) = ATMOS_PHY_SF_Const_SH * modulation
       SFLX_LH(i,j) = ATMOS_PHY_SF_Const_LH
    enddo
    enddo

    !-----< mass flux >-----

    SFLX_QTRC(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       SFLX_QTRC(i,j,I_QV) = SFLX_LH(i,j) / LHV
    enddo
    enddo

    !-----< U10, T2, q2 >-----

    do j = JS, JE
    do i = IS, IE
       R10 = 10.0_RP / ATM_Z1(i,j)

       Uabs10(i,j) = R10 * ATM_Uabs(i,j)
       U10   (i,j) = R10 * ATM_U   (i,j)
       V10   (i,j) = R10 * ATM_V   (i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       T2(i,j) = ATM_TEMP(i,j)
       Q2(i,j) = ATM_QTRC(i,j,I_QV)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_const

end module scale_atmos_phy_sf_const
