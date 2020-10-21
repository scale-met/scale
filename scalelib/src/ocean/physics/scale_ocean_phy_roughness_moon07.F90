!-------------------------------------------------------------------------------
!> module ocean / physics / surface roughness length / moon07
!!
!! @par Description
!!          surface roughness length Moon07 scheme
!!      IL-JU MOON, ISAAC GINIS, TETSU HARA, AND BIJU THOMAS, 2007:
!!      A Physics-Based Parameterization of Air-Sea Momentum Flux at High Wind Speeds
!!      and Its Impact on Hurricane Intensity Predictions, Mon. Wea. Rev., 135, 2869-2878
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_roughness_moon07
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
  public :: OCEAN_PHY_ROUGHNESS_moon07_setup
  public :: OCEAN_PHY_ROUGHNESS_moon07

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
  integer,  private :: OCEAN_PHY_ROUGHNESS_moon07_itelim  = 10        ! maximum iteration number

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_moon07_setup
    use scale_prc, only: &
       PRC_abort
    use scale_ocean_phy_roughness, only: &
       OCEAN_PHY_ROUGHNESS_setup
    implicit none

    namelist / PARAM_OCEAN_PHY_ROUGHNESS_MOON07 / &
       OCEAN_PHY_ROUGHNESS_moon07_itelim

    integer :: ierr
    !---------------------------------------------------------------------------

    ! common setup
    call OCEAN_PHY_ROUGHNESS_setup

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ROUGHNESS_MOON07,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ROUGHNESS_moon07_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ROUGHNESS_moon07_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ROUGHNESS_MOON07. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ROUGHNESS_MOON07)

    return
  end subroutine OCEAN_PHY_ROUGHNESS_moon07_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_moon07( &
       OIA, OIS, OIE, OJA, OJS, OJE, &
       Uabs, Z1,     &
       mask,         &
       Z0M, Z0H, Z0E )
    use scale_const, only: &
       UNDEF  => CONST_UNDEF,  &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN
    use scale_ocean_phy_roughness, only: &
       OCEAN_PHY_ROUGHNESS_visck, &
       OCEAN_PHY_ROUGHNESS_Ustar_min, &
       OCEAN_PHY_ROUGHNESS_Z0M_min, &
       OCEAN_PHY_ROUGHNESS_Z0H_min, &
       OCEAN_PHY_ROUGHNESS_Z0E_min
    implicit none
    integer, intent(in) :: OIA, OIS, OIE
    integer, intent(in) :: OJA, OJS, OJE

    real(RP), intent(in) :: Uabs(OIA,OJA) ! velocity at the lowest atomspheric layer [m/s]
    real(RP), intent(in) :: Z1  (OIA,OJA) ! cell center height at the lowest atmospheric layer [m]
    logical,  intent(in) :: mask(OIA,OJA)

    real(RP), intent(inout) :: Z0M(OIA,OJA) ! roughness length for momentum [m]
    real(RP), intent(out)   :: Z0H(OIA,OJA) ! roughness length for heat [m]
    real(RP), intent(out)   :: Z0E(OIA,OJA) ! roughness length for vapor [m]

    ! works
    real(RP) :: Ustar
    real(RP) :: U10M

    integer  :: ite
    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ &
    !$omp shared(OJS,OJE,OIS,OIE, &
    !$omp        GRAV,UNDEF, &
    !$omp        OCEAN_PHY_ROUGHNESS_Z0M_min,OCEAN_PHY_ROUGHNESS_Z0H_min,OCEAN_PHY_ROUGHNESS_Z0E_min, &
    !$omp        OCEAN_PHY_ROUGHNESS_Ustar_min,OCEAN_PHY_ROUGHNESS_visck,OCEAN_PHY_ROUGHNESS_moon07_itelim, &
    !$omp        Z0M,Z0H,Z0E,Uabs,Z1,mask) &
    !$omp private(i,j,ite,U10M,Ustar)
    do j = OJS, OJE
    do i = OIS, OIE
       if ( mask(i,j) ) then

          Z0M(i,j) = max( Z0M(i,j), OCEAN_PHY_ROUGHNESS_Z0M_min )

          do ite = 1, OCEAN_PHY_ROUGHNESS_moon07_itelim
             Ustar = max( KARMAN * Uabs(i,j) / log( Z1(i,j)/Z0M(i,j) ), OCEAN_PHY_ROUGHNESS_Ustar_min )
             U10M = Ustar / KARMAN * log( 10.0_RP/Z0M(i,j) )

             if ( U10M <= 12.5_RP ) then
                Z0M(i,j) = 0.0185_RP * Ustar**2 / GRAV
             else
                Z0M(i,j) = ( 0.085_RP * ( -  0.56_RP  * Ustar**2   &
                                          + 20.255_RP * Ustar      &
                                          +  2.458_RP            ) &
                           - 0.58_RP ) *  1.0E-3_RP
             end if
             Z0M(i,j) = max( Z0M(i,j), OCEAN_PHY_ROUGHNESS_Z0M_min )
          enddo

          !  Fairall et al. TOGA V3.0
          !  Fairall et al. (2003) JCLI, vol. 16, 571-591. Eq. (28)
          Z0H(i,j) = min( 5.5E-5_RP / ( Z0M(i,j) * Ustar / OCEAN_PHY_ROUGHNESS_visck )**0.6_RP, &
                          1.1E-4_RP )
          Z0E(i,j) = Z0H(i,j)
          Z0H(i,j) = max( Z0H(i,j), OCEAN_PHY_ROUGHNESS_Z0H_min )
          Z0E(i,j) = max( Z0E(i,j), OCEAN_PHY_ROUGHNESS_Z0E_min )

       else
          Z0H(i,j) = UNDEF
          Z0E(i,j) = UNDEF
       end if
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ROUGHNESS_moon07

end module scale_ocean_phy_roughness_moon07
