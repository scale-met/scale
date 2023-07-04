!-------------------------------------------------------------------------------
!> module ocean / physics / surface albedo / nakajima00
!!
!! @par Description
!!          Ocean surface albedo nakajima00 model
!!
!! @par Reference
!!        Nakajima, T., M. Tsukamoto, Y. Tsushima, A. Numaguti, and T. Kimura, 2000: Modeling of the radiative process in an atmospheric general circulation model, Applied Optics, 39, 4869-4878, doi:10.1364/AO.39.004869
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_albedo_nakajima00
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
  public :: OCEAN_PHY_ALBEDO_nakajima00_setup
  public :: OCEAN_PHY_ALBEDO_nakajima00

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
  !-----------------------------------------------------------------------------
contains

  subroutine OCEAN_PHY_ALBEDO_nakajima00_setup

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ALBEDO_nakajima00_setup",*) 'Setup'

    return
  end subroutine OCEAN_PHY_ALBEDO_nakajima00_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ALBEDO_nakajima00( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       cosSZA,        &
!      tau,           &
       SFC_albedo     )
    use scale_cpl_sfc_index
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(in)  :: cosSZA    (OIA,OJA)           ! cos(solar zenith angle)         (0-1)
!   real(RP), intent(in)  :: tau       (OIA,OJA)           ! cloud optical thickness
    real(RP), intent(out) :: SFC_albedo(OIA,OJA,N_RAD_DIR) ! sea surface albedo (short wave) (0-1)

    real(RP), parameter :: c_ocean_albedo(5,3) = reshape( &
         (/ -2.8108_RP   , -1.3651_RP,  2.9210E1_RP, -4.3907E1_RP,  1.8125E1_RP, &
             6.5626E-1_RP, -8.7253_RP, -2.7749E1_RP,  4.9486E1_RP, -1.8345E1_RP, &
            -6.5423E-1_RP,  9.9967_RP,  2.7769_RP  , -1.7620E1_RP,  7.0838_RP    /), &
         (/ 5,3 /) )

    real(RP) :: am1
!   real(RP) :: tr1
    real(RP) :: s

    integer :: i, j, n
    !---------------------------------------------------------------------------
    !$acc data copyin(cosSZA) copyout(SFC_albedo)

    !$acc kernels
    !$omp parallel do private(am1,s)
    do j = OJS, OJE
    do i = OIS, OIE
       am1 = max( min( cosSZA(i,j), 0.961_RP ), 0.0349_RP )
!      tr1 = max( min( cosSZA(i,j) / ( 4.0_RP * tau(i,j) ), 1.0_RP ), 0.05_RP )

       s = 0.0_RP
!OCL UNROLL
       !$acc loop seq
       do n = 1, 5
          ! power of am1 is different form the paper, in which the number is wrong.
!         s = s + c_ocean_albedo(n,1) * tr1**(n-1) &
!               + c_ocean_albedo(n,2) * tr1**(n-1) * am1 &
!               + c_ocean_albedo(n,3) * tr1**(n-1) * am1**2
          s = s + c_ocean_albedo(n,1) &
                + c_ocean_albedo(n,2) * am1 &
                + c_ocean_albedo(n,3) * am1**2
       enddo

       SFC_albedo(i,j,I_R_direct ) = exp(s)
       SFC_albedo(i,j,I_R_diffuse) = SFC_albedo(i,j,I_R_direct)
    enddo
    enddo
    !$acc end kernels 

    !$acc end data
    return
  end subroutine OCEAN_PHY_ALBEDO_nakajima00

end module scale_ocean_phy_albedo_nakajima00
