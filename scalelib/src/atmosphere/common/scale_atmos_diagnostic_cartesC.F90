!-------------------------------------------------------------------------------
!> module atmosphere / diagnostic / CartesianC
!!
!! @par Description
!!          Calculate diagnostic variables for Cartesian-C grid
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_diagnostic_cartesC
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
  public :: ATMOS_DIAGNOSTIC_CARTESC_get_vel

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
  !-----------------------------------------------------------------------------
  !> ATMOS_DIAGNOSTIC_CARTESC_get_vel
  !! W, U, V
  !<
  subroutine ATMOS_DIAGNOSTIC_CARTESC_get_vel( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       W,    &
       U,    &
       V     )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_grid_cartesC_metric, only: &
       GSQRT => ATMOS_GRID_CARTESC_METRIC_GSQRT, &
       J13G  => ATMOS_GRID_CARTESC_METRIC_J13G,  &
       J23G  => ATMOS_GRID_CARTESC_METRIC_J23G
    use scale_atmos_grid_cartesC_index, only: &
       I_XYZ, &
       I_XYW
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)

    real(RP), intent(out) :: W   (KA,IA,JA)
    real(RP), intent(out) :: U   (KA,IA,JA)
    real(RP), intent(out) :: V   (KA,IA,JA)

    real(RP) :: momws

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! Note: W(KS,:,:) is filled because the values at i=1 or j=1 are not calculated below.
!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       W(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,momws)
    do j = max(2,JS), JE
    do i = max(2,IS), IE
       ! at KS+1/2
       momws = MOMZ(KS,i,j) &
             + ( J13G(KS,i,j,I_XYW) * ( MOMX(KS,i,j) + MOMX(KS,i-1,j) + MOMX(KS+1,i,j) + MOMX(KS+1,i-1,j) ) &
             + J23G(KS,i,j,I_XYW) * ( MOMY(KS,i,j) + MOMY(KS,i,j-1) + MOMY(KS+1,i,j) + MOMY(KS+1,i,j-1) ) ) &
             * 0.25_RP / GSQRT(KS,i,j,I_XYW)
       ! at KS
       ! momws at the surface is assumed to be zero
       W(KS,i,j) = ( momws * 0.5_RP                                                                &
                    - ( J13G(KS,i,j,I_XYZ) * ( MOMX(KS,i,j) + MOMX(KS,i-1,j) )                     &
                      + J23G(KS,i,j,I_XYZ) * ( MOMY(KS,i,j) + MOMY(KS,i,j-1) ) )                   &
                      * 0.5_RP / GSQRT(KS,i,j,I_XYZ)                              ) / DENS(KS,i,j)
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       W(KE,i,j) = 0.5_RP * ( MOMZ(KE-1,i,j) ) / DENS(KE,i,j)
    enddo
    enddo
 
!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = max(2,IS), IE
    do k = KS, KE
       U(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do k = KS, KE
       U(k,1,j) = MOMX(k,1,j) / DENS(k,1,j)
    enddo
    enddo

   !OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = max(2,JS), JE
    do i = IS, IE
    do k = KS, KE
       V(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
    do i = IS, IE
    do k = KS, KE
       V(k,i,1) = MOMY(k,i,1) / DENS(k,i,1)
    enddo
    enddo

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = JS, JE
    do i  = IS, IE
       W(   1:KS-1,i,j) = W(KS,i,j)
       U(   1:KS-1,i,j) = U(KS,i,j)
       V(   1:KS-1,i,j) = V(KS,i,j)
       W(KE+1:KA,  i,j) = W(KE,i,j)
       U(KE+1:KA,  i,j) = U(KE,i,j)
       V(KE+1:KA,  i,j) = V(KE,i,j)
    enddo
    enddo

    call COMM_vars8( W(:,:,:), 1 )
    call COMM_vars8( U(:,:,:), 2 )
    call COMM_vars8( V(:,:,:), 3 )
    call COMM_wait ( W(:,:,:), 1, .false. )
    call COMM_wait ( U(:,:,:), 2, .false. )
    call COMM_wait ( V(:,:,:), 3, .false. )

    return
  end subroutine ATMOS_DIAGNOSTIC_CARTESC_get_vel

end module scale_atmos_diagnostic_cartesC
