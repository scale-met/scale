!-------------------------------------------------------------------------------
!> Module generic tool
!!
!! @par Description
!!         This module is for the generic subroutine, e.g., global mean.
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_gtl
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
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
  public :: GTL_clip_region
  public :: GTL_clip_region_1layer
  public :: GTL_clip_region_1layer_k

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region( v, v_clip, kmin, kmax )
    use mod_adm, only: &
       ADM_gall,        &
       ADM_kall,        &
       ADM_lall,        &
       ADM_GIoJo,       &
       ADM_IopJop_nmax, &
       ADM_IopJop
    implicit none

    integer, intent(in)  :: kmin
    integer, intent(in)  :: kmax
    real(RP), intent(in)  :: v     (ADM_gall,       ADM_kall,       ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_IopJop_nmax,1:(kmax-kmin+1),ADM_lall)

    integer :: n, k, l
    !---------------------------------------------------------------------------

    do l = 1,    ADM_lall
    do k = kmin, kmax
    do n = 1,    ADM_IopJop_nmax

       v_clip(n,k-kmin+1,l) = v(ADM_IopJop(n,ADM_GIoJo),k,l)

    enddo
    enddo
    enddo

    return
  end subroutine GTL_clip_region

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer( v, v_clip )
    use mod_adm, only: &
       ADM_gall,        &
       ADM_lall,        &
       ADM_GIoJo,       &
       ADM_IopJop_nmax, &
       ADM_IopJop
    implicit none

    real(RP), intent(in)  :: v     (ADM_gall,       ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_IopJop_nmax,ADM_lall)

    integer :: n, l
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do n = 1, ADM_IopJop_nmax

       v_clip(n,l) = v(ADM_IopJop(n,ADM_GIoJo),l)

    enddo
    enddo

    return
  end subroutine GTL_clip_region_1layer

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer_k(v,v_clip,ksize,k)
    use mod_adm, only: &
       ADM_gall,        &
       ADM_lall,        &
       ADM_GIoJo,       &
       ADM_IopJop_nmax, &
       ADM_IopJop
    implicit none

    integer, intent(in)  :: ksize
    real(RP), intent(in)  :: v     (ADM_gall,ksize, ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_IopJop_nmax,ADM_lall)
    integer, intent(in)  :: k

    integer :: n, l
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do n = 1, ADM_IopJop_nmax

       v_clip(n,l) = v(ADM_IopJop(n,ADM_GIoJo),k,l)

    enddo
    enddo

    return
  end subroutine GTL_clip_region_1layer_k

end module mod_gtl
