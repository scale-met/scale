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
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_kall,    &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_imin,    &
       ADM_imax
    implicit none

    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(RP), intent(in)  :: v     (ADM_gall   ,ADM_kall       ,ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_gall_in,1:(kmax-kmin+1),ADM_lall)

    integer :: i, j, k, l, n
    !---------------------------------------------------------------------------

    do l = 1,    ADM_lall
    do k = kmin, kmax
       n = 1
       do j = ADM_jmin, ADM_jmax+1
       do i = ADM_imin, ADM_imax+1
          v_clip(n,k-kmin+1,l) = v(suf(i,j),k,l)

          n = n + 1
       enddo
       enddo
    enddo
    enddo

    return
  end subroutine GTL_clip_region

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer( v, v_clip )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_imin,    &
       ADM_imax
    implicit none

    real(RP), intent(in)  :: v     (ADM_gall   ,ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_gall_in,ADM_lall)

    integer :: i, j, l, n
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       n = 1
       do j = ADM_jmin, ADM_jmax+1
       do i = ADM_imin, ADM_imax+1
          v_clip(n,l) = v(suf(i,j),l)

          n = n + 1
       enddo
       enddo
    enddo

    return
  end subroutine GTL_clip_region_1layer

  !-----------------------------------------------------------------------------
  subroutine GTL_clip_region_1layer_k(v,v_clip,ksize,k)
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall,    &
       ADM_gall_in, &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_imin,    &
       ADM_imax
    implicit none

    integer,  intent(in)  :: ksize
    real(RP), intent(in)  :: v     (ADM_gall   ,ksize, ADM_lall)
    real(RP), intent(out) :: v_clip(ADM_gall_in,ADM_lall)
    integer,  intent(in)  :: k

    integer :: i, j, l, n
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       n = 1
       do j = ADM_jmin, ADM_jmax+1
       do i = ADM_imin, ADM_imax+1
          v_clip(n,l) = v(suf(i,j),k,l)

          n = n + 1
       enddo
       enddo
    enddo

    return
  end subroutine GTL_clip_region_1layer_k

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_gtl
