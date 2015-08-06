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
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GTL_global_sum
  public :: GTL_global_sum_srf
  public :: GTL_global_sum_eachlayer
  public :: GTL_global_mean
  public :: GTL_max
  public :: GTL_max_k
  public :: GTL_min
  public :: GTL_min_k

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
  function GTL_global_sum( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_kall,        &
       ADM_kmin,        &
       ADM_kmax,        &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_comm, only: &
       COMM_Stat_sum
    use mod_vmtr, only: &
       VMTR_VOLUME,   &
       VMTR_VOLUME_pl
    implicit none

    real(RP), intent(in) :: var   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP)             :: sum_g

    real(RP) :: sum
    integer :: n, k, l
    !---------------------------------------------------------------------------

    sum = 0.0_RP
    do l = 1,        ADM_lall
    do k = ADM_kmin, ADM_kmax
    do n = 1,        ADM_IooJoo_nmax
       sum = sum + var(ADM_IooJoo(n,ADM_GIoJo),k,l) * VMTR_VOLUME(ADM_IooJoo(n,ADM_GIoJo),k,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,        ADM_lall_pl
       do k = ADM_kmin, ADM_kmax
          sum = sum + var_pl(ADM_GSLF_PL,k,l) * VMTR_VOLUME_pl(ADM_GSLF_PL,k,l)
       enddo
       enddo
    endif

    call COMM_Stat_sum( sum, sum_g )

    return
  end function GTL_global_sum

  !-----------------------------------------------------------------------------
  function GTL_global_sum_srf( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_KNONE,       &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_comm, only: &
       COMM_Stat_sum
    use mod_gmtr, only: &
       GMTR_area,   &
       GMTR_area_pl
    implicit none

    real(RP), intent(in) :: var   (ADM_gall,   ADM_KNONE,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP)             :: sum_g

    real(RP) :: sum
    integer :: n, l
    !---------------------------------------------------------------------------

    sum = 0.0_RP
    do l = 1, ADM_lall
    do n = 1, ADM_IooJoo_nmax
       sum = sum + var      (ADM_IooJoo(n,ADM_GIoJo),ADM_KNONE,l) &
                 * GMTR_area(ADM_IooJoo(n,ADM_GIoJo),l)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          sum = sum + var_pl      (ADM_GSLF_PL,ADM_KNONE,l) &
                    * GMTR_area_pl(ADM_GSLF_PL,l)
       enddo
    endif

    call COMM_Stat_sum( sum, sum_g )

    return
  end function GTL_global_sum_srf

  !-----------------------------------------------------------------------------
  subroutine GTL_global_sum_eachlayer( var, var_pl, sum_g )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_kall,        &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_comm, only: &
       COMM_Stat_sum_eachlayer
    use mod_gmtr, only: &
       GMTR_area,   &
       GMTR_area_pl
    use mod_vmtr, only: &
       VMTR_GAM2,   &
       VMTR_GAM2_pl
    implicit none

    real(RP), intent(in)  :: var   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: sum_g (ADM_kall)

    real(RP) :: sum(ADM_kall)
    integer :: n, k, l
    !---------------------------------------------------------------------------

    sum(:) = 0.0_RP
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_IooJoo_nmax
       sum(k) = sum(k) + var      (ADM_IooJoo(n,ADM_GIoJo),k,l) &
                       * GMTR_area(ADM_IooJoo(n,ADM_GIoJo),l)   &
                       * VMTR_GAM2(ADM_IooJoo(n,ADM_GIoJo),k,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          sum(k) = sum(k) + var_pl      (ADM_GSLF_PL,k,l) &
                          * GMTR_area_pl(ADM_GSLF_PL,l)   &
                          * VMTR_GAM2_pl(ADM_GSLF_PL,k,l)
       enddo
       enddo
    endif

    call COMM_Stat_sum_eachlayer( ADM_kall, sum(:), sum_g(:) )

    return
  end subroutine GTL_global_sum_eachlayer

  !-----------------------------------------------------------------------------
  function GTL_global_mean( var, var_pl ) result( sum_g )
    use mod_adm, only: &
       ADM_gall,    &
       ADM_lall,    &
       ADM_gall_pl, &
       ADM_lall_pl, &
       ADM_kall
    use mod_comm, only: &
       COMM_Stat_sum
    implicit none

    real(RP), intent(in) :: var   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP)             :: sum_g

    real(RP)       :: one   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP)       :: one_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    logical, save :: first = .true.
    real(RP), save :: volume_g
    !---------------------------------------------------------------------------

    if ( first ) then
       !--- calc global volume at first time
       one   (:,:,:) = 1.0_RP
       one_pl(:,:,:) = 1.0_RP

       volume_g = GTL_global_sum( one(:,:,:), one_pl(:,:,:) )

       first = .false.
    endif

    sum_g = GTL_global_sum( var(:,:,:), var_pl(:,:,:) )

    sum_g = sum_g / volume_g

    return
  end function GTL_global_mean

  !-----------------------------------------------------------------------------
  function GTL_max( var, var_pl, kdim, kstart, kend ) result( vmax_g )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_comm, only: &
       COMM_Stat_max
    implicit none

    integer, intent(in) :: kdim
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    real(RP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(RP)             :: vmax_g

    real(RP) :: vmax
    integer :: n, k, l
    !---------------------------------------------------------------------------

    vmax = -1.E+30_RP
    do l = 1,      ADM_lall
    do k = kstart, kend
    do n = 1,      ADM_IooJoo_nmax
       vmax = max( vmax, var(ADM_IooJoo(n,ADM_GIoJo),k,l) )
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,      ADM_lall_pl
       do k = kstart, kend
          vmax = max( vmax, var_pl(ADM_GSLF_PL,k,l) )
       enddo
       enddo
    endif

    call COMM_Stat_max( vmax, vmax_g )

    return
  end function GTL_max

  !-----------------------------------------------------------------------------
  function GTL_max_k( var, var_pl, k ) result( vmax_g )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_kall,        &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_comm, only: &
       COMM_Stat_max
    implicit none

    real(RP), intent(in) :: var   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in) :: k
    real(RP)             :: vmax_g

    real(RP) :: vmax
    integer :: n, l
    !---------------------------------------------------------------------------

    vmax = -1.E+30_RP
    do l = 1,        ADM_lall
    do n = 1,        ADM_IooJoo_nmax
       vmax = max( vmax, var(ADM_IooJoo(n,ADM_GIoJo),k,l) )
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,        ADM_lall_pl
          vmax = max( vmax, var_pl(ADM_GSLF_PL,k,l) )
       enddo
    endif

    call COMM_Stat_max( vmax, vmax_g )

    return
  end function GTL_max_k

  !-----------------------------------------------------------------------------
  function GTL_min( var, var_pl, kdim, kstart, kend, nonzero ) result( vmin_g )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_comm, only: &
       COMM_Stat_min
    implicit none

    integer, intent(in) :: kdim
    integer, intent(in) :: kstart
    integer, intent(in) :: kend
    real(RP), intent(in) :: var   (ADM_gall,   kdim,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,kdim,ADM_lall_pl)
    real(RP)             :: vmin_g
    logical, optional, intent(in) :: nonzero

    real(RP) :: vmin
    integer :: n, k, l
    !---------------------------------------------------------------------------

    if ( present(nonzero) ) then
       if ( nonzero ) then

          vmin = 1.E+30_RP
          do l = 1,      ADM_lall
          do k = kstart, kend
          do n = 1,      ADM_IooJoo_nmax
             if (       var(ADM_IooJoo(n,ADM_GIoJo),k,l) > 0.0_RP &
                  .AND. var(ADM_IooJoo(n,ADM_GIoJo),k,l) < vmin ) then

                vmin = var(ADM_IooJoo(n,ADM_GIoJo),k,l)

             endif
          enddo
          enddo
          enddo

       if ( ADM_have_pl ) then
          do l = 1,      ADM_lall_pl
          do k = kstart, kend
             if (       var_pl(ADM_GSLF_PL,k,l) > 0.0_RP &
                  .AND. var_pl(ADM_GSLF_PL,k,l) < vmin ) then

                vmin = var_pl(ADM_GSLF_PL,k,l)

             endif
          enddo
          enddo
       endif

       call COMM_Stat_min( vmin, vmin_g )
       return

       endif
    endif

    vmin = 1.E+30_RP
    do l = 1,      ADM_lall
    do k = kstart, kend
    do n = 1,      ADM_IooJoo_nmax
       vmin = min( vmin, var(ADM_IooJoo(n,ADM_GIoJo),k,l) )
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,      ADM_lall_pl
       do k = kstart, kend
          vmin = min( vmin, var_pl(ADM_GSLF_PL,k,l) )
       enddo
       enddo
    endif

    call COMM_Stat_min( vmin, vmin_g )

    return
  end function GTL_min

  !-----------------------------------------------------------------------------
  function GTL_min_k( var, var_pl, k ) result( vmin_g )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_kall,        &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_comm, only: &
       COMM_Stat_min
    implicit none

    real(RP), intent(in) :: var   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in) :: k
    real(RP)             :: vmin_g

    real(RP) :: vmin
    integer :: n, l
    !---------------------------------------------------------------------------

    vmin = +1.E+30_RP
    do l = 1,        ADM_lall
    do n = 1,        ADM_IooJoo_nmax
       vmin = min( vmin, var(ADM_IooJoo(n,ADM_GIoJo),k,l) )
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1,        ADM_lall_pl
          vmin = min( vmin, var_pl(ADM_GSLF_PL,k,l) )
       enddo
    endif

    call COMM_Stat_min( vmin, vmin_g )

    return
  end function GTL_min_k

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
