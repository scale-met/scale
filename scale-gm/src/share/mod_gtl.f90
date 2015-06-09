!-------------------------------------------------------------------------------
!>
!! Generic tool module
!!
!! @par Description
!!         This module is for the generic subroutine, e.g., global mean.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2004-06-23 (        )  [add] GTL_input_var2]
!! @li      2011-03-02 (        )  [add] GTL_clip_region_1layer_k
!! @li      2011-07-22 (T.Ohno)    MPI_Bcasts of poles data are suppressed when poles are excluded from communication.
!! @li      2012-06-06 (M.Terai)   Modification to reduce communication by FJSE
!!
!<
!-------------------------------------------------------------------------------
module mod_gtl
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS,    &
     ADM_MAXFNAME
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

  public :: GTL_input_var2_da    ! direct access
  public :: GTL_output_var2_da   ! direct access

  public :: GTL_generate_vxvyvz
  public :: GTL_generate_uv
  public :: GTL_mk_rigidrotation

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
  subroutine GTL_input_var2_da( basename, var, k_start, k_end, recnum, input_size )
    use mod_misc, only: &
       MISC_get_available_fid, &
       MISC_make_idstr
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_gall,    &
       ADM_lall
    implicit none

    character(len=*), intent(in)  :: basename
    integer,          intent(in)  :: k_start
    integer,          intent(in)  :: k_end
    real(RP),          intent(out) :: var(ADM_gall,k_start:k_end,ADM_lall)
    integer,          intent(in)  :: recnum
    integer,          intent(in)  :: input_size

    real(SP) :: var4(ADM_gall,k_start:k_end)
    real(DP) :: var8(ADM_gall,k_start:k_end)

    character(len=ADM_MAXFNAME) :: fname

    integer :: fid
    integer :: l, rgnid
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       call MISC_make_idstr(fname,trim(basename),'rgn',rgnid)
       fid = MISC_get_available_fid()
       open( unit   = fid,           &
             file   = trim(fname),   &
             form   = 'unformatted', &
             access = 'direct',      &
             recl   = ADM_gall*(k_end-k_start+1)*input_size, &
             status = 'old'          )

       if ( input_size == 4 ) then
          read(fid,rec=recnum) var4(:,:)
          var(:,k_start:k_end,l) = real(var4(:,k_start:k_end),kind=RP)
       elseif( input_size == 8 ) then
          read(fid,rec=recnum) var8(:,:)
          var(:,k_start:k_end,l) = real(var8(:,k_start:k_end),kind=RP)
       endif

       close(fid)
    enddo

    return
  end subroutine GTL_input_var2_da

  !-----------------------------------------------------------------------------
  subroutine GTL_output_var2_da( basename, var, k_start, k_end, recnum, output_size )
    use mod_misc, only: &
       MISC_get_available_fid, &
       MISC_make_idstr
    use mod_adm, only: &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_gall,    &
       ADM_lall
    use mod_cnst, only: &
       CNST_UNDEF4
    implicit none

    character(len=ADM_MAXFNAME), intent(in)  :: basename
    integer,                     intent(in)  :: k_start
    integer,                     intent(in)  :: k_end
    real(RP),                     intent(in)  :: var(:,:,:)
    integer,                     intent(in)  :: recnum
    integer,                     intent(in)  :: output_size

    real(4) :: var4(ADM_gall,k_start:k_end)
    real(RP) :: var8(ADM_gall,k_start:k_end)

    character(len=ADM_MAXFNAME) :: fname

    integer :: fid
    integer :: l, rgnid
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)

       call MISC_make_idstr(fname,trim(basename),'rgn',rgnid)
       fid = MISC_get_available_fid()
       open( unit   = fid,           &
             file   = trim(fname),   &
             form   = 'unformatted', &
             access = 'direct',      &
             recl   = ADM_gall*(k_end-k_start+1)*output_size, &
             status = 'unknown'      )

       if ( output_size == 4 ) then
         var4(:,k_start:k_end) = real(var(:,k_start:k_end,l),kind=4)

         where( var4(:,k_start:k_end) < CNST_UNDEF4+1.0_RP )
            var4(:,k_start:k_end) = CNST_UNDEF4
         end where

         write(fid,rec=recnum) var4(:,:)
       elseif( output_size == 8 ) then
         var8(:,k_start:k_end) = var(:,k_start:k_end,l)
         write(fid,rec=recnum) var8(:,:)
       endif

       close(fid)
    enddo

    return
  end subroutine GTL_output_var2_da

  !-----------------------------------------------------------------------------
  subroutine GTL_generate_vxvyvz( &
       ucos, ucos_pl, &
       vcos, vcos_pl, &
       vx,   vx_pl,   &
       vy,   vy_pl,   &
       vz,   vz_pl    )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_KNONE
    use mod_gmtr, only: &
       P_IX  => GMTR_P_IX,  &
       P_IY  => GMTR_P_IY,  &
       P_IZ  => GMTR_P_IZ,  &
       P_JX  => GMTR_P_JX,  &
       P_JY  => GMTR_P_JY,  &
       P_JZ  => GMTR_P_JZ,  &
       P_LAT => GMTR_P_LAT, &
       GMTR_P_var,          &
       GMTR_P_var_pl
    implicit none

    real(RP), intent(in)  :: ucos   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: ucos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vcos   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vcos_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vx     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vy     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: vz     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: vz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: u, v, coslat, sw

    integer :: n, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_gall
       coslat = cos(GMTR_P_var(n,k0,l,P_LAT))

       sw = 0.5_RP + sign(0.5_RP,-abs(coslat)) ! if (coslat == 0), u=v=0

       u = ucos(n,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )
       v = vcos(n,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )

       vx(n,k,l) = u * GMTR_P_var(n,k0,l,P_IX) &
                 + v * GMTR_P_var(n,k0,l,P_JX)
       vy(n,k,l) = u * GMTR_P_var(n,k0,l,P_IY) &
                 + v * GMTR_P_var(n,k0,l,P_JY)
       vz(n,k,l) = u * GMTR_P_var(n,k0,l,P_IZ) &
                 + v * GMTR_P_var(n,k0,l,P_JZ)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do n = 1, ADM_gall_pl
          coslat = cos(GMTR_P_var_pl(n,k0,l,P_LAT))

          sw = 0.5_RP + sign(0.5_RP,-abs(coslat)) ! if (coslat == 0), u=v=0

          u = ucos_pl(n,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )
          v = vcos_pl(n,k,l) * ( 1.0_RP - sw ) / ( coslat - sw )

          vx_pl(n,k,l) = u * GMTR_P_var_pl(n,k0,l,P_IX) &
                       + v * GMTR_P_var_pl(n,k0,l,P_JX)
          vy_pl(n,k,l) = u * GMTR_P_var_pl(n,k0,l,P_IY) &
                       + v * GMTR_P_var_pl(n,k0,l,P_JY)
          vz_pl(n,k,l) = u * GMTR_P_var_pl(n,k0,l,P_IZ) &
                       + v * GMTR_P_var_pl(n,k0,l,P_JZ)
       enddo
       enddo
       enddo
    endif

    return
  end subroutine GTL_generate_vxvyvz

  !-----------------------------------------------------------------------------
  subroutine GTL_generate_uv( &
       u,  u_pl,  &
       v,  v_pl,  &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl, &
       icos       )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_KNONE
    use mod_gmtr, only: &
       P_IX  => GMTR_P_IX,  &
       P_IY  => GMTR_P_IY,  &
       P_IZ  => GMTR_P_IZ,  &
       P_JX  => GMTR_P_JX,  &
       P_JY  => GMTR_P_JY,  &
       P_JZ  => GMTR_P_JZ,  &
       P_LAT => GMTR_P_LAT, &
       GMTR_P_var,          &
       GMTR_P_var_pl
    implicit none

    real(RP), intent(out) :: u    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: u_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: v    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: v_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer, optional, intent(in) :: icos

    integer :: n, k, l, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    if (       present(icos) &
         .AND. icos /= 0     ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          u(n,k,l) = ( vx(n,k,l) * GMTR_P_var(n,k0,l,P_IX) &
                     + vy(n,k,l) * GMTR_P_var(n,k0,l,P_IY) &
                     + vz(n,k,l) * GMTR_P_var(n,k0,l,P_IZ) ) * cos(GMTR_P_var(n,k0,l,P_LAT))
          v(n,k,l) = ( vx(n,k,l) * GMTR_P_var(n,k0,l,P_JX) &
                     + vy(n,k,l) * GMTR_P_var(n,k0,l,P_JY) &
                     + vz(n,k,l) * GMTR_P_var(n,k0,l,P_JZ) ) * cos(GMTR_P_var(n,k0,l,P_LAT))
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             u_pl(n,k,l) = ( vx_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_IX) &
                           + vy_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_IY) &
                           + vz_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_IZ) ) * cos(GMTR_P_var_pl(n,k0,l,P_LAT))
             v_pl(n,k,l) = ( vx_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_JX) &
                           + vy_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_JY) &
                           + vz_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_JZ) ) * cos(GMTR_P_var_pl(n,k0,l,P_LAT))
          enddo
          enddo
          enddo
       endif

    else

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do n = 1, ADM_gall
          u(n,k,l) = ( vx(n,k,l) * GMTR_P_var(n,k0,l,P_IX) &
                     + vy(n,k,l) * GMTR_P_var(n,k0,l,P_IY) &
                     + vz(n,k,l) * GMTR_P_var(n,k0,l,P_IZ) )
          v(n,k,l) = ( vx(n,k,l) * GMTR_P_var(n,k0,l,P_JX) &
                     + vy(n,k,l) * GMTR_P_var(n,k0,l,P_JY) &
                     + vz(n,k,l) * GMTR_P_var(n,k0,l,P_JZ) )
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do n = 1, ADM_gall_pl
             u_pl(n,k,l) = ( vx_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_IX) &
                           + vy_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_IY) &
                           + vz_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_IZ) )
             v_pl(n,k,l) = ( vx_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_JX) &
                           + vy_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_JY) &
                           + vz_pl(n,k,l) * GMTR_P_var_pl(n,k0,l,P_JZ) )
          enddo
          enddo
          enddo
       endif

    endif

    return
  end subroutine GTL_generate_uv

  !-----------------------------------------------------------------------------
  subroutine GTL_mk_rigidrotation( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl, &
       alpha,     &
       vmax       )
    use mod_adm, only: &
       ADM_have_pl,     &
       ADM_gall,        &
       ADM_gall_pl,     &
       ADM_lall,        &
       ADM_lall_pl,     &
       ADM_kall,        &
       ADM_KNONE,       &
       ADM_IooJoo_nmax, &
       ADM_IooJoo,      &
       ADM_GIoJo,       &
       ADM_GSLF_PL
    use mod_gmtr, only: &
       P_IX  => GMTR_P_IX,  &
       P_IY  => GMTR_P_IY,  &
       P_IZ  => GMTR_P_IZ,  &
       P_JX  => GMTR_P_JX,  &
       P_JY  => GMTR_P_JY,  &
       P_JZ  => GMTR_P_JZ,  &
       P_LON => GMTR_P_LON, &
       P_LAT => GMTR_P_LAT, &
       GMTR_P_var,          &
       GMTR_P_var_pl
    implicit none

    real(RP), intent(inout) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: alpha
    real(RP), intent(in)    :: vmax

    real(RP) :: u, v
    integer :: n, k, l, ij, k0
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do n = 1, ADM_IooJoo_nmax
       ij = ADM_IooJoo(n,ADM_GIoJo)

       u =  vmax * ( cos(GMTR_P_var(ij,k0,l,P_LAT)) * cos(alpha) &
                   + sin(GMTR_P_var(ij,k0,l,P_LAT))              &
                   * cos(GMTR_P_var(ij,k0,l,P_LON)) * sin(alpha) )
       v = -vmax * ( sin(GMTR_P_var(ij,k0,l,P_LON)) * sin(alpha) )

       vx(ij,k,l) = u * GMTR_P_var(ij,k0,l,P_IX) &
                  + v * GMTR_P_var(ij,k0,l,P_JX)
       vy(ij,k,l) = u * GMTR_P_var(ij,k0,l,P_IY) &
                  + v * GMTR_P_var(ij,k0,l,P_JY)
       vz(ij,k,l) = u * GMTR_P_var(ij,k0,l,P_IZ) &
                  + v * GMTR_P_var(ij,k0,l,P_JZ)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       ij = ADM_GSLF_PL

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          u =  vmax * ( cos(GMTR_P_var_pl(ij,k0,l,P_LAT)) * cos(alpha) &
                      + sin(GMTR_P_var_pl(ij,k0,l,P_LAT))              &
                      * cos(GMTR_P_var_pl(ij,k0,l,P_LON)) * sin(alpha) )
          v = -vmax * ( sin(GMTR_P_var_pl(ij,k0,l,P_LON)) * sin(alpha) )

          vx_pl(ij,k,l) = u * GMTR_P_var_pl(ij,k0,l,P_IX) &
                        + v * GMTR_P_var_pl(ij,k0,l,P_JX)
          vy_pl(ij,k,l) = u * GMTR_P_var_pl(ij,k0,l,P_IY) &
                        + v * GMTR_P_var_pl(ij,k0,l,P_JY)
          vz_pl(ij,k,l) = u * GMTR_P_var_pl(ij,k0,l,P_IZ) &
                        + v * GMTR_P_var_pl(ij,k0,l,P_JZ)
       enddo
       enddo
    endif

    return
  end subroutine GTL_mk_rigidrotation

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
  ! 2011/03/02 NEC [Add]
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
!-------------------------------------------------------------------------------
