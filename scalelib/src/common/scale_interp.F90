!-------------------------------------------------------------------------------
!> module INTERPOLATION
!!
!! @par Description
!!          INTERPOLATION module for nesting system
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_interp
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
  public :: INTERP_setup
  public :: INTERP_domain_compatibility
  public :: INTERP_factor1d
  public :: INTERP_factor2d
  public :: INTERP_factor3d
  public :: INTERP_factor2d_linear_latlon
  public :: INTERP_factor2d_linear_xy
  public :: INTERP_factor2d_weight
  public :: INTERP_factor3d_linear_latlon
  public :: INTERP_factor3d_linear_xy
  public :: INTERP_factor3d_weight

  public :: INTERP_interp1d
  public :: INTERP_interp2d
  public :: INTERP_interp3d

  interface INTERP_factor2d
     procedure INTERP_factor2d_linear_latlon
     procedure INTERP_factor2d_linear_xy
     procedure INTERP_factor2d_weight
  end interface INTERP_factor2d
  interface INTERP_factor3d
     procedure INTERP_factor3d_linear_latlon
     procedure INTERP_factor3d_linear_xy
     procedure INTERP_factor3d_weight
  end interface INTERP_factor3d

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: INTERP_search_horiz
  private :: INTERP_insert

  private :: haversine

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: large_number = 9.999E+15_RP

  integer,  private :: INTERP_weight_order = 2
  real(RP), private :: INTERP_search_limit
  real(RP), private :: INTERP_buffer_size_fact = 2.0_RP
  logical,  private :: INTERP_use_spline_vert = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTERP_setup( &
       weight_order, &
       search_limit  )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer,  intent(in) :: weight_order
    real(RP), intent(in), optional :: search_limit

    namelist /PARAM_INTERP/ &
         INTERP_buffer_size_fact, &
         INTERP_use_spline_vert

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("INTERP_setup",*) 'Setup'

    INTERP_weight_order = weight_order

    INTERP_search_limit = large_number
    if ( present(search_limit) ) then
       INTERP_search_limit = search_limit
       LOG_INFO("INTERP_setup",*) 'search limit [m] : ', INTERP_search_limit
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_INTERP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("INTERP_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("INTERP_setup",*) 'Not appropriate names in namelist PARAM_INTERP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_INTERP)

    return
  end subroutine INTERP_setup

  !-----------------------------------------------------------------------------
  subroutine INTERP_domain_compatibility( &
       lon_org,  &
       lat_org,  &
       topc_org, &
       lon_loc,  &
       lat_loc,  &
       topc_loc, &
       topf_loc, &
       skip_x,   &
       skip_y,   &
       skip_z    )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in) :: lon_org (:,:)
    real(RP), intent(in) :: lat_org (:,:)
    real(RP), intent(in) :: topc_org(:,:) ! full level
    real(RP), intent(in) :: lon_loc (:,:)
    real(RP), intent(in) :: lat_loc (:,:)
    real(RP), intent(in) :: topc_loc(:,:) ! full level
    real(RP), intent(in) :: topf_loc(:,:) ! half level (ceil)

    logical,  intent(in), optional :: skip_z
    logical,  intent(in), optional :: skip_x
    logical,  intent(in), optional :: skip_y

    real(RP) :: max_ref, min_ref
    real(RP) :: max_loc, min_loc

    logical  :: skip_z_
    logical  :: skip_x_
    logical  :: skip_y_

    intrinsic size
    !---------------------------------------------------------------------------

    skip_z_ = .false.
    if ( present(skip_z) ) then
       skip_z_ = skip_z
    endif

    skip_x_ = .false.
    if ( present(skip_x) ) then
       skip_x_ = skip_x
    endif

    skip_y_ = .false.
    if ( present(skip_y) ) then
       skip_y_ = skip_y
    endif

    if ( .NOT. skip_z_ ) then
       max_ref = maxval( topc_org(:,:) ) + minval( topf_loc(:,:) - topc_loc(:,:) ) ! not real top boundary (only for check)
       max_loc = maxval( topc_loc(:,:) )

       if ( max_ref < max_loc ) then
          LOG_ERROR("INTERP_domain_compatibility",*) 'REQUESTED DOMAIN IS TOO MUCH BROAD'
          LOG_ERROR_CONT(*) '-- VERTICAL direction over the limit'
          LOG_ERROR_CONT(*) '-- reference max: ', max_ref
          LOG_ERROR_CONT(*) '--     local max: ', max_loc
          call PRC_abort
       endif
    endif

    if ( .NOT. skip_x_ ) then
       max_ref = maxval( lon_org(:,:) / D2R )
       min_ref = minval( lon_org(:,:) / D2R )
       max_loc = maxval( lon_loc(:,:) / D2R )
       min_loc = minval( lon_loc(:,:) / D2R )

       if    ( (min_ref+360.0_RP-max_ref) < 360.0_RP / size(lon_org,1) * 2.0_RP ) then
          ! cyclic OK
       elseif( max_ref < max_loc .OR. min_ref > min_loc ) then
          LOG_ERROR("INTERP_domain_compatibility",*) 'REQUESTED DOMAIN IS TOO MUCH BROAD'
          LOG_ERROR_CONT(*) '-- LONGITUDINAL direction over the limit'
          LOG_ERROR_CONT(*) '-- reference max: ', max_ref
          LOG_ERROR_CONT(*) '-- reference min: ', min_ref
          LOG_ERROR_CONT(*) '--     local max: ', max_loc
          LOG_ERROR_CONT(*) '--     local min: ', min_loc
          call PRC_abort
       endif
    endif

    if ( .NOT. skip_y_ ) then
       max_ref = maxval( lat_org(:,:) / D2R )
       min_ref = minval( lat_org(:,:) / D2R )
       max_loc = maxval( lat_loc(:,:) / D2R )
       min_loc = minval( lat_loc(:,:) / D2R )

       if ( max_ref < max_loc .OR. min_ref > min_loc ) then
          LOG_ERROR("INTERP_domain_compatibility",*) 'REQUESTED DOMAIN IS TOO MUCH BROAD'
          LOG_ERROR_CONT(*) '-- LATITUDINAL direction over the limit'
          LOG_ERROR_CONT(*) '-- reference max: ', max_ref
          LOG_ERROR_CONT(*) '-- reference min: ', min_ref
          LOG_ERROR_CONT(*) '--     local max: ', max_loc
          LOG_ERROR_CONT(*) '--     local min: ', min_loc
          call PRC_abort
       endif
    endif

    return
  end subroutine INTERP_domain_compatibility

  !-----------------------------------------------------------------------------
  ! vertical search of interpolation points for two-points
!OCL SERIAL
  subroutine INTERP_factor1d( &
       KA_ref, KS_ref, KE_ref, &
       KA, KS, KE,             &
       hgt_ref,                &
       hgt,                    &
       idx_k,                  &
       vfact,                  &
       flag_extrap             )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       EPS   => CONST_EPS,  &
       UNDEF => CONST_UNDEF
    implicit none
    integer,  intent(in)  :: KA_ref, KS_ref, KE_ref        ! number of z-direction    (reference)
    integer,  intent(in)  :: KA, KS, KE                    ! number of z-direction    (target)

    real(RP), intent(in)  :: hgt_ref(KA_ref)               ! height [m]               (reference)
    real(RP), intent(in)  :: hgt    (KA)                   ! height [m]               (target)

    integer,  intent(out) :: idx_k  (KA,2)                 ! k-index in reference     (target)
    real(RP), intent(out) :: vfact  (KA)                   ! horizontal interp factor (target)

    logical,  intent(in), optional :: flag_extrap          ! when true, extrapolation will be executed (just copy)

    logical :: flag_extrap_

    integer  :: k, kk
    !---------------------------------------------------------------------------

    if ( present(flag_extrap) ) then
       flag_extrap_ = flag_extrap
    else
       flag_extrap_ = .true.
    end if

    do k = KS, KE
       idx_k(k,1) = -1
       idx_k(k,2) = -1

       if    ( hgt(k) <  hgt_ref(KS_ref) - EPS ) then
          if ( flag_extrap_ ) then
             idx_k(k,1) = KS_ref
             idx_k(k,2) = -1
             vfact(k) = 1.0_RP
          else
             idx_k(k,1) = -1
             idx_k(k,2) = -1
             vfact(k) = UNDEF
          end if
       elseif( hgt(k) < hgt_ref(KS_ref) ) then
          idx_k(k,1) = KS_ref
          idx_k(k,2) = -1
          vfact(k) = 1.0_RP
       elseif( hgt(k) > hgt_ref(KE_ref) + EPS ) then
          if ( flag_extrap_ ) then
             idx_k(k,1) = KE_ref
             idx_k(k,2) = -1
             vfact(k) = 1.0_RP
          else
             idx_k(k,1) = -1
             idx_k(k,2) = -1
             vfact(k) = UNDEF
          end if
       elseif( hgt(k) >= hgt_ref(KE_ref) ) then
          idx_k(k,1) = KE_ref
          idx_k(k,2) = -1
          vfact(k) = 1.0_RP
       else
          do kk = KS_ref, KE_ref-1
             if (       hgt(k) >= hgt_ref(kk  ) &
                  .AND. hgt(k) <  hgt_ref(kk+1) ) then
                idx_k(k,1) = kk
                idx_k(k,2) = kk + 1
                vfact(k) = ( hgt_ref(kk+1) - hgt    (k)  ) &
                         / ( hgt_ref(kk+1) - hgt_ref(kk) )

                exit
             endif
          enddo
       endif

    enddo ! k-loop

    return
  end subroutine INTERP_factor1d


  !-----------------------------------------------------------------------------
  ! make interpolation factor using bi-linear method on the latlon coordinate
  ! This can be used only for structured grid
  subroutine INTERP_factor2d_linear_latlon( &
       IA_ref, JA_ref,     &
       IA, JA,             &
       lon_ref, lat_ref,   &
       lon, lat,           &
       idx_i, idx_j, hfact )
    use scale_sort, only: &
       SORT_exec
    implicit none
    integer,  intent(in)  :: IA_ref, JA_ref
    integer,  intent(in)  :: IA, JA

    real(RP), intent(in)  :: lon_ref(IA_ref)
    real(RP), intent(in)  :: lat_ref(JA_ref)
    real(RP), intent(in)  :: lon(IA,JA)
    real(RP), intent(in)  :: lat(IA,JA)
    integer,  intent(out) :: idx_i(IA,JA,4)
    integer,  intent(out) :: idx_j(IA,JA,4)
    real(RP), intent(out) :: hfact(IA,JA,4)

    real(RP) :: f1, f2
    integer  :: i, j, ii, jj

    call PROF_rapstart('INTERP_fact',3)

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(f1,f2)
    do j = 1, JA
    do i = 1, IA

       ! longitude
       if ( lon(i,j) < lon_ref(1) ) then
          idx_i(i,j,:) = 1
          hfact(i,j,:) = 0.0_RP
       else if ( lon(i,j) > lon_ref(IA_ref) ) then
          idx_i(i,j,:) = IA_ref
          hfact(i,j,:) = 0.0_RP
       else if ( lon(i,j) == lon_ref(IA_ref) ) then
          idx_i(i,j,:) = IA_ref
          hfact(i,j,1) = 0.0_RP
          hfact(i,j,2) = 1.0_RP
          hfact(i,j,3) = 0.0_RP
          hfact(i,j,4) = 1.0_RP
       else
          do ii = 1, IA_ref
             if ( lon(i,j) < lon_ref(ii) ) then
                f1 = ( lon_ref(ii) - lon(i,j)      ) / ( lon_ref(ii) - lon_ref(ii-1) )
                f2 = ( lon(i,j)    - lon_ref(ii-1) ) / ( lon_ref(ii) - lon_ref(ii-1) )
                hfact(i,j,1) = f1
                hfact(i,j,2) = f2
                hfact(i,j,3) = f1
                hfact(i,j,4) = f2
                idx_i(i,j,1) = ii-1
                idx_i(i,j,2) = ii
                idx_i(i,j,3) = ii-1
                idx_i(i,j,4) = ii
                exit
             end if
          end do
       end if

       ! latitude
       if ( lat(i,j) < lat_ref(1) ) then
          idx_j(i,j,:) = 1
          hfact(i,j,:) = 0.0_RP
       else if ( lat(i,j) > lat_ref(JA_ref) ) then
          idx_j(i,j,:) = JA_ref
          hfact(i,j,:) = 0.0_RP
       else if ( lat(i,j) == lat_ref(JA_ref) ) then
          idx_j(i,j,:) = JA_ref
          hfact(i,j,1) = 0.0_RP
          hfact(i,j,2) = 0.0_RP
          !hfact(i,j,3) = hfact(i,j,3)
          !hfact(i,j,4) = hfact(i,j,4)
       else
          do jj = 1, JA_ref
             if ( lat(i,j) < lat_ref(jj) ) then
                f1 = ( lat_ref(jj) - lat(i,j)      ) / ( lat_ref(jj) - lat_ref(jj-1) )
                f2 = ( lat(i,j)    - lat_ref(jj-1) ) / ( lat_ref(jj) - lat_ref(jj-1) )
                hfact(i,j,1) = hfact(i,j,1) * f1
                hfact(i,j,2) = hfact(i,j,2) * f1
                hfact(i,j,3) = hfact(i,j,3) * f2
                hfact(i,j,4) = hfact(i,j,4) * f2
                idx_j(i,j,1) = jj-1
                idx_j(i,j,2) = jj-1
                idx_j(i,j,3) = jj
                idx_j(i,j,4) = jj
                exit
             end if
          end do
       end if

       call SORT_exec( 4,                                        & ! [IN]
                       hfact(i,j,:), idx_i(i,j,:), idx_j(i,j,:), & ! [INOUT]
                       reverse = .true.                          ) ! [IN]

    end do
    end do

    call PROF_rapend  ('INTERP_fact',3)

    return
  end subroutine INTERP_factor2d_linear_latlon

  !-----------------------------------------------------------------------------
  ! make interpolation factor using bi-linear method on the xy coordinate
  ! This can be used only for structured grid
  subroutine INTERP_factor2d_linear_xy( &
       IA_ref, JA_ref,     &
       IA, JA,             &
       x_ref, y_ref,       &
       x, y,               &
       idx_i, idx_j,       &
       hfact,              &
       zonal, pole         )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_sort, only: &
       SORT_exec
    implicit none
    integer,  intent(in)  :: IA_ref, JA_ref
    integer,  intent(in)  :: IA, JA

    real(RP), intent(in)  :: x_ref(IA_ref,JA_ref)
    real(RP), intent(in)  :: y_ref(IA_ref,JA_ref)
    real(RP), intent(in)  :: x(IA)
    real(RP), intent(in)  :: y(JA)

    integer,  intent(out) :: idx_i(IA,JA,4)
    integer,  intent(out) :: idx_j(IA,JA,4)
    real(RP), intent(out) :: hfact(IA,JA,4)

    logical, intent(in), optional :: zonal
    logical, intent(in), optional :: pole

    real(RP) :: u, v
    integer  :: inc_i, inc_j
    logical  :: error, err
    logical  :: zonal_, pole_

    integer :: ii0, jj0
    integer :: ii, jj
    integer :: i1, i2, i3, i4
    integer :: j1, j2, j3, j4
    integer :: i, j
    integer :: ite, ite_max

    call PROF_rapstart('INTERP_fact',3)

    if ( present(zonal) ) then
       zonal_ = zonal
    else
       zonal_ = .false.
    end if
    if ( present(pole) ) then
       pole_ = pole
    else
       pole_ = .false.
    end if

    ii0 = IA_ref * 0.5_RP
    jj0 = JA_ref * 0.5_RP
    error = .false.

    ite_max = IA_ref * JA_ref / 4

    !$omp parallel do &
    !$omp private(inc_i,inc_j,ii,jj,i1,i2,i3,i4,j1,j2,j3,j4,u,v,err) &
    !$omp firstprivate(ii0,jj0)
    do j = 1, JA
    do i = 1, IA

       if ( i==1 ) then
          ii = ii0
          jj = jj0
       end if
       do ite = 1, ite_max
          i1 = ii
          i2 = ii + 1
          i3 = ii + 1
          i4 = ii
          j1 = jj
          j2 = jj
          j3 = jj + 1
          j4 = jj + 1
          if ( ii == IA_ref ) then ! zonal_ must be .true.
             i2 = 1
             i3 = 1
          end if
          if ( jj == 0 ) then ! pole_ must be .true.
             j1 = 1
             j2 = 1
             i2 = IA_ref / 2 + ii
             if ( i2 > IA_ref ) i2 = i2 - IA_ref
             i1 = i2 + 1
             if ( i1 > IA_ref ) i1 = i1 - IA_ref
          end if
          if ( jj == JA_ref ) then ! pole_ must be .true.
             j3 = JA_ref
             j4 = JA_ref
             i3 = IA_ref / 2 + ii
             if ( i3 > IA_ref ) i3 = i3 - IA_ref
             i4 = i3 + 1
             if ( i4 > IA_ref ) i4 = i4 - IA_ref
          end if
          if ( x_ref(i1,j1)==UNDEF .or. x_ref(i2,j2)==UNDEF .or. x_ref(i3,j3)==UNDEF .or. x_ref(i4,j4)==UNDEF ) then
             if ( ii == IA_ref-1 ) then
                ii = 1
                if ( jj == JA_ref-1 ) then
                   jj = 1
                else
                   jj = jj + 2
                end if
             else
                ii = ii + 2
             end if
             if ( ii == IA_ref ) then
                ii = 1
                jj = jj + 2
             end if
             if ( jj == JA_ref ) then
                jj = 1
                ii = ii + 2
                if ( ii == IA_ref ) ii = 1
             end if
          else
             call INTERP_check_inside( &
                  x_ref(i1,j1), x_ref(i2,j2), x_ref(i3,j3), x_ref(i4,j4), & ! (in)
                  y_ref(i1,j1), y_ref(i2,j2), y_ref(i3,j3), y_ref(i4,j4), & ! (in)
                  x(i), y(j),                                             & ! (in)
                  inc_i, inc_j,                                           & ! (out)
                  err                                                     ) ! (out)
             if ( err ) error = .true.
             if ( error ) exit
             if ( inc_i == 0 .and. inc_j == 0 ) then ! inside the quadrilateral
                call INTERP_bilinear_inv( &
                     x_ref(i1,j1), x_ref(i2,j2), x_ref(i3,j3), x_ref(i4,j4), & ! (in)
                     y_ref(i1,j1), y_ref(i2,j2), y_ref(i3,j3), y_ref(i4,j4), & ! (in)
                     x(i), y(j),                                             & ! (in)
                     u, v,                                                   & ! (out)
                     err                                                     ) ! (out)
                if ( err ) error = .true.
                if ( error ) exit
                idx_i(i,j,1) = i1
                idx_i(i,j,2) = i2
                idx_i(i,j,3) = i3
                idx_i(i,j,4) = i4
                idx_j(i,j,1) = j1
                idx_j(i,j,2) = j2
                idx_j(i,j,3) = j3
                idx_j(i,j,4) = j4
                hfact(i,j,1) = ( 1.0_RP - u ) * ( 1.0_RP - v )
                hfact(i,j,2) = ( u          ) * ( 1.0_RP - v )
                hfact(i,j,3) = ( u          ) * ( v          )
                hfact(i,j,4) = ( 1.0_RP - u ) * ( v          )
                exit
             end if
             ii = ii + inc_i
             jj = jj + inc_j
             if ( zonal_ .or. pole_ ) then
                if ( ii == 0        ) ii = IA_ref
                if ( ii == IA_ref+1 ) ii = 1
                if ( pole_ ) then
                   jj = max( jj, 0 )
                   jj = min( jj, JA_ref )
                else
                   jj = max( jj, 1 )
                   jj = min( jj, JA_ref-1 )
                end if
             else
                if ( ii == 0 .and. jj == 0 ) then
                   ii = IA_ref - 1
                   jj = JA_ref - 1
                end if
                if ( ii == 0 .and. jj == JA_ref ) then
                   ii = IA_ref - 1
                   jj = 1
                end if
                if ( ii == IA_ref .and. jj == 0 ) then
                   ii = 1
                   jj = JA_ref - 1
                end if
                if ( ii == IA_ref .and. jj == JA_ref  ) then
                   ii = 1
                   jj = 1
                end if
                ii = max( min( ii, IA_ref-1 ), 1 )
                jj = max( min( jj, JA_ref-1 ), 1 )
             end if
          end if
       end do
       if ( ite == ite_max+1 ) then
          LOG_ERROR("INTERP_factor2d_linear_xy",*) 'iteration max has been reached', i, j, x(i), y(j)
          LOG_ERROR_CONT(*) minval(x_ref), maxval(x_ref), minval(y_ref), maxval(y_ref)
          LOG_ERROR_CONT(*) x_ref(1,1), x_ref(IA_ref,1),x_ref(1,JA_ref)
          LOG_ERROR_CONT(*) y_ref(1,1), y_ref(IA_ref,1),y_ref(1,JA_ref)
          idx_i(i,j,:) = 1
          idx_j(i,j,:) = 1
          hfact(i,j,:) = 0.0_RP
          call PRC_abort
       end if

       if ( i==1 ) then
          ii0 = ii
          jj0 = jj
       end if

       if ( error ) exit

       call SORT_exec( 4,                                        & ! [IN]
                       hfact(i,j,:), idx_i(i,j,:), idx_j(i,j,:), & ! [INOUT]
                       reverse = .true.                          ) ! [IN]

    end do
    end do

    if ( error ) call PRC_abort

    call PROF_rapend  ('INTERP_fact',3)

    return
  end subroutine INTERP_factor2d_linear_xy

  !-----------------------------------------------------------------------------
  ! make interpolation factor using Lat-Lon with weight based on the distance
  subroutine INTERP_factor2d_weight( &
       npoints,             &
       IA_ref, JA_ref,      &
       IA, JA,              &
       lon_ref,lat_ref,     &
       lon, lat,            &
       idx_i, idx_j, hfact, &
       search_limit,        &
       latlon_structure,    &
       lon_1d, lat_1d,      &
       weight_order         )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,  intent(in)  :: npoints                ! number of interpolation point for horizontal
    integer,  intent(in)  :: IA_ref                 ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                 ! number of y-direction    (reference)
    integer,  intent(in)  :: IA                     ! number of x-direction    (target)
    integer,  intent(in)  :: JA                     ! number of y-direction    (target)
    real(RP), intent(in)  :: lon_ref(IA_ref,JA_ref) ! longitude [rad]          (reference)
    real(RP), intent(in)  :: lat_ref(IA_ref,JA_ref) ! latitude  [rad]          (reference)
    real(RP), intent(in)  :: lon  (IA,JA)           ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat  (IA,JA)           ! latitude  [rad]          (target)
    integer,  intent(out) :: idx_i(IA,JA,npoints)   ! i-index in reference     (target)
    integer,  intent(out) :: idx_j(IA,JA,npoints)   ! j-index in reference     (target)
    real(RP), intent(out) :: hfact(IA,JA,npoints)   ! horizontal interp factor (target)

    real(RP), intent(in), optional :: search_limit
    logical,  intent(in), optional :: latlon_structure
    real(RP), intent(in), optional :: lon_1d(IA_ref), lat_1d(JA_ref)
    integer,  intent(in), optional :: weight_order

    logical :: ll_struct_

    real(RP) :: lon_min, lon_max
    real(RP) :: lat_min, lat_max
    real(RP) :: dlon, dlat

    ! for structure grid
    integer  :: psizex, psizey
    real(RP) :: lon0, lat0
    integer, allocatable :: i0(:), i1(:), j0(:), j1(:)

    ! for unstructure grid
    integer :: nsize, psize, nidx_max
    integer, allocatable :: idx_blk(:,:,:), nidx(:,:)
    integer  :: idx_ref(npoints)


    integer  :: i, j, ii, jj, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_fact',3)

    if ( present(latlon_structure) ) then
       ll_struct_ = latlon_structure
    else
       ll_struct_ = .false.
    end if

    hfact(:,:,:) = 0.0_RP


    if ( ll_struct_ ) then

       if ( IA_ref > 10 ) then
          psizex = int( 2.0_RP*sqrt(real(IA_ref+1,RP)) )
       else
          psizex = 1
       end if
       if ( JA_ref > 10 ) then
          psizey = int( 2.0_RP*sqrt(real(JA_ref+1,RP)) )
       else
          psizey = 1
       end if

       allocate( i0(psizex), i1(psizex) )
       allocate( j0(psizey), j1(psizey) )

       dlon = ( lon_max - lon_min ) / psizex
       dlat = ( lat_max - lat_min ) / psizey

       do ii = 1, psizex
          lon0 = lon_min + dlon * (ii-1)
          do i = 1, IA_ref
             if ( lon_1d(i) >= lon0 ) then
                i0(ii) = i
                exit
             end if
          end do
       end do
       do ii = 1, psizex-1
          i1(ii) = i0(ii+1) - 1
       end do
       i1(psizex) = IA_ref

       do jj = 1, psizey
          lat0 = lat_min + dlat * (jj-1)
          do j = 1, JA_ref
             if ( lat_1d(j) >= lat0 ) then
                j0(jj) = j
                exit
             end if
          end do
       end do
       do jj = 1, psizey-1
          j1(jj) = j0(jj+1) - 1
       end do
       j1(psizey) = JA_ref

       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       do j = 1, JA
       do i = 1, IA
          ! main search
          call INTERP_search_horiz_struct( npoints,                     & ! [IN]
                                           psizex, psizey,              & ! [IN]
                                           IA_ref, JA_ref,              & ! [IN]
                                           lon_1d(:), lat_1d(:),        & ! [IN]
                                           lon_min, lat_min,            & ! [IN]
                                           dlon, dlat,                  & ! [IN]
                                           i0(:), i1(:), j0(:), j1(:),  & ! [IN]
                                           lon(i,j), lat(i,j),          & ! [IN]
                                           idx_i(i,j,:), idx_j(i,j,:),  & ! [OUT]
                                           hfact(i,j,:),                & ! [OUT]
                                           search_limit = search_limit, & ! [IN]
                                           weight_order = weight_order  ) ! [IN]
       enddo
       enddo

       deallocate( i0, i1, j0, j1 )

    else

       nsize = IA_ref * JA_ref
       if ( nsize > 100 ) then
          psize = int( sqrt(2.0_RP*sqrt(real(nsize,RP))) )
          nidx_max = nsize / psize * INTERP_buffer_size_fact
       else
          psize = 1
          nidx_max = nsize
       end if

       allocate(idx_blk(nidx_max,psize,psize))
       allocate(nidx   (         psize,psize))

       call INTERP_div_block(nsize, psize, nidx_max,     & ! [IN]
                             lon_ref(:,:), lat_ref(:,:), & ! [IN]
                             idx_blk(:,:,:), nidx(:,:),  & ! [OUT]
                             lon_min, lon_max,           & ! [OUT]
                             lat_min, lat_max,           & ! [OUT]
                             dlon, dlat                  ) ! [OUT]

       !$omp parallel do OMP_SCHEDULE_ collapse(2) &
       !$omp private(idx_ref)
       do j = 1, JA
       do i = 1, IA
          ! main search
          call INTERP_search_horiz( npoints,                     & ! [IN]
                                    nsize,                       & ! [IN]
                                    lon_ref(:,:), lat_ref(:,:),  & ! [IN]
                                    lon_min, lon_max,            & ! [IN]
                                    lat_min, lat_max,            & ! [IN]
                                    psize, nidx_max,             & ! [IN]
                                    dlon, dlat,                  & ! [IN]
                                    idx_blk(:,:,:), nidx(:,:),   & ! [IN]
                                    lon(i,j), lat(i,j),          & ! [IN]
                                    idx_ref(:),                  & ! [OUT]
                                    hfact(i,j,:),                & ! [OUT]
                                    search_limit = search_limit, & ! [IN]
                                    weight_order = weight_order  ) ! [IN]
          do n = 1, npoints
             idx_i(i,j,n) = mod(idx_ref(n) - 1, IA_ref) + 1
             idx_j(i,j,n) = ( idx_ref(n) - 1 ) / IA_ref + 1
          end do
       enddo
       enddo

       deallocate(idx_blk, nidx)

    end if

    call PROF_rapend  ('INTERP_fact',3)

    return
  end subroutine INTERP_factor2d_weight

  !-----------------------------------------------------------------------------
  ! make interpolation factor using bi-linear method for the horizontal direction on latlon grid
  ! This can be used only for structured grid
  subroutine INTERP_factor3d_linear_latlon( &
       KA_ref, KS_ref, KE_ref, &
       IA_ref, JA_ref,         &
       KA, KS, KE,             &
       IA, JA,                 &
       lon_ref, lat_ref,       &
       hgt_ref,                &
       lon, lat,               &
       hgt,                    &
       idx_i, idx_j,           &
       hfact,                  &
       idx_k,                  &
       vfact,                  &
       flag_extrap             )
    implicit none

    integer,  intent(in)  :: KA_ref, KS_ref, KE_ref        ! number of z-direction    (reference)
    integer,  intent(in)  :: IA_ref                        ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                        ! number of y-direction    (reference)
    integer,  intent(in)  :: KA, KS, KE                    ! number of z-direction    (target)
    integer,  intent(in)  :: IA                            ! number of x-direction    (target)
    integer,  intent(in)  :: JA                            ! number of y-direction    (target)
    real(RP), intent(in)  :: lon_ref(IA_ref)               ! longitude                (reference)
    real(RP), intent(in)  :: lat_ref(JA_ref)               ! latitude                 (reference)
    real(RP), intent(in)  :: hgt_ref(KA_ref,IA_ref,JA_ref) ! height    [m]            (reference)
    real(RP), intent(in)  :: lon(IA,JA)                    ! longitude                (target)
    real(RP), intent(in)  :: lat(IA,JA)                    ! latitude                 (target)
    real(RP), intent(in)  :: hgt(KA,IA,JA)                 ! longitude [m]            (target)

    integer,  intent(out) :: idx_i(IA,JA,4)                ! i-index in reference     (target)
    integer,  intent(out) :: idx_j(IA,JA,4)                ! j-index in reference     (target)
    real(RP), intent(out) :: hfact(IA,JA,4)                ! horizontal interp factor (target)
    integer,  intent(out) :: idx_k(KA,2,IA,JA,4)           ! k-index in reference     (target)
    real(RP), intent(out) :: vfact(KA,  IA,JA,4)           ! vertical interp factor   (target)

    logical,  intent(in), optional :: flag_extrap          ! when true, vertical extrapolation will be executed (just copy)

    integer :: i, j, ii, jj, n

    call INTERP_factor2d_linear_latlon( &
         IA_ref, JA_ref,                          & ! [IN]
         IA, JA,                                  & ! [IN]
         lon_ref(:), lat_ref(:),                  & ! [IN]
         lon(:,:), lat(:,:),                      & ! [IN]
         idx_i(:,:,:), idx_j(:,:,:), hfact(:,:,:) ) ! [OUT]

    call PROF_rapstart('INTERP_fact',3)

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(ii,jj)
    do j = 1, JA
    do i = 1, IA
       do n = 1, 4
          ii = idx_i(i,j,n)
          jj = idx_j(i,j,n)
          call INTERP_factor1d( KA_ref, KS_ref, KE_ref,   & ! [IN]
                                KA,     KS,     KE,       & ! [IN]
                                hgt_ref(:,ii,jj),         & ! [IN]
                                hgt    (:,i,j),           & ! [IN]
                                idx_k  (:,:,i,j,n),       & ! [OUT]
                                vfact  (:,  i,j,n),       & ! [OUT]
                                flag_extrap = flag_extrap ) ! [IN, optional]

       enddo
    enddo
    enddo

    call PROF_rapend('INTERP_fact',3)

    return
  end subroutine INTERP_factor3d_linear_latlon

  !-----------------------------------------------------------------------------
  ! make interpolation factor using bi-linear method for the horizontal direction on the xy coordinate
  ! This can be used only for structured grid
  subroutine INTERP_factor3d_linear_xy( &
       KA_ref, KS_ref, KE_ref, &
       IA_ref, JA_ref,         &
       KA, KS, KE,             &
       IA, JA,                 &
       x_ref, y_ref,           &
       hgt_ref,                &
       x, y,                   &
       hgt,                    &
       idx_i, idx_j,           &
       hfact,                  &
       idx_k,                  &
       vfact,                  &
       flag_extrap,            &
       zonal, pole             )
    implicit none
    integer,  intent(in)  :: KA_ref, KS_ref, KE_ref        ! number of z-direction    (reference)
    integer,  intent(in)  :: IA_ref                        ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                        ! number of y-direction    (reference)
    integer,  intent(in)  :: KA, KS, KE                    ! number of z-direction    (target)
    integer,  intent(in)  :: IA                            ! number of x-direction    (target)
    integer,  intent(in)  :: JA                            ! number of y-direction    (target)
    real(RP), intent(in)  :: x_ref(IA_ref,JA_ref)          ! x point                  (reference)
    real(RP), intent(in)  :: y_ref(IA_ref,JA_ref)          ! y point                  (reference)
    real(RP), intent(in)  :: hgt_ref(KA_ref,IA_ref,JA_ref) ! height    [m]            (reference)
    real(RP), intent(in)  :: x  (IA)                       ! x point                  (target)
    real(RP), intent(in)  :: y  (JA)                       ! y point                  (target)
    real(RP), intent(in)  :: hgt(KA,IA,JA)                 ! longitude [m]            (target)

    integer,  intent(out) :: idx_i(IA,JA,4)                ! i-index in reference
    integer,  intent(out) :: idx_j(IA,JA,4)                ! j-index in reference
    real(RP), intent(out) :: hfact(IA,JA,4)                ! horizontal interp factor
    integer,  intent(out) :: idx_k(KA,2,IA,JA,4)           ! k-index in reference
    real(RP), intent(out) :: vfact(KA,  IA,JA,4)           ! vertical interp factor

    logical,  intent(in), optional :: flag_extrap          ! when true, vertical extrapolation will be executed (just copy)
    logical,  intent(in), optional :: zonal
    logical,  intent(in), optional :: pole

    integer :: i, j, ii, jj, n

    call INTERP_factor2d_linear_xy( &
         IA_ref, JA_ref,                           & ! [IN]
         IA, JA,                                   & ! [IN]
         x_ref(:,:), y_ref(:,:),                   & ! [IN]
         x(:), y(:),                               & ! [IN]
         idx_i(:,:,:), idx_j(:,:,:), hfact(:,:,:), & ! [OUT]
         zonal = zonal, pole = pole                ) ! [IN]

    call PROF_rapstart('INTERP_fact',3)

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(ii,jj)
    do j = 1, JA
    do i = 1, IA
       do n = 1, 4
          ii = idx_i(i,j,n)
          jj = idx_j(i,j,n)
          call INTERP_factor1d( KA_ref, KS_ref, KE_ref,   & ! [IN]
                                KA,     KS,     KE,       & ! [IN]
                                hgt_ref(:,ii,jj),         & ! [IN]
                                hgt    (:,i,j),           & ! [IN]
                                idx_k  (:,:,i,j,n),       & ! [OUT]
                                vfact  (:  ,i,j,n),       & ! [OUT]
                                flag_extrap = flag_extrap ) ! [IN, optional]

       enddo
    enddo
    enddo

    call PROF_rapend('INTERP_fact',3)

    return
  end subroutine INTERP_factor3d_linear_xy

  !-----------------------------------------------------------------------------
  ! make interpolation factor using Lat-Lon and Z-Height information
  subroutine INTERP_factor3d_weight( &
       npoints, &
       KA_ref, KS_ref, KE_ref, &
       IA_ref, JA_ref,         &
       KA, KS, KE,             &
       IA, JA,                 &
       lon_ref, lat_ref,       &
       hgt_ref,                &
       lon, lat,               &
       hgt,                    &
       idx_i, idx_j,           &
       hfact,                  &
       idx_k,                  &
       vfact,                  &
       flag_extrap             )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,  intent(in)  :: npoints                       ! number of interpolation point for horizontal
    integer,  intent(in)  :: KA_ref, KS_ref, KE_ref        ! number of z-direction    (reference)
    integer,  intent(in)  :: IA_ref                        ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                        ! number of y-direction    (reference)
    integer,  intent(in)  :: KA, KS, KE                    ! number of z-direction    (target)
    integer,  intent(in)  :: IA                            ! number of x-direction    (target)
    integer,  intent(in)  :: JA                            ! number of y-direction    (target)

    real(RP), intent(in)  :: lon_ref(IA_ref,JA_ref)        ! longitude [rad]          (reference)
    real(RP), intent(in)  :: lat_ref(IA_ref,JA_ref)        ! latitude  [rad]          (reference)
    real(RP), intent(in)  :: hgt_ref(KA_ref,IA_ref,JA_ref) ! height    [m]            (reference)
    real(RP), intent(in)  :: lon  (IA,JA)                  ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat  (IA,JA)                  ! latitude  [rad]          (target)
    real(RP), intent(in)  :: hgt  (KA,IA,JA)               ! longitude [m]            (target)

    integer,  intent(out) :: idx_i(IA,JA,npoints)          ! i-index in reference
    integer,  intent(out) :: idx_j(IA,JA,npoints)          ! j-index in reference
    real(RP), intent(out) :: hfact(IA,JA,npoints)          ! horizontal interp factor
    integer,  intent(out) :: idx_k(KA,2,IA,JA,npoints)     ! k-index in reference
    real(RP), intent(out) :: vfact(KA,  IA,JA,npoints)     ! vertical interp factor

    logical,  intent(in), optional :: flag_extrap          ! when true, vertical extrapolation will be executed (just copy)

    integer :: nsize, psize, nidx_max
    integer, allocatable :: idx_blk(:,:,:), nidx(:,:)
    real(RP) :: lon_min, lon_max
    real(RP) :: lat_min, lat_max
    real(RP) :: dlon, dlat
    integer  :: idx_ref(npoints)

    integer  :: i, j, ii, jj, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_fact',3)

    nsize = IA_ref * JA_ref
    if ( nsize > 100 ) then
       psize = int( sqrt(2.0_RP*sqrt(real(nsize,RP))) )
       nidx_max = nsize / psize * INTERP_buffer_size_fact
    else
       psize = 1
       nidx_max = nsize
    end if

    allocate(idx_blk(nidx_max,psize,psize))
    allocate(nidx   (         psize,psize))

    call INTERP_div_block(nsize, psize, nidx_max,     & ! [IN]
                          lon_ref(:,:), lat_ref(:,:), & ! [IN]
                          idx_blk(:,:,:), nidx(:,:),  & ! [OUT]
                          lon_min, lon_max,           & ! [OUT]
                          lat_min, lat_max,           & ! [OUT]
                          dlon, dlat                  ) ! [OUT]

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(ii,jj,idx_ref)
    do j = 1, JA
    do i = 1, IA

       ! main search
       call INTERP_search_horiz( npoints,                     & ! [IN]
                                 nsize,                       & ! [IN]
                                 lon_ref(:,:), lat_ref(:,:),  & ! [IN]
                                 lon_min, lon_max,            & ! [IN]
                                 lat_min, lat_max,            & ! [IN]
                                 psize, nidx_max,             & ! [IN]
                                 dlon, dlat,                  & ! [IN]
                                 idx_blk(:,:,:), nidx(:,:),   & ! [IN]
                                 lon(i,j), lat(i,j),          & ! [IN]
                                 idx_ref(:), hfact(i,j,:)     ) ! [OUT]

       do n = 1, npoints
          ii = mod(idx_ref(n) - 1, IA_ref) + 1
          jj = ( idx_ref(n) - 1 ) / IA_ref + 1
          idx_i(i,j,n) = ii
          idx_j(i,j,n) = jj
          call INTERP_factor1d( KA_ref, KS_ref, KE_ref,   & ! [IN]
                                KA,     KS,     KE,       & ! [IN]
                                hgt_ref(:,ii,jj),         & ! [IN]
                                hgt    (:,i,j),           & ! [IN]
                                idx_k  (:,:,i,j,n),       & ! [OUT]
                                vfact  (:  ,i,j,n),       & ! [OUT]
                                flag_extrap = flag_extrap ) ! [IN, optional]

       enddo
    enddo
    enddo

    deallocate(idx_blk, nidx)

    call PROF_rapend  ('INTERP_fact',3)

    return
  end subroutine INTERP_factor3d_weight

  !-----------------------------------------------------------------------------
  ! interpolation for 1D data
!OCL SERIAL
  subroutine INTERP_interp1d( &
       KA_ref, KS_ref, KE_ref, &
       KA, KS, KE, &
       idx_k,         &
       vfact,         &
       hgt_ref,       &
       hgt,           &
       val_ref,       &
       val,           &
       logwgt         )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    integer,  intent(in) :: KA_ref, KS_ref, KE_ref ! number of z-direction (reference)
    integer,  intent(in) :: KA, KS, KE             ! number of z-direction (target)

    integer,  intent(in)         :: idx_k(KA,2) ! k-index in reference
    real(RP), intent(in)         :: vfact(KA  ) ! vertical interp factor
    real(RP), intent(in)         :: hgt_ref(KA_ref) ! height (reference)
    real(RP), intent(in)         :: hgt    (KA)     ! height  (target)
    real(RP), intent(in), target :: val_ref(KA_ref) ! value  (reference)

    real(RP), intent(out)        :: val    (KA)     ! value  (target)

    logical,  intent(in), optional :: logwgt !> use logarithmic weighted interpolation?

    logical :: logwgt_

    real(RP), pointer :: work(:)

    integer  :: idx  (KA_ref)
    integer  :: idx_r(KA_ref)
    real(RP) :: FDZ  (KA_ref)
    real(RP) :: U    (KA_ref)

    integer :: kmax
    integer :: k
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_interp',3)

    logwgt_ = .false.
    if ( present(logwgt) ) then
       logwgt_ = logwgt
    endif

    if ( logwgt_ ) then
       allocate( work(KA_ref) )
       do k = KS_ref, KE_ref
          if ( val_ref(k) == UNDEF ) then
             work(k) = UNDEF
          else
             work(k) = log( val_ref(k) )
          end if
       enddo
    else
       work => val_ref
    endif

    call spline_coef( KA_ref, KS_ref, KE_ref, &
                      hgt_ref(:), work(:), & ! (in)
                      kmax,                & ! (out)
                      idx(:), idx_r(:),    & ! (out)
                      U(:), FDZ(:)         ) ! (out)

    call spline_exec( KA_ref, kmax, KA, KS, KE, &
                      idx_k(:,:), vfact(:), & ! (in)
                      hgt_ref(:), hgt(:),   & ! (in)
                      work(:),              & ! (in)
                      idx(:), idx_r(:),     & ! (in)
                      U(:), FDZ(:),         & ! (in)
                      val(:)                ) ! (out)

    if ( logwgt_ ) then
       deallocate( work )
       do k = KS, KE
          if ( val(k) /= UNDEF ) then
             val(k) = exp( val(k) )
          endif
       end do
    endif

    call PROF_rapend  ('INTERP_interp',3)

    return
  end subroutine INTERP_interp1d

  !-----------------------------------------------------------------------------
  ! interpolation for 2D data (nearest-neighbor)
  subroutine INTERP_interp2d( &
       npoints,        &
       IA_ref, JA_ref, &
       IA,  JA,        &
       idx_i, idx_j,   &
       hfact,          &
       val_ref,        &
       val,            &
       threshold_undef )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
    implicit none

    integer,  intent(in)  :: npoints                ! number of interpolation point for horizontal
    integer,  intent(in)  :: IA_ref                 ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                 ! number of y-direction    (reference)
    integer,  intent(in)  :: IA                     ! number of x-direction    (target)
    integer,  intent(in)  :: JA                     ! number of y-direction    (target)
    integer,  intent(in)  :: idx_i  (IA,JA,npoints) ! i-index in reference
    integer,  intent(in)  :: idx_j  (IA,JA,npoints) ! j-index in reference
    real(RP), intent(in)  :: hfact  (IA,JA,npoints) ! horizontal interp factor
    real(RP), intent(in)  :: val_ref(IA_ref,JA_ref) ! value                    (reference)

    real(RP), intent(out) :: val    (IA,JA)         ! value                    (target)

    real(RP), intent(in), optional :: threshold_undef !> return UNDEF if sum of the weight factor is undef the shreshold

    real(RP) :: th_undef

    real(RP) :: fact, valn, f, w
    real(RP) :: sw

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_interp',3)

    th_undef = 0.0_RP
    if ( present(threshold_undef) ) then
       th_undef = threshold_undef
    end if
    th_undef = max( th_undef, EPS * 2.0_RP )

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(fact,valn,f,w,sw)
!OCL PREFETCH
    do j = 1, JA
    do i = 1, IA
       fact = 0.0_RP
       valn = 0.0_RP
       do n = 1, npoints
          f = hfact(i,j,n)
          w = val_ref(idx_i(i,j,n),idx_j(i,j,n))
          if ( f > EPS .and. w .ne. UNDEF ) then
             fact = fact + f
             valn = valn + f * w
          else
             sw = 0.5_RP - sign( 0.5_RP, fact - th_undef + EPS ) ! 1.0 when fact < threshold
             valn = valn / ( fact + sw ) * ( 1.0_RP - sw ) + UNDEF * sw
             exit
          end if
       end do
       val(i,j) = valn
    enddo
    enddo

    call PROF_rapend  ('INTERP_interp',3)

    return
  end subroutine INTERP_interp2d

  !-----------------------------------------------------------------------------
  ! interpolation using one-points for 3D data (nearest-neighbor)
  subroutine INTERP_interp3d( &
       npoints,                &
       KA_ref, KS_ref, KE_ref, &
       IA_ref, JA_ref,         &
       KA, KS, KE,             &
       IA, JA,                 &
       idx_i, idx_j,           &
       hfact,                  &
       idx_k,                  &
       vfact,                  &
       hgt_ref,                &
       hgt,                    &
       val_ref,                &
       val,                    &
       logwgt, threshold_undef )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
    implicit none
    integer,  intent(in) :: npoints                ! number of interpolation point for horizontal
    integer,  intent(in) :: KA_ref, KS_ref, KE_ref ! number of z-direction    (reference)
    integer,  intent(in) :: IA_ref                 ! number of x-direction    (reference)
    integer,  intent(in) :: JA_ref                 ! number of y-direction    (reference)
    integer,  intent(in) :: KA, KS, KE             ! number of z-direction    (target)
    integer,  intent(in) :: IA                     ! number of x-direction    (target)
    integer,  intent(in) :: JA                     ! number of y-direction    (target)

    integer,  intent(in)         :: idx_i  (IA,JA,npoints)        ! i-index in reference
    integer,  intent(in)         :: idx_j  (IA,JA,npoints)        ! j-index in reference
    real(RP), intent(in)         :: hfact  (IA,JA,npoints)        ! horizontal interp factor
    integer,  intent(in)         :: idx_k  (KA,2,IA,JA,npoints)   ! k-index in reference
    real(RP), intent(in)         :: vfact  (KA,  IA,JA,npoints)   ! vertical interp factor
    real(RP), intent(in)         :: hgt_ref(KA_ref,IA_ref,JA_ref) ! height (reference)
    real(RP), intent(in)         :: hgt    (KA,IA,JA)             ! height  (target)
    real(RP), intent(in), target :: val_ref(KA_ref,IA_ref,JA_ref) ! value (reference)

    real(RP), intent(out)        :: val    (KA,IA,JA)             ! value (target)

    logical,  intent(in), optional :: logwgt          !> use logarithmic weighted interpolation?
    real(RP), intent(in), optional :: threshold_undef !> return UNDEF if sum of the weight factor is undef the shreshold

    real(RP) :: th_undef
    logical  :: logwgt_

    real(RP), pointer :: work(:,:,:)

    integer,  allocatable :: kmax(:,:)
    integer,  allocatable :: idx(:,:,:), idx_r(:,:,:)
    real(RP), allocatable :: U(:,:,:), FDZ(:,:,:)

    real(RP) :: valn
    real(RP) :: fact
    real(RP) :: f
    real(RP) :: sw
    real(RP) :: w(KA,npoints)

    integer :: imin, imax
    integer :: jmin, jmax
    integer :: ii, jj
    integer :: k, i, j, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_interp',3)

    th_undef = 0.0_RP
    if ( present(threshold_undef) ) then
       th_undef = threshold_undef
    end if
    th_undef = max( th_undef, EPS * 2.0_RP )


    logwgt_ = .false.
    if ( present(logwgt) ) then
       logwgt_ = logwgt
    endif

    imin = IA_ref
    jmin = JA_ref
    imax = 1
    jmax = 1
    !$omp parallel do OMP_SCHEDULE_ collapse(3) &
    !$omp reduction(min: imin,jmin) &
    !$omp reduction(max: imax,jmax)
    do n = 1, npoints
    do j = 1, JA
    do i = 1, IA
       imin = min(imin, idx_i(i,j,n))
       imax = max(imax, idx_i(i,j,n))
       jmin = min(jmin, idx_j(i,j,n))
       jmax = max(jmax, idx_j(i,j,n))
    end do
    end do
    end do

    allocate( kmax(imin:imax,jmin:jmax) )
    allocate( idx(KA_ref,imin:imax,jmin:jmax), idx_r(KA_ref,imin:imax,jmin:jmax) )
    allocate( U(KA_ref,imin:imax,jmin:jmax), FDZ(KA_ref,imin:imax,jmin:jmax) )

    if ( logwgt_ ) then
       allocate( work(KA_ref,imin:imax,jmin:jmax) )
       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       do j = jmin, jmax
       do i = imin, imax
       do k = KS_ref, KE_ref
          if ( val_ref(k,i,j) == UNDEF ) then
             work(k,i,j) = UNDEF
          else
             work(k,i,j) = log( val_ref(k,i,j) )
          end if
       enddo
       enddo
       enddo
    else
       work => val_ref
    endif

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = jmin, jmax
    do i = imin, imax
       call spline_coef( KA_ref, KS_ref, KE_ref, &
                         hgt_ref(:,i,j), work(:,i,j), & ! (in)
                         kmax(i,j),                   & ! (out)
                         idx(:,i,j), idx_r(:,i,j),    & ! (out)
                         U(:,i,j), FDZ(:,i,j)         ) ! (out)
    end do
    end do

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(valn,fact,w,f,ii,jj,sw)
    do j = 1, JA
    do i = 1, IA

       do n = 1, npoints
          if ( hfact(i,j,n) < EPS ) exit
          ii = idx_i(i,j,n)
          jj = idx_j(i,j,n)
          call spline_exec( KA_ref, kmax(ii,jj), KA, KS, KE, &
                            idx_k(:,:,i,j,n), vfact(:,  i,j,n), & ! (in)
                            hgt_ref(:,ii,jj), hgt(:,i,j),       & ! (in)
                            work(:,ii,jj),                      & ! (in)
                            idx(:,ii,jj), idx_r(:,ii,jj),       & ! (in)
                            U(:,ii,jj), FDZ(:,ii,jj),           & ! (in)
                            w(:,n)                              ) ! (out)
       end do

       do k = KS, KE
          fact = 0.0_RP
          valn = 0.0_RP
          do n = 1, npoints
             f = hfact(i,j,n)
             if ( f > EPS .and. w(k,n) .ne. UNDEF ) then
                fact = fact + f
                valn = valn + f * w(k,n)
             else
                sw = 0.5_RP - sign( 0.5_RP, fact - th_undef ) ! 1.0 when fact < threshold
                valn = valn / ( fact + sw ) * ( 1.0_RP - sw ) + UNDEF * sw
                exit
             endif
          enddo
          val(k,i,j) = valn
       enddo

    enddo
    enddo

    deallocate( kmax, idx, idx_r )
    deallocate( U, FDZ )

    if ( logwgt_ ) then
       deallocate( work )
       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          if ( val(k,i,j) /= UNDEF ) then
             val(k,i,j) = exp( val(k,i,j) )
          end if
       end do
       end do
       end do
    end if

    call PROF_rapend  ('INTERP_interp',3)

    return
  end subroutine INTERP_interp3d


  ! private


  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points
!OCL SERIAL
  subroutine INTERP_search_horiz( &
       npoints,          &
       nsize,            &
       lon_ref, lat_ref, &
       lon_min, lon_max, &
       lat_min, lat_max, &
       psize, nidx_max,  &
       dlon, dlat,       &
       idx_blk, nidx,    &
       lon, lat,         &
       idx_ref,          &
       hfact,            &
       search_limit,     &
       weight_order      )
    use scale_const, only: &
       EPS => CONST_EPS, &
       RADIUS => CONST_RADIUS
    implicit none

    integer,  intent(in)  :: npoints        ! number of interpolation point for horizontal
    integer,  intent(in)  :: nsize          ! number of grids          (reference)
    real(RP), intent(in)  :: lon_ref(nsize) ! longitude [rad]          (reference)
    real(RP), intent(in)  :: lat_ref(nsize) ! latitude  [rad]          (reference)
    real(RP), intent(in)  :: lon_min        ! minmum  longitude [rad]  (reference)
    real(RP), intent(in)  :: lon_max        ! maximum longitude [rad]  (reference)
    real(RP), intent(in)  :: lat_min        ! minmum  latitude  [rad]  (reference)
    real(RP), intent(in)  :: lat_max        ! maximum latitude  [rad]  (reference)
    integer,  intent(in)  :: psize          ! number of blocks for each dimension
    integer,  intent(in)  :: nidx_max       ! maximum number of index in the block
    real(RP), intent(in)  :: dlon           ! block longitude difference
    real(RP), intent(in)  :: dlat           ! block latitude  difference
    integer,  intent(in)  :: idx_blk(nidx_max,psize,psize) ! index of the reference in the block
    integer,  intent(in)  :: nidx            (psize,psize) ! number of indexes in the block
    real(RP), intent(in)  :: lon              ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat              ! latitude  [rad]          (target)
    integer,  intent(out) :: idx_ref(npoints) ! index in reference     (target)
    real(RP), intent(out) :: hfact(npoints)   ! horizontal interp factor (target)

    real(RP), intent(in), optional :: search_limit
    integer,  intent(in), optional :: weight_order

    real(RP) :: drad(npoints)
    real(RP) :: sum
    real(RP) :: search_limit_
    integer  :: weight_order_

    real(RP) :: lon0, lon1, lat0, lat1
    real(RP) :: dlon_sl, dlat_sl
    integer  :: i, n
    integer  :: ii, jj, ii0, jj0
    !---------------------------------------------------------------------------

    if ( present(search_limit) ) then
       search_limit_ = search_limit
    else
       search_limit_ = INTERP_search_limit
    end if
    search_limit_ = search_limit_ / RADIUS ! m to radian

    if ( present(weight_order) ) then
       weight_order_ = weight_order
    else
       weight_order_ = INTERP_weight_order
    end if



    drad   (:) = large_number
    idx_ref(:) = -1

    ! find k-nearest points in the nearest block
    ii0 = max( min( int(( lon - lon_min ) / dlon) + 1, psize ), 1 )
    jj0 = max( min( int(( lat - lat_min ) / dlat) + 1, psize ), 1 )
    do i = 1, nidx(ii0,jj0)
       n = idx_blk(i,ii0,jj0)
       call INTERP_insert( npoints, &
                           lon, lat, &
                           lon_ref(n), lat_ref(n), & ! [IN]
                           n,                      & ! [IN]
                           drad(:), idx_ref(:)     ) ! [INOUT]
    end do

    if ( abs(drad(1)) < EPS ) then
       hfact(:) = 0.0_RP
       hfact(1) = 1.0_RP

       return
    else if ( drad(1) > search_limit_ ) then
       hfact(:) = 0.0_RP
       idx_ref(:) = 1 ! dummy

       return
    end if

    dlon_sl = max(dlon * 0.5_RP, drad(npoints))
    dlat_sl = max(dlat * 0.5_RP, drad(npoints))
    do jj = 1, psize
       lat0 = lat_min + dlat * (jj-1)
       lat1 = lat_min + dlat * jj
       if (     lat <  lat0 - dlat_sl &
           .or. lat >= lat1 + dlat_sl ) cycle
       do ii = 1, psize
          if ( ii==ii0 .and. jj==jj0 ) cycle
          lon0 = lon_min + dlon * (ii-1)
          lon1 = lon_min + dlon * ii
          if (     lon <  lon0 - dlon_sl &
              .or. lon >= lon1 + dlon_sl ) cycle
          do i = 1, nidx(ii,jj)
             n = idx_blk(i,ii,jj)
             call INTERP_insert( npoints, &
                                 lon, lat, &
                                 lon_ref(n), lat_ref(n), & ! [IN]
                                 n,                      & ! [IN]
                                 drad(:), idx_ref(:)     ) ! [INOUT]
          end do
          dlon_sl = max(dlon * 0.5_RP, drad(npoints))
          dlat_sl = max(dlat * 0.5_RP, drad(npoints))
       end do
    end do

    ! factor = 1 / dradian
    if ( weight_order_ < 0 ) then
       do n = 1, npoints
          hfact(n) = 1.0_RP / drad(n)**(-1.0_RP/weight_order_)
       enddo
    else
       do n = 1, npoints
          hfact(n) = 1.0_RP / drad(n)**weight_order_
       enddo
    end if

    ! ignore far point
    do n = 1, npoints
       if ( drad(n) >= search_limit_ ) then
          hfact(n) = 0.0_RP
       endif
    enddo

    ! normalize factor
    sum = 0.0_RP
    do n = 1, npoints
       sum = sum + hfact(n)
    enddo

    if ( sum > 0.0_RP ) then
       do n = 1, npoints
          hfact(n) = hfact(n) / sum
       enddo
    endif

    return
  end subroutine INTERP_search_horiz

  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points for structure grid
!OCL SERIAL
  subroutine INTERP_search_horiz_struct( &
       npoints,          &
       psizex, psizey,   &
       IA_ref, JA_ref,   &
       lon_ref, lat_ref, &
       lon_min, lat_min, &
       dlon, dlat,       &
       i0, i1, j0, j1,   &
       lon, lat,         &
       idx_i, idx_j,     &
       hfact,            &
       search_limit,     &
       weight_order      )
    use scale_const, only: &
       EPS    => CONST_EPS,   &
       RADIUS => CONST_RADIUS
    implicit none

    integer,  intent(in)  :: npoints         ! number of interpolation point for horizontal
    integer,  intent(in)  :: psizex          ! number of block in x-direction (reference)
    integer,  intent(in)  :: psizey          ! number of block in y-direction (reference)
    integer,  intent(in)  :: IA_ref          ! number of grids in x-direction (reference)
    integer,  intent(in)  :: JA_ref          ! number of grids in y-direction (reference)
    real(RP), intent(in)  :: lon_ref(IA_ref) ! longitude [rad]                (reference)
    real(RP), intent(in)  :: lat_ref(JA_ref) ! latitude  [rad]                (reference)
    real(RP), intent(in)  :: lon_min         ! minimum longitude
    real(RP), intent(in)  :: lat_min         ! minimum latitude
    real(RP), intent(in)  :: dlon
    real(RP), intent(in)  :: dlat
    integer,  intent(in)  :: i0(psizex)      ! start index in the block
    integer,  intent(in)  :: i1(psizex)      ! end   index in the block
    integer,  intent(in)  :: j0(psizey)      ! start index in the block
    integer,  intent(in)  :: j1(psizey)      ! start index in the block
    real(RP), intent(in)  :: lon             ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat             ! latitude  [rad]          (target)
    integer,  intent(out) :: idx_i(npoints)  ! x-index in reference     (target)
    integer,  intent(out) :: idx_j(npoints)  ! y-index in reference     (target)
    real(RP), intent(out) :: hfact(npoints)  ! horizontal interp factor (target)

    real(RP), intent(in), optional :: search_limit
    integer,  intent(in), optional :: weight_order

    real(RP) :: drad(npoints)
    real(RP) :: sum
    real(RP) :: search_limit_
    integer  :: weight_order_

    real(RP) :: lon0, lon1, lat0, lat1
    real(RP) :: dlon_sl, dlat_sl

    integer  :: i, j, n
    integer  :: ii, jj, ii0, jj0
    !---------------------------------------------------------------------------

    if ( present(search_limit) ) then
       search_limit_ = search_limit
    else
       search_limit_ = INTERP_search_limit
    end if
    search_limit_ = search_limit_ / RADIUS ! m to radian

    if ( present(weight_order) ) then
       weight_order_ = weight_order
    else
       weight_order_ = INTERP_weight_order
    end if


    drad (:) = large_number
    idx_i(:) = 1
    idx_j(:) = 1

    ! find k-nearest points in the nearest block
    ii0 = max( min( int(( lon - lon_min ) / dlon) + 1, psizex ), 1)
    jj0 = max( min( int(( lat - lat_min ) / dlat) + 1, psizey ), 1)
    do j = j0(jj0), j1(jj0)
    do i = i0(ii0), i1(ii0)
       call INTERP_insert_2d( npoints, &
                              lon, lat, &
                              lon_ref(i), lat_ref(j),     & ! [IN]
                              i, j,                       & ! [IN]
                              drad(:), idx_i(:), idx_j(:) ) ! [INOUT]
    end do
    end do

    if ( abs(drad(1)) < EPS ) then
       hfact(:) = 0.0_RP
       hfact(1) = 1.0_RP

       return
    else if ( drad(1) > search_limit_ ) then
       hfact(:) = 0.0_RP
       idx_i(:) = 1 ! dummy
       idx_j(:) = 1 ! dummy

       return
    end if

    dlon_sl = max(dlon * 0.5_RP, drad(npoints))
    dlat_sl = max(dlat * 0.5_RP, drad(npoints))
    do jj = 1, psizey
       lat0 = lat_min + dlat * (jj-1)
       lat1 = lat_min + dlat * jj
       if (     lat <  lat0 - dlat_sl &
           .or. lat >= lat1 + dlat_sl ) cycle
       do ii = 1, psizex
          if ( ii==ii0 .and. jj==jj0 ) cycle
          lon0 = lon_min + dlon * (ii-1)
          lon1 = lon_min + dlon * ii
          if (     lon <  lon0 - dlon_sl &
              .or. lon >= lon1 + dlon_sl ) cycle
          do j = j0(jj0), j1(jj0)
          do i = i0(ii0), i1(ii0)
             call INTERP_insert_2d( npoints, &
                                    lon, lat, &
                                    lon_ref(i), lat_ref(j),     & ! [IN]
                                    i, j,                       & ! [IN]
                                    drad(:), idx_i(:), idx_j(:) ) ! [INOUT]
          end do
          end do
          dlon_sl = max(dlon * 0.5_RP, drad(npoints))
          dlat_sl = max(dlat * 0.5_RP, drad(npoints))
       end do
    end do


    ! factor = 1 / drad
    if ( weight_order_ < 0 ) then
       do n = 1, npoints
          hfact(n) = 1.0_RP / drad(n)**(-1.0_RP/weight_order_)
       enddo
    else
       do n = 1, npoints
          hfact(n) = 1.0_RP / drad(n)**weight_order_
       enddo
    end if

    ! ignore far point
    do n = 1, npoints
       if ( drad(n) >= search_limit_ ) then
          hfact(n) = 0.0_RP
       endif
    enddo

    ! normalize factor
    sum = 0.0_RP
    do n = 1, npoints
       sum = sum + hfact(n)
    enddo

    if ( sum > 0.0_RP ) then
       do n = 1, npoints
          hfact(n) = hfact(n) / sum
       enddo
    endif

    return
  end subroutine INTERP_search_horiz_struct


!OCL SERIAL
  subroutine INTERP_insert( npoints, &
                            lon, lat, &
                            lon_ref, lat_ref, &
                            i,                &
                            drad, idx_i       )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_sort, only: &
       SORT_exec

    integer,  intent(in) :: npoints
    real(RP), intent(in) :: lon, lat
    real(RP), intent(in) :: lon_ref, lat_ref
    integer,  intent(in) :: i

    real(RP), intent(inout) :: drad(npoints)
    integer,  intent(inout) :: idx_i(npoints)

    real(RP) :: dradian

    if ( lon_ref == UNDEF ) return

    dradian = haversine( lon, lat, lon_ref, lat_ref )

    if ( dradian <= drad(npoints) ) then
       ! replace last(=longest) value
       drad (npoints) = dradian
       idx_i(npoints) = i

       ! sort by ascending order
       call SORT_exec( npoints,   & ! [IN]
                       drad (:),  & ! [INOUT]
                       idx_i(:)   ) ! [INOUT]
    endif

    return
  end subroutine INTERP_insert

!OCL SERIAL
  subroutine INTERP_insert_2d( npoints, &
                               lon, lat, &
                               lon_ref, lat_ref,  &
                               i, j,              &
                               drad, idx_i, idx_j )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_sort, only: &
       SORT_exec

    integer,  intent(in) :: npoints
    real(RP), intent(in) :: lon, lat
    real(RP), intent(in) :: lon_ref, lat_ref
    integer,  intent(in) :: i, j

    real(RP), intent(inout) :: drad(npoints)
    integer,  intent(inout) :: idx_i(npoints)
    integer,  intent(inout) :: idx_j(npoints)

    real(RP) :: dradian

    if ( lon_ref == UNDEF ) return

    dradian = haversine( lon, lat, lon_ref, lat_ref )

    if ( dradian <= drad(npoints) ) then
       ! replace last(=longest) value
       drad (npoints) = dradian
       idx_i(npoints) = i
       idx_j(npoints) = j

       ! sort by ascending order
       call SORT_exec( npoints,   & ! [IN]
                       drad (:),  & ! [INOUT]
                       idx_i(:),  & ! [INOUT]
                       idx_j(:)   ) ! [INOUT]
    endif

    return
  end subroutine INTERP_insert_2d


  subroutine INTERP_div_block(nsize, psize, nidx_max,  &
                              lon_ref, lat_ref, &
                              idx, nidx,        &
                              lon_min, lon_max, &
                              lat_min, lat_max, &
                              dlon, dlat        )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_prc, only: &
       PRC_abort
    integer,  intent(in) :: nsize
    integer,  intent(in) :: psize
    integer,  intent(in) :: nidx_max
    real(RP), intent(in) :: lon_ref(nsize)
    real(RP), intent(in) :: lat_ref(nsize)

    integer,  intent(out) :: idx (nidx_max,psize,psize)
    integer,  intent(out) :: nidx(psize,psize)
    real(RP), intent(out) :: lon_min, lon_max
    real(RP), intent(out) :: lat_min, lat_max
    real(RP), intent(out) :: dlon, dlat

    integer :: i, ii, jj, n

    lon_min = minval(lon_ref(:), mask=lon_ref.ne.UNDEF)
    lon_max = maxval(lon_ref(:), mask=lon_ref.ne.UNDEF)
    lat_min = minval(lat_ref(:), mask=lat_ref.ne.UNDEF)
    lat_max = maxval(lat_ref(:), mask=lat_ref.ne.UNDEF)

    dlon = ( lon_max - lon_min ) / psize
    dlat = ( lat_max - lat_min ) / psize

    nidx(:,:) = 0
    do i = 1, nsize
       if ( lon_ref(i) == UNDEF ) cycle
       ii = min(int((lon_ref(i) - lon_min) / dlon) + 1, psize)
       jj = min(int((lat_ref(i) - lat_min) / dlat) + 1, psize)
       n = nidx(ii,jj) + 1
       nidx(ii,jj) = n
       if ( n <= nidx_max ) idx(n,ii,jj) = i
    end do

    if ( maxval(nidx) > nidx_max ) then
       LOG_ERROR("INTERP_search_horiz",*) 'Buffer size is not enough'
       LOG_ERROR_CONT(*)                  '   Use larger INTERP_buffer_size_fact: ', INTERP_buffer_size_fact * maxval(nidx)/nidx_max
       call PRC_abort
    end if

    return
  end subroutine INTERP_div_block

  !-----------------------------------------------------------------------------
  subroutine INTERP_bilinear_inv( &
       x_ref0, x_ref1, x_ref2, x_ref3, &
       y_ref0, y_ref1, y_ref2, y_ref3, &
       x, y,                           &
       u, v,                           &
       error                           )
    implicit none
    real(RP), intent(in) :: x_ref0, x_ref1, x_ref2, x_ref3
    real(RP), intent(in) :: y_ref0, y_ref1, y_ref2, y_ref3
    real(RP), intent(in) :: x, y

    real(RP), intent(out) :: u, v
    logical,  intent(out) :: error

    real(RP), parameter :: EPS = 1E-6_RP

    real(RP) :: e_x, e_y
    real(RP) :: f_x, f_y
    real(RP) :: g_x, g_y
    real(RP) :: h_x, h_y
    real(RP) :: k0, k1, k2
    real(RP) :: w
    real(RP) :: sig

    e_x = x_ref1 - x_ref0
    e_y = y_ref1 - y_ref0

    f_x = x_ref3 - x_ref0
    f_y = y_ref3 - y_ref0

    g_x = x_ref0 - x_ref1 + x_ref2 - x_ref3
    g_y = y_ref0 - y_ref1 + y_ref2 - y_ref3

    h_x = x - x_ref0
    h_y = y - y_ref0

    k0 = cross(h_x, h_y, e_x, e_y)
    k1 = cross(e_x, e_y, f_x, f_y) + cross(h_x, h_y, g_x, g_y)
    k2 = cross(g_x, g_y, f_x, f_y)

    w = k1**2 - 4.0_RP * k0 * k2
    if ( w < 0.0_RP ) then
       LOG_ERROR("INTERP_bilinear_inv",*) 'Outside of the quadrilateral'
       error = .true.
       return
    end if

    if ( abs(k1) < EPS ) then
       LOG_ERROR("INTERP_bilinear_inv",*) 'Unexpected error occured', k1
       LOG_ERROR_CONT(*) x_ref0, x_ref1, x_ref2, x_ref3
       LOG_ERROR_CONT(*) y_ref0, y_ref1, y_ref2, y_ref3
       LOG_ERROR_CONT(*) x, y
       error = .true.
       return
    end if
    if ( abs(k2) < EPS * sqrt( (x_ref2-x_ref0)**2+(y_ref2-y_ref0)**2 ) ) then
       v = - k0 / k1
    else
       sig = sign( 1.0_RP, cross(x_ref1-x_ref0, y_ref1-y_ref0, x_ref3-x_ref0, y_ref3-y_ref0) )
       v = ( - k1 + sig * sqrt(w) ) / ( 2.0_RP * k2 )
    end if
    u = ( h_x - f_x * v ) / ( e_x + g_x * v )

    if ( u < -EPS .or. u > 1.0_RP+EPS .or. v < -EPS .or. v > 1.0_RP+EPS ) then
       LOG_ERROR("INTERP_bilinear_inv",*) 'Unexpected error occured', u, v
       LOG_ERROR_CONT(*) x_ref0, x_ref1, x_ref2, x_ref3
       LOG_ERROR_CONT(*) y_ref0, y_ref1, y_ref2, y_ref3
       LOG_ERROR_CONT(*) x, y
       LOG_ERROR_CONT(*) k0, k1, k2, sig
       error = .true.
       return
    end if
    u = min( max( u, 0.0_RP ), 1.0_RP )
    v = min( max( v, 0.0_RP ), 1.0_RP )

    error = .false.

    return
  end subroutine INTERP_bilinear_inv

  !-----------------------------------------------------------------------------
  ! check whether the point is inside of the quadrilateral
  subroutine INTERP_check_inside( &
       x_ref0, x_ref1, x_ref2, x_ref3, &
       y_ref0, y_ref1, y_ref2, y_ref3, &
       x, y,                           &
       inc_i, inc_j,                   &
       error                           )
    implicit none
    real(RP), intent(in) :: x_ref0, x_ref1, x_ref2, x_ref3
    real(RP), intent(in) :: y_ref0, y_ref1, y_ref2, y_ref3
    real(RP), intent(in) :: x, y

    integer,  intent(out) :: inc_i, inc_j
    logical,  intent(out) :: error

    real(RP) :: sig
    real(RP) :: c1, c2, c3, c4
    logical :: fx, fy

    error = .false.

    sig = sign( 1.0_RP, cross(x_ref1-x_ref0, y_ref1-y_ref0, x_ref3-x_ref0, y_ref3-y_ref0) )

    c1 = sig * cross(x_ref1-x_ref0, y_ref1-y_ref0, x-x_ref0, y-y_ref0)
    c2 = sig * cross(x_ref2-x_ref1, y_ref2-y_ref1, x-x_ref1, y-y_ref1)
    c3 = sig * cross(x_ref3-x_ref2, y_ref3-y_ref2, x-x_ref2, y-y_ref2)
    c4 = sig * cross(x_ref0-x_ref3, y_ref0-y_ref3, x-x_ref3, y-y_ref3)

    ! if all the c1 - c4 are positive, the point is inside the quadrilateral
    inc_i = 0
    inc_j = 0
    if ( c1 < 0.0_RP ) inc_j = -1
    if ( c2 < 0.0_RP ) inc_i =  1
    if ( c3 < 0.0_RP ) inc_j =  1
    if ( c4 < 0.0_RP ) inc_i = -1

    fx = c2 < 0.0_RP .and. c4 < 0.0_RP
    fy = c1 < 0.0_RP .and. c3 < 0.0_RP
    if ( fx .and. fy ) then
       LOG_ERROR("INTERP_check_inside",*) 'Unexpected error occured', c1, c2, c3, c4
       LOG_ERROR_CONT(*) x_ref0, x_ref1, x_ref2, x_ref3
       LOG_ERROR_CONT(*) y_ref0, y_ref1, y_ref2, y_ref3
       LOG_ERROR_CONT(*) x, y
       error = .true.
       return
    else if ( fx ) then
       inc_i = 0
    else if ( fy ) then
       inc_j = 0
    end if

    return
  end subroutine INTERP_check_inside

  !-----------------------------------------------------------------------------
  !> cross product
  function cross( &
       x0, y0, &
       x1, y1 )
    implicit none
    real(RP), intent(in) :: x0, y0
    real(RP), intent(in) :: x1, y1
    real(RP) :: cross

    cross = x0 * y1 - x1 * y0

  end function cross

  !-----------------------------------------------------------------------------
  ! Haversine Formula (from R.W. Sinnott, "Virtues of the Haversine",
  ! Sky and Telescope, vol. 68, no. 2, 1984, p. 159):
  function haversine( &
       lon0, lat0, &
       lon1, lat1  ) &
       result( d )
    implicit none

    real(RP), intent(in) :: lon0   ! [rad]
    real(RP), intent(in) :: lat0   ! [rad]
    real(RP), intent(in) :: lon1   ! [rad]
    real(RP), intent(in) :: lat1   ! [rad]
    real(RP)             :: d      ! [rad]

    real(RP) :: dlonh, dlath
    real(RP) :: work1
    !---------------------------------------------------------------------------

    dlonh = 0.5_RP * ( lon0 - lon1 )
    dlath = 0.5_RP * ( lat0 - lat1 )

    work1 = sin(dlath)**2 + cos(lat0) * cos(lat1) * sin(dlonh)**2

    d = 2.0_RP * asin( min(sqrt(work1),1.0_RP) )

  end function haversine


!OCL SERIAL
  subroutine spline_coef( &
       KA_ref, KS_ref, KE_ref, &
       hgt_ref, val_ref, &
       kmax,             &
       idx, idx_r,       &
       U, FDZ            )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none
    integer,  intent(in) :: KA_ref, KS_ref, KE_ref

    real(RP), intent(in) :: hgt_ref(KA_ref)
    real(RP), intent(in) :: val_ref(KA_ref)

    integer,  intent(out) :: kmax
    integer,  intent(out) :: idx  (KA_ref)
    integer,  intent(out) :: idx_r(KA_ref)
    real(RP), intent(out) :: U  (KA_ref)
    real(RP), intent(out) :: FDZ(KA_ref)

    real(RP) :: MD(KA_ref)
    real(RP) :: V (KA_ref)
    real(RP) :: dz
    integer  :: k

    if ( INTERP_use_spline_vert ) then

       do k = KS_ref, KE_ref-1
          if ( val_ref(k) .ne. UNDEF ) then
             idx(1) = k
             idx_r(k) = 1
             exit
          end if
       end do
       kmax = 1
       FDZ(1) = 1e10 ! dummy
       do k = idx(1)+1, KE_ref
          dz = hgt_ref(k) - hgt_ref(idx(kmax))
          if ( val_ref(k) .ne. UNDEF .and. dz > EPS ) then
             do while ( kmax > 1 .and. FDZ(kmax) < dz * 0.1_RP )
                kmax = kmax - 1 ! marge
             end do
             kmax = kmax + 1
             idx(kmax) = k
             if ( idx(kmax-1)+1 <= k-1 ) idx_r(idx(kmax-1)+1:k-1) = kmax-1
             idx_r(k) = kmax
             FDZ(kmax) = hgt_ref(k) - hgt_ref(idx(kmax-1))
          end if
       end do

       if ( kmax > 3 ) then

          MD(2) = 2.0_RP * ( FDZ(2) + FDZ(3) ) + FDZ(2)
          do k = 3, kmax-2
             MD(k) = 2.0_RP * ( FDZ(k) + FDZ(k+1) )
          end do
          MD(kmax-1) = 2.0_RP * ( FDZ(kmax-1) + FDZ(kmax) ) + FDZ(kmax)

          do k = 2, kmax-1
             V(k) = ( val_ref(idx(k+1)) - val_ref(idx(k  )) ) / FDZ(k+1) &
                  - ( val_ref(idx(k  )) - val_ref(idx(k-1)) ) / FDZ(k  )
          end do

          call MATRIX_SOLVER_tridiagonal( kmax, 2, kmax-1, &
                                          FDZ(2:), MD(:), FDZ(:), & ! (in)
                                          V(:),                   & ! (in)
                                          U(:)                    ) ! (out)
!          U(1) = 0.0_RP
!          U(kmax) = 0.0_RP
          U(1) = U(2)
          U(kmax) = U(kmax-1)

       else

          idx(kmax) = idx(1) ! force linear interpolateion

       end if

    else

       kmax = 1
       idx(1) = -999 ! dummy

    end if

    return
  end subroutine spline_coef

!OCL SERIAL
  subroutine spline_exec( &
       KA_ref, kmax, KA, KS, KE, &
       idx_k, vfact, &
       hgt_ref, hgt, &
       val_ref,      &
       idx, idx_r,   &
       U, FDZ,       &
       val           )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    integer,  intent(in) :: KA_ref, kmax
    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: idx_k(KA,2)
    real(RP), intent(in) :: vfact(KA)
    real(RP), intent(in) :: hgt_ref(KA_ref)
    real(RP), intent(in) :: hgt    (KA)
    real(RP), intent(in) :: val_ref(KA_ref)
    integer,  intent(in) :: idx    (KA_ref)
    integer,  intent(in) :: idx_r  (KA_ref)
    real(RP), intent(in) :: U      (KA_ref)
    real(RP), intent(in) :: FDZ    (KA_ref)

    real(RP), intent(out) :: val(KA)

    real(RP) :: c1, c2, c3, d
    integer  :: k, kk, kk2

    do k = KS, KE
       kk = idx_k(k,1)
       if ( kk == -1 ) then
          val(k) = UNDEF
       else if ( idx_k(k,2) == -1 ) then
          val(k) = val_ref(kk)
       else if ( kk < idx(1) .or. kk >= idx(kmax) ) then ! linear interpolation
          if ( val_ref(kk) == UNDEF .or. val_ref(kk+1) == UNDEF ) then
             val(k) = UNDEF
          else
             val(k) = val_ref(kk) * vfact(k) + val_ref(kk+1) * ( 1.0_RP - vfact(k) )
          end if
       else
          kk2 = idx_r(kk)
          kk = idx(kk2)
          c3 = ( U(kk2+1) - U(kk2) ) / FDZ(kk2+1)
          c2 = 3.0_RP * U(kk2)
          c1 = ( val_ref(idx(kk2+1)) - val_ref(kk) ) / FDZ(kk2+1) - ( 2.0_RP * U(kk2) + U(kk2+1) ) * FDZ(kk2+1)
          d = hgt(k) - hgt_ref(kk)

          val(k) = ( ( c3 * d + c2 ) * d + c1 ) * d + val_ref(kk)
       end if
    end do

    return
  end subroutine spline_exec

end module scale_interp
