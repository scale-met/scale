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
  use scale_debug
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: INTERP_setup
  public :: INTERP_domain_compatibility
  public :: INTERP_factor2d
  public :: INTERP_factor3d
  public :: INTERP_interp2d
  public :: INTERP_interp3d

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: INTERP_search_nearest_block
  private :: INTERP_search_horiz
  private :: INTERP_search_vert

  private :: haversine

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: large_number = 9.999E+15_RP

  integer,  private :: INTERP_divnum       = 10
  integer,  private :: INTERP_weight_order = 2
  real(RP), private :: INTERP_search_limit

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTERP_setup( &
       divnum,       &
       weight_order, &
       search_limit  )
    implicit none

    integer,  intent(in) :: divnum
    integer,  intent(in) :: weight_order
    real(RP), intent(in), optional :: search_limit
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("INTERP_setup",*) 'Setup'

    INTERP_divnum       = divnum
    INTERP_weight_order = weight_order

    INTERP_search_limit = large_number
    if ( present(search_limit) ) then
       INTERP_search_limit = search_limit
       LOG_INFO("INTERP_setup",*) 'search limit [m] : ', INTERP_search_limit
    endif

    return
  end subroutine INTERP_setup

  !-----------------------------------------------------------------------------
  subroutine INTERP_domain_compatibility( &
       lon_org, &
       lat_org, &
       lev_org, &
       lon_loc, &
       lat_loc, &
       lev_loc, &
       skip_x,  &
       skip_y,  &
       skip_z   )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none

    real(RP), intent(in) :: lon_org(:,:)
    real(RP), intent(in) :: lat_org(:,:)
    real(RP), intent(in) :: lev_org(:,:,:)
    real(RP), intent(in) :: lon_loc(:,:)
    real(RP), intent(in) :: lat_loc(:,:)
    real(RP), intent(in) :: lev_loc(:,:,:)

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
       max_ref = maxval( lev_org(:,:,:) )
       max_loc = maxval( lev_loc(:,:,:) ) ! HALO + 1

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
  ! make interpolation factor using Lat-Lon
  subroutine INTERP_factor2d( &
       npoints, &
       IA_ref,  &
       JA_ref,  &
       lon_ref, &
       lat_ref, &
       IA,      &
       JA,      &
       lon,     &
       lat,     &
       idx_i,   &
       idx_j,   &
       hfact,   &
       divnum   )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer,  intent(in)  :: npoints                ! number of interpolation point for horizontal
    integer,  intent(in)  :: IA_ref                 ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                 ! number of y-direction    (reference)
    real(RP), intent(in)  :: lon_ref(IA_ref,JA_ref) ! longitude [rad]          (reference)
    real(RP), intent(in)  :: lat_ref(IA_ref,JA_ref) ! latitude  [rad]          (reference)
    integer,  intent(in)  :: IA                     ! number of x-direction    (target)
    integer,  intent(in)  :: JA                     ! number of y-direction    (target)
    real(RP), intent(in)  :: lon  (IA,JA)           ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat  (IA,JA)           ! latitude  [rad]          (target)
    integer,  intent(out) :: idx_i(IA,JA,npoints)   ! i-index in reference     (target)
    integer,  intent(out) :: idx_j(IA,JA,npoints)   ! j-index in reference     (target)
    real(RP), intent(out) :: hfact(IA,JA,npoints)   ! horizontal interp factor (target)
    integer,  intent(in), optional :: divnum

    integer  :: IS, IE ! [start,end] index for x-direction
    integer  :: JS, JE ! [start,end] index for y-direction

    integer  :: i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_fact',3)

    hfact(:,:,:) = 0.0_RP

    !$omp parallel do OMP_SCHEDULE_ collapse(2), &
    !$omp private(IS,IE,JS,JE)
    do j = 1, JA
    do i = 1, IA
       ! nearest block search
       call INTERP_search_nearest_block( IA_ref, JA_ref, & ! [IN]
                                         lon_ref(:,:),   & ! [IN]
                                         lat_ref(:,:),   & ! [IN]
                                         lon    (i,j),   & ! [IN]
                                         lat    (i,j),   & ! [IN]
                                         IS, IE,         & ! [OUT]
                                         JS, JE,         & ! [OUT]
                                         divnum = divnum ) ! [IN]

       ! main search
       call INTERP_search_horiz( npoints,        & ! [IN]
                                 IA_ref, JA_ref, & ! [IN]
                                 IS, IE,         & ! [IN]
                                 JS, JE,         & ! [IN]
                                 lon_ref(:,:),   & ! [IN]
                                 lat_ref(:,:),   & ! [IN]
                                 lon    (i,j),   & ! [IN]
                                 lat    (i,j),   & ! [IN]
                                 idx_i  (i,j,:), & ! [OUT]
                                 idx_j  (i,j,:), & ! [OUT]
                                 hfact  (i,j,:)  ) ! [OUT]
    enddo
    enddo

    call PROF_rapend  ('INTERP_fact',3)

    return
  end subroutine INTERP_factor2d

  !-----------------------------------------------------------------------------
  ! make interpolation factor using Lat-Lon and Z-Height information
  subroutine INTERP_factor3d( &
       npoints, &
       KA_ref,  &
       KS_ref,  &
       KE_ref,  &
       IA_ref,  &
       JA_ref,  &
       lon_ref, &
       lat_ref, &
       hgt_ref, &
       KA,      &
       KS,      &
       KE,      &
       IA,      &
       JA,      &
       lon,     &
       lat,     &
       hgt,     &
       idx_i,   &
       idx_j,   &
       hfact,   &
       idx_k,   &
       vfact,   &
       divnum   )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer,  intent(in)  :: npoints                       ! number of interpolation point for horizontal
    integer,  intent(in)  :: KA_ref, KS_ref, KE_ref        ! number of z-direction    (reference)
    integer,  intent(in)  :: IA_ref                        ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                        ! number of y-direction    (reference)
    real(RP), intent(in)  :: lon_ref(IA_ref,JA_ref)        ! longitude [rad]          (reference)
    real(RP), intent(in)  :: lat_ref(IA_ref,JA_ref)        ! latitude  [rad]          (reference)
    real(RP), intent(in)  :: hgt_ref(KA_ref,IA_ref,JA_ref) ! height    [m]            (reference)
    integer,  intent(in)  :: KA, KS, KE                    ! number of z-direction    (target)
    integer,  intent(in)  :: IA                            ! number of x-direction    (target)
    integer,  intent(in)  :: JA                            ! number of y-direction    (target)
    real(RP), intent(in)  :: lon  (IA,JA)                  ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat  (IA,JA)                  ! latitude  [rad]          (target)
    real(RP), intent(in)  :: hgt  (KA,IA,JA)               ! longitude [m]            (target)
    integer,  intent(out) :: idx_i(IA,JA,npoints)          ! i-index in reference     (target)
    integer,  intent(out) :: idx_j(IA,JA,npoints)          ! j-index in reference     (target)
    real(RP), intent(out) :: hfact(IA,JA,npoints)          ! horizontal interp factor (target)
    integer,  intent(out) :: idx_k(KA,2,IA,JA,npoints)     ! i-index in reference     (target)
    real(RP), intent(out) :: vfact(KA,2,IA,JA,npoints)     ! horizontal interp factor (target)
    integer,  intent(in), optional :: divnum

    integer  :: IS, IE ! [start,end] index for x-direction
    integer  :: JS, JE ! [start,end] index for y-direction

    integer  :: k, i, j, ii, jj, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_fact',3)

    hfact(:,:,:)     = 0.0_RP
    vfact(:,:,:,:,:) = 0.0_RP

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(ii,jj,IS,IE,JS,JE)
    do j = 1, JA
    do i = 1, IA
       ! nearest block search
       call INTERP_search_nearest_block( IA_ref, JA_ref, & ! [IN]
                                         lon_ref(:,:),   & ! [IN]
                                         lat_ref(:,:),   & ! [IN]
                                         lon    (i,j),   & ! [IN]
                                         lat    (i,j),   & ! [IN]
                                         IS, IE,         & ! [OUT]
                                         JS, JE,         & ! [OUT]
                                         divnum = divnum ) ! [IN]

       ! main search
       call INTERP_search_horiz( npoints,        & ! [IN]
                                 IA_ref, JA_ref, & ! [IN]
                                 IS, IE,         & ! [IN]
                                 JS, JE,         & ! [IN]
                                 lon_ref(:,:),   & ! [IN]
                                 lat_ref(:,:),   & ! [IN]
                                 lon    (i,j),   & ! [IN]
                                 lat    (i,j),   & ! [IN]
                                 idx_i  (i,j,:), & ! [OUT]
                                 idx_j  (i,j,:), & ! [OUT]
                                 hfact  (i,j,:)  ) ! [OUT]

       do n = 1, npoints
          ii = idx_i(i,j,n)
          jj = idx_j(i,j,n)

          call INTERP_search_vert( KA_ref, KS_ref, KE_ref, & ! [IN]
                                   KA,     KS,     KE,     & ! [IN]
                                   hgt_ref(:,ii,jj),       & ! [IN]
                                   hgt    (:,i,j),         & ! [IN]
                                   idx_k  (:,:,i,j,n),     & ! [OUT]
                                   vfact  (:,:,i,j,n)      ) ! [OUT]
       enddo
    enddo
    enddo

    call PROF_rapend  ('INTERP_fact',3)

    return
  end subroutine INTERP_factor3d

  !-----------------------------------------------------------------------------
  ! interpolation using one-points for 2D data (nearest-neighbor)
  subroutine INTERP_interp2d( &
       npoints, &
       IA_ref,  &
       JA_ref,  &
       IA,      &
       JA,      &
       idx_i,   &
       idx_j,   &
       hfact,   &
       val_ref, &
       val      )
    use scale_const, only: &
       EPS   => CONST_EPS,  &
       UNDEF => CONST_UNDEF
    implicit none

    integer,  intent(in)  :: npoints                ! number of interpolation point for horizontal
    integer,  intent(in)  :: IA_ref                 ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                 ! number of y-direction    (reference)
    integer,  intent(in)  :: IA                     ! number of x-direction    (target)
    integer,  intent(in)  :: JA                     ! number of y-direction    (target)
    integer,  intent(in)  :: idx_i  (IA,JA,npoints) ! i-index in reference     (target)
    integer,  intent(in)  :: idx_j  (IA,JA,npoints) ! j-index in reference     (target)
    real(RP), intent(in)  :: hfact  (IA,JA,npoints) ! horizontal interp factor (target)
    real(RP), intent(in)  :: val_ref(IA_ref,JA_ref) ! value                    (reference)
    real(RP), intent(out) :: val    (IA,JA)         ! value                    (target)

    real(RP) :: sw
    integer  :: i, j, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_interp',3)

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(sw)
!OCL PREFETCH
    do j = 1, JA
    do i = 1, IA
       sw = 0.5_RP + sign(0.5_RP,hfact(i,j,1)-EPS)

       val(i,j) = (        sw ) * hfact(i,j,1) * val_ref(idx_i(i,j,1),idx_j(i,j,1)) &
                + ( 1.0_RP-sw ) * UNDEF
    enddo
    enddo

!OCL SERIAL
    do n = 2, npoints
    !$omp parallel do OMP_SCHEDULE_ collapse(2)
!OCL PREFETCH
    do j = 1, JA
    do i = 1, IA
       val(i,j) = val(i,j) &
                + hfact(i,j,n) * val_ref(idx_i(i,j,n),idx_j(i,j,n))
    enddo
    enddo
    enddo

    call PROF_rapend  ('INTERP_interp',3)

    return
  end subroutine INTERP_interp2d

  !-----------------------------------------------------------------------------
  ! interpolation using one-points for 3D data (nearest-neighbor)
  subroutine INTERP_interp3d( &
       npoints,    &
       KA_ref,     &
       IA_ref,     &
       JA_ref,     &
       KA, KS, KE, &
       IA,         &
       JA,         &
       idx_i,      &
       idx_j,      &
       hfact,      &
       idx_k,      &
       vfact,      &
       val_ref,    &
       val,        &
       logwgt      )
    implicit none

    integer,  intent(in) :: npoints                       ! number of interpolation point for horizontal
    integer,  intent(in) :: KA_ref                        ! number of x-direction    (reference)
    integer,  intent(in) :: IA_ref                        ! number of x-direction    (reference)
    integer,  intent(in) :: JA_ref                        ! number of y-direction    (reference)
    integer,  intent(in) :: KA, KS, KE                    ! number of z-direction    (target)
    integer,  intent(in) :: IA                            ! number of x-direction    (target)
    integer,  intent(in) :: JA                            ! number of y-direction    (target)
    integer,  intent(in) :: idx_i  (IA,JA,npoints)        ! i-index in reference     (target)
    integer,  intent(in) :: idx_j  (IA,JA,npoints)        ! j-index in reference     (target)
    real(RP), intent(in) :: hfact  (IA,JA,npoints)        ! horizontal interp factor (target)
    integer,  intent(in) :: idx_k  (KA,2,IA,JA,npoints)   ! k-index in reference     (target)
    real(RP), intent(in) :: vfact  (KA,2,IA,JA,npoints)   ! vertical interp factor   (target)

    real(RP), intent(in), target :: val_ref(KA_ref,IA_ref,JA_ref) ! value (reference)

    real(RP), intent(out)        :: val    (KA,IA,JA)             ! value (target)

    logical,  intent(in), optional:: logwgt                  ! use logarithmic weighted interpolation?

    logical :: logwgt_

    real(RP), pointer :: work(:,:,:)

    integer :: k, i, j, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_interp',3)

    logwgt_ = .false.
    if ( present(logwgt) ) then
       logwgt_ = logwgt
    endif

    if ( logwgt_ ) then
       allocate( work(KA_ref,IA_ref,JA_ref) )
       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       do j = 1, JA_ref
       do i = 1, IA_ref
       do k = 1, KA_ref
          work(k,i,j) = log( val_ref(k,i,j) )
       enddo
       enddo
       enddo
    else
       work => val_ref
    endif

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       val(k,i,j) = hfact(i,j,1) * vfact(k,1,i,j,1) * work(idx_k(k,1,i,j,1),idx_i(i,j,1),idx_j(i,j,1)) &
                  + hfact(i,j,1) * vfact(k,2,i,j,1) * work(idx_k(k,2,i,j,1),idx_i(i,j,1),idx_j(i,j,1))
    enddo
    enddo
    enddo

!OCL SERIAL
    do n = 2, npoints
    !$omp parallel do OMP_SCHEDULE_ collapse(2)
!OCL PREFETCH
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       val(k,i,j) = val(k,i,j) &
                  + hfact(i,j,n) * vfact(k,1,i,j,n) * work(idx_k(k,1,i,j,n),idx_i(i,j,n),idx_j(i,j,n)) &
                  + hfact(i,j,n) * vfact(k,2,i,j,n) * work(idx_k(k,2,i,j,n),idx_i(i,j,n),idx_j(i,j,n))
    enddo
    enddo
    enddo
    enddo

    if ( logwgt_ ) then
       deallocate( work )
       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          val(k,i,j) = exp( val(k,i,j) )
       end do
       end do
       end do
    endif

    call PROF_rapend  ('INTERP_interp',3)

    return
  end subroutine INTERP_interp3d

  !-----------------------------------------------------------------------------
  ! search of nearest region for speed up of interpolation
!OCL SERIAL
  subroutine INTERP_search_nearest_block( &
       IA_ref,  &
       JA_ref,  &
       lon_ref, &
       lat_ref, &
       lon,     &
       lat,     &
       IS,      &
       IE,      &
       JS,      &
       JE,      &
       divnum   )
    use scale_const, only: &
       CONST_RADIUS
    implicit none

    integer,  intent(in)  :: IA_ref                 ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                 ! number of y-direction    (reference)
    real(RP), intent(in)  :: lon_ref(IA_ref,JA_ref) ! longitude [rad]          (reference)
    real(RP), intent(in)  :: lat_ref(IA_ref,JA_ref) ! latitude  [rad]          (reference)
    real(RP), intent(in)  :: lon                    ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat                    ! latitude  [rad]          (target)
    integer,  intent(out) :: IS                     ! start index for x-direction
    integer,  intent(out) :: IE                     ! end   index for x-direction
    integer,  intent(out) :: JS                     ! start index for y-direction
    integer,  intent(out) :: JE                     ! end   index for y-direction
    integer,  intent(in), optional :: divnum

    real(RP) :: distance, dist
    integer  :: iskip, jskip
    integer  :: i_bulk, j_bulk
    integer  :: divnum_

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( present(divnum) ) then
       divnum_ = divnum
    else
       divnum_ = INTERP_divnum
    end if

    iskip = max( (IA_ref+1) / divnum_, 1 )
    jskip = max( (JA_ref+1) / divnum_, 1 )

    dist  = large_number

    j = 1 + (jskip/2)
    do while (j <= JA_ref)
       i = 1 + (iskip/2)
       do while (i <= IA_ref)

          distance = haversine( lon, lat, lon_ref(i,j), lat_ref(i,j), CONST_RADIUS )

          if ( distance < dist ) then
             dist   = distance
             i_bulk = i
             j_bulk = j
          endif

          i = i + iskip
       enddo
       j = j + jskip
    enddo

    ! +- 3 is buffer for 12 points
    IS = max( i_bulk - (iskip/2) - 3, 1      )
    IE = min( i_bulk + (iskip/2) + 3, IA_ref )
    JS = max( j_bulk - (jskip/2) - 3, 1      )
    JE = min( j_bulk + (jskip/2) + 3, JA_ref )

    return
  end subroutine INTERP_search_nearest_block

  !-----------------------------------------------------------------------------
  ! horizontal search of interpolation points for three-points
  subroutine INTERP_search_horiz( &
       npoints, &
       IA_ref,  &
       JA_ref,  &
       IS,      &
       IE,      &
       JS,      &
       JE,      &
       lon_ref, &
       lat_ref, &
       lon,     &
       lat,     &
       idx_i,   &
       idx_j,   &
       hfact    )
    use scale_const, only: &
       CONST_RADIUS, &
       CONST_EPS
    use scale_sort, only: &
       SORT_exec
    implicit none

    integer,  intent(in)  :: npoints                ! number of interpolation point for horizontal
    integer,  intent(in)  :: IA_ref                 ! number of x-direction    (reference)
    integer,  intent(in)  :: JA_ref                 ! number of y-direction    (reference)
    integer,  intent(in)  :: IS                     ! start index for x-direction
    integer,  intent(in)  :: IE                     ! end   index for x-direction
    integer,  intent(in)  :: JS                     ! start index for y-direction
    integer,  intent(in)  :: JE                     ! end   index for y-direction
    real(RP), intent(in)  :: lon_ref(IA_ref,JA_ref) ! longitude [rad]          (reference)
    real(RP), intent(in)  :: lat_ref(IA_ref,JA_ref) ! latitude  [rad]          (reference)
    real(RP), intent(in)  :: lon                    ! longitude [rad]          (target)
    real(RP), intent(in)  :: lat                    ! latitude  [rad]          (target)
    integer,  intent(out) :: idx_i(npoints)         ! i-index in reference     (target)
    integer,  intent(out) :: idx_j(npoints)         ! j-index in reference     (target)
    real(RP), intent(out) :: hfact(npoints)         ! horizontal interp factor (target)

    real(RP) :: distance
    real(RP) :: dist(npoints)
    real(RP) :: sum

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    dist (:) = large_number
    idx_i(:) = -1
    idx_j(:) = -1

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(distance)
    do j = JS, JE
    do i = IS, IE
       distance = haversine( lon, lat, lon_ref(i,j), lat_ref(i,j), CONST_RADIUS )

       if ( distance <= dist(npoints) ) then
          ! replace last(=longest) value
          dist (npoints) = distance
          idx_i(npoints) = i
          idx_j(npoints) = j

          ! sort by ascending order
          call SORT_exec( npoints,           & ! [IN]
                          dist (:),          & ! [INOUT]
                          idx_i(:), idx_j(:) ) ! [INOUT]
       endif
    enddo
    enddo

    if ( abs(dist(1)) < CONST_EPS ) then
       hfact(:) = 0.0_RP
       hfact(1) = 1.0_RP
    else
       ! factor = 1 / distance
       do n = 1, npoints
          hfact(n) = 1.0_RP / dist(n)**INTERP_weight_order
       enddo

       ! ignore far point
       do n = 1, npoints
          if ( dist(n) >= INTERP_search_limit ) then
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
    endif

    return
  end subroutine INTERP_search_horiz

  !-----------------------------------------------------------------------------
  ! vertical search of interpolation points for two-points
!OCL SERIAL
  subroutine INTERP_search_vert( &
       KA_ref, KS_ref, KE_ref, &
       KA, KS, KE,             &
       hgt_ref,                &
       hgt,                    &
       idx_k,                  &
       vfact                   )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer,  intent(in)    :: KA_ref, KS_ref, KE_ref        ! number of z-direction    (reference)
    integer,  intent(in)    :: KA, KS, KE                    ! number of z-direction    (target)
    real(RP), intent(in)    :: hgt_ref(KA_ref)               ! height [m]               (reference)
    real(RP), intent(in)    :: hgt    (KA)                   ! height [m]               (target)
    integer,  intent(inout) :: idx_k  (KA,2)                 ! k-index in reference     (target)
    real(RP), intent(inout) :: vfact  (KA,2)                 ! horizontal interp factor (target)

    real(RP) :: weight
    integer  :: k, kk
    !---------------------------------------------------------------------------

    do k = KS, KE
       idx_k(k,1) = -1
       idx_k(k,2) = -1

       if    ( hgt(k) <  hgt_ref(KS_ref) ) then
          idx_k(k,1) = KS_ref
          idx_k(k,2) = KS_ref
          vfact(k,1) = 1.0_RP
          vfact(k,2) = 0.0_RP
       elseif( hgt(k) >= hgt_ref(KE_ref) ) then
          idx_k(k,1) = KE_ref
          idx_k(k,2) = KE_ref
          vfact(k,1) = 1.0_RP
          vfact(k,2) = 0.0_RP
       else
          do kk = KS_ref, KE_ref-1
             if (       hgt(k) >= hgt_ref(kk  ) &
                  .AND. hgt(k) <  hgt_ref(kk+1) ) then

                weight = ( hgt    (k)    - hgt_ref(kk) ) &
                       / ( hgt_ref(kk+1) - hgt_ref(kk) )

                idx_k(k,1) = kk
                idx_k(k,2) = kk + 1
                vfact(k,1) = 1.0_RP - weight
                vfact(k,2) =          weight
                exit
             endif
          enddo
       endif

       if ( idx_k(k,1) < 0 ) then
          LOG_ERROR("INTERP_search_vert",*) 'data for interpolation was not found.'
          LOG_ERROR_CONT(*) 'k=', k, ', hgt(k)=', hgt(k), ', hgt_ref(:)=', hgt_ref(:)
          call PRC_abort
       endif

    enddo ! k-loop

    return
  end subroutine INTERP_search_vert

  !-----------------------------------------------------------------------------
  ! Haversine Formula (from R.W. Sinnott, "Virtues of the Haversine",
  ! Sky and Telescope, vol. 68, no. 2, 1984, p. 159):
  function haversine( &
       lon0,  &
       lat0,  &
       lon1,  &
       lat1,  &
       r_in_m ) &
       result( d )
    implicit none

    real(RP), intent(in) :: lon0   ! [rad]
    real(RP), intent(in) :: lat0   ! [rad]
    real(RP), intent(in) :: lon1   ! [rad]
    real(RP), intent(in) :: lat1   ! [rad]
    real(RP), intent(in) :: r_in_m ! [m]
    real(RP)             :: d      ! [m]

    real(RP) :: dlonh, dlath
    real(RP) :: work1
    !---------------------------------------------------------------------------

    dlonh = 0.5_RP * ( lon0 - lon1 )
    dlath = 0.5_RP * ( lat0 - lat1 )

    work1 = sin(dlath)**2 + cos(lat0) * cos(lat1) * sin(dlonh)**2

    d = r_in_m * 2.0_RP * asin( min(sqrt(work1),1.0_RP) )

  end function haversine

end module scale_interp
