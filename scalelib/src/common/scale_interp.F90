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
  private :: INTERP_search_horiz
  private :: INTERP_search_vert
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
         INTERP_buffer_size_fact

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
  ! make interpolation factor using Lat-Lon
  subroutine INTERP_factor2d( &
       npoints,             &
       IA_ref, JA_ref,      &
       lon_ref,lat_ref,     &
       IA, JA,              &
       lon, lat,            &
       idx_i, idx_j, hfact, &
       search_limit,        &
       latlon_structure,    &
       weight_order         )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
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

    real(RP), intent(in), optional :: search_limit
    logical,  intent(in), optional :: latlon_structure
    integer,  intent(in), optional :: weight_order

    logical :: ll_struct_

    real(RP) :: lon_min, lon_max
    real(RP) :: lat_min, lat_max
    real(RP) :: dlon, dlat

    ! for structure grid
    integer  :: is, ie, js, je
    integer  :: psizex, psizey
    real(RP) :: lon1d(IA_ref), lat1d(JA_ref)
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

       lon1d(:) = UNDEF
       lon_min  = UNDEF
       do i = 1, IA_ref
          do j = 1, JA_ref
             if ( lon_ref(i,j) .ne. UNDEF ) then
                lon1d(i) = lon_ref(i,j)
                exit
             end if
          end do
          if ( lon_min == UNDEF .and. lon1d(i) .ne. UNDEF ) then
             lon_min = lon1d(i)
             is = i
          endif
          if ( lon_min .ne. UNDEF ) then
             if ( lon1d(i) .ne. UNDEF ) then
                lon_max = lon1d(i)
                ie = i
             end if
          end if
       end do

       lat1d(:) = UNDEF
       lat_min  = UNDEF
       do j = 1, JA_ref
          do i = 1, IA_ref
             if ( lat_ref(i,j) .ne. UNDEF ) then
                lat1d(j) = lat_ref(i,j)
                exit
             end if
          end do
          if ( lat_min == UNDEF .and. lat1d(j) .ne. UNDEF ) then
             lat_min = lat1d(j)
             js = j
          endif
          if ( lat_min .ne. UNDEF ) then
             if ( lat1d(j) .ne. UNDEF ) then
                lat_max = lat1d(j)
                je = j
             end if
          end if
       end do

       ! fill undef
       dlon = ( lon_max - lon_min ) / ( ie - is )
       do i = is, ie
          if ( lon1d(i) == UNDEF ) lon1d(i) = lon_min + dlon * ( i - is )
       end do
       dlat = ( lat_max - lat_min ) / ( je - js )
       do j = js, je
          if ( lat1d(j) == UNDEF ) lat1d(j) = lat_min + dlat * ( j - js )
       end do



       if ( ie-is > 10 ) then
          psizex = int( 2.0_RP*sqrt(real(ie-is+1,RP)) )
       else
          psizex = 1
       end if
       if ( je-js > 10 ) then
          psizey = int( 2.0_RP*sqrt(real(je-js+1,RP)) )
       else
          psizey = 1
       end if

       allocate( i0(psizex), i1(psizex) )
       allocate( j0(psizey), j1(psizey) )

       dlon = ( lon_max - lon_min ) / psizex
       dlat = ( lat_max - lat_min ) / psizey

       do ii = 1, psizex
          lon0 = lon_min + dlon * (ii-1)
          do i = is, ie
             if ( lon1d(i) >= lon0 ) then
                i0(ii) = i
                exit
             end if
          end do
       end do
       do ii = 1, psizex-1
          i1(ii) = i0(ii+1) - 1
       end do
       i1(psizex) = ie

       do jj = 1, psizey
          lat0 = lat_min + dlat * (jj-1)
          do j = js, je
             if ( lat1d(j) >= lat0 ) then
                j0(jj) = j
                exit
             end if
          end do
       end do
       do jj = 1, psizey-1
          j1(jj) = j0(jj+1) - 1
       end do
       j1(psizey) = je

       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       do j = 1, JA
       do i = 1, IA
          ! main search
          call INTERP_search_horiz_struct( npoints,                     & ! [IN]
                                           psizex, psizey,              & ! [IN]
                                           IA_ref, JA_ref,              & ! [IN]
                                           lon1d(:), lat1d(:),          & ! [IN]
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
  end subroutine INTERP_factor2d

  !-----------------------------------------------------------------------------
  ! make interpolation factor using Lat-Lon and Z-Height information
  subroutine INTERP_factor3d( &
       npoints, &
       KA_ref, KS_ref, KE_ref, &
       IA_ref, JA_ref,         &
       lon_ref, lat_ref,       &
       hgt_ref,                &
       KA, KS, KE,             &
       IA, JA,                 &
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
    integer,  intent(out) :: idx_k(KA,2,IA,JA,npoints)     ! k-index in reference     (target)
    real(RP), intent(out) :: vfact(KA,2,IA,JA,npoints)     ! vertical interp factor   (target)

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

          call INTERP_search_vert( KA_ref, KS_ref, KE_ref,   & ! [IN]
                                   KA,     KS,     KE,       & ! [IN]
                                   hgt_ref(:,ii,jj),         & ! [IN]
                                   hgt    (:,i,j),           & ! [IN]
                                   idx_k  (:,:,i,j,n),       & ! [OUT]
                                   vfact  (:,:,i,j,n),       & ! [OUT]
                                   flag_extrap = flag_extrap ) ! [IN, optional]

       enddo
    enddo
    enddo

    deallocate(idx_blk, nidx)

    call PROF_rapend  ('INTERP_fact',3)

    return
  end subroutine INTERP_factor3d

  !-----------------------------------------------------------------------------
  ! interpolation using one-points for 2D data (nearest-neighbor)
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
    integer,  intent(in)  :: idx_i  (IA,JA,npoints) ! i-index in reference     (target)
    integer,  intent(in)  :: idx_j  (IA,JA,npoints) ! j-index in reference     (target)
    real(RP), intent(in)  :: hfact  (IA,JA,npoints) ! horizontal interp factor (target)
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
       KA_ref, IA_ref, JA_ref, &
       KA, KS, KE,             &
       IA, JA,                 &
       idx_i, idx_j,           &
       hfact,                  &
       idx_k,                  &
       vfact,                  &
       val_ref,                &
       val,                    &
       logwgt, threshold_undef )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
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

    logical,  intent(in), optional :: logwgt          !> use logarithmic weighted interpolation?
    real(RP), intent(in), optional :: threshold_undef !> return UNDEF if sum of the weight factor is undef the shreshold


    logical  :: logwgt_
    real(RP) :: th_undef

    real(RP), pointer :: work(:,:,:)
    real(RP)          :: valn
    real(RP)          :: fact
    real(RP)          :: w1, w2
    real(RP)          :: f1, f2
    real(RP)          :: sw

    integer :: k, i, j, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('INTERP_interp',3)

    logwgt_ = .false.
    if ( present(logwgt) ) then
       logwgt_ = logwgt
    endif

    th_undef = 0.0_RP
    if ( present(threshold_undef) ) then
       th_undef = threshold_undef
    end if
    th_undef = max( th_undef, EPS * 2.0_RP )

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

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(valn,fact,w1,w2,f1,f2,sw)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       fact = 0.0_RP
       valn = 0.0_RP
       do n = 1, npoints
          w1 = work(idx_k(k,1,i,j,n),idx_i(i,j,n),idx_j(i,j,n))
          w2 = work(idx_k(k,2,i,j,n),idx_i(i,j,n),idx_j(i,j,n))
          f1 = hfact(i,j,n) * vfact(k,1,i,j,n)
          f2 = hfact(i,j,n) * vfact(k,2,i,j,n)
          if ( ( f1 + f2 ) > EPS .and. w1 .ne. UNDEF .and. w2 .ne. UNDEF ) then
             fact = fact + f1 + f2
             valn = valn + f1 * w1 + f2 * w2
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

    if ( logwgt_ ) then
       deallocate( work )
       !$omp parallel do OMP_SCHEDULE_ collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          if ( val(k,i,j) /= UNDEF ) then
             val(k,i,j) = exp( val(k,i,j) )
          endif
       end do
       end do
       end do
    endif

    call PROF_rapend  ('INTERP_interp',3)

    return
  end subroutine INTERP_interp3d

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

  !-----------------------------------------------------------------------------
  ! vertical search of interpolation points for two-points
!OCL SERIAL
  subroutine INTERP_search_vert( &
       KA_ref, KS_ref, KE_ref, &
       KA, KS, KE,             &
       hgt_ref,                &
       hgt,                    &
       idx_k,                  &
       vfact,                  &
       flag_extrap             )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer,  intent(in)  :: KA_ref, KS_ref, KE_ref        ! number of z-direction    (reference)
    integer,  intent(in)  :: KA, KS, KE                    ! number of z-direction    (target)
    real(RP), intent(in)  :: hgt_ref(KA_ref)               ! height [m]               (reference)
    real(RP), intent(in)  :: hgt    (KA)                   ! height [m]               (target)

    integer,  intent(out) :: idx_k  (KA,2)                 ! k-index in reference     (target)
    real(RP), intent(out) :: vfact  (KA,2)                 ! horizontal interp factor (target)

    logical,  intent(in), optional :: flag_extrap          ! when true, extrapolation will be executed (just copy)

    logical :: flag_extrap_

    real(RP) :: weight
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

       if    ( hgt(k) <  hgt_ref(KS_ref) ) then
          if ( flag_extrap_ ) then
             idx_k(k,1) = KS_ref
             idx_k(k,2) = KS_ref ! dummy
             vfact(k,1) = 1.0_RP
             vfact(k,2) = 0.0_RP
          else
             idx_k(k,1) = KS_ref ! dummy
             idx_k(k,2) = KS_ref ! dummy
             vfact(k,1) = 0.0_RP
             vfact(k,2) = 0.0_RP
          end if
       elseif( hgt(k) >= hgt_ref(KE_ref) ) then
          if ( flag_extrap_ ) then
             idx_k(k,1) = KE_ref
             idx_k(k,2) = KE_ref ! dummy
             vfact(k,1) = 1.0_RP
             vfact(k,2) = 0.0_RP
          else
             idx_k(k,1) = KE_ref ! dummy
             idx_k(k,2) = KE_ref ! dummy
             vfact(k,1) = 0.0_RP
             vfact(k,2) = 0.0_RP
          end if
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

  ! private

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

end module scale_interp
