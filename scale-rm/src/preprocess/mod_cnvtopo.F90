!-------------------------------------------------------------------------------
!> module Convert topography
!!
!! @par Description
!!          subroutines for preparing topography data (convert from external file)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cnvtopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CNVTOPO_setup
  public :: CNVTOPO

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: CNVTOPO_DoNothing
  logical, public :: CNVTOPO_UseGTOPO30   = .false.
  logical, public :: CNVTOPO_UseGMTED2010 = .false.
  logical, public :: CNVTOPO_UseDEM50M    = .false.
  logical, public :: CNVTOPO_UseUSERFILE  = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNVTOPO_GTOPO30
  private :: CNVTOPO_GMTED2010
  private :: CNVTOPO_DEM50M
  private :: CNVTOPO_USERFILE
  private :: CNVTOPO_smooth

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: CNVTOPO_smooth_type           = 'LAPLACIAN' ! smoothing type
                                                                   ! 'OFF'         Do not apply smoothing
                                                                   ! 'LAPLACIAN'   Laplacian filter
                                                                   ! 'GAUSSIAN'    Gaussian filter

  integer,                private :: CNVTOPO_smooth_hypdiff_niter  = 20
  logical,                private :: CNVTOPO_smooth_local          = .true.
  integer,                private :: CNVTOPO_smooth_itelim         = 10000

  logical,                private :: CNVTOPO_copy_parent           = .false.

  real(RP),               private :: CNVTOPO_unittile_ddeg         =  0.0_RP ! dx for unit tile [deg]
  real(RP),               private :: CNVTOPO_oversampling_factor   =  2.0_RP ! factor of min. dx against the unit tile
  real(RP),               private :: CNVTOPO_smooth_maxslope_ratio =  1.0_RP ! ratio of DZDX, DZDY
  real(RP),               private :: CNVTOPO_smooth_maxslope       = -1.0_RP ! [deg]

  real(RP),               private :: CNVTOPO_smooth_maxslope_limit

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNVTOPO_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       D2R  => CONST_D2R, &
       HUGE => CONST_HUGE
    use scale_comm_cartesC, only: &
       COMM_horizontal_min
    use scale_atmos_grid_cartesC, only: &
       CDZ => ATMOS_GRID_CARTESC_CDZ, &
       FDX => ATMOS_GRID_CARTESC_FDX, &
       FDY => ATMOS_GRID_CARTESC_FDY
    use scale_atmos_grid_cartesC_real, only: &
       DLAT => ATMOS_GRID_CARTESC_REAL_DLAT, &
       DLON => ATMOS_GRID_CARTESC_REAL_DLON
    implicit none

    character(len=H_SHORT) :: CNVTOPO_name = 'NONE' ! keep backward compatibility

    namelist / PARAM_CNVTOPO / &
       CNVTOPO_name,                  &
       CNVTOPO_UseGTOPO30,            &
       CNVTOPO_UseGMTED2010,          &
       CNVTOPO_UseDEM50M,             &
       CNVTOPO_UseUSERFILE,           &
       CNVTOPO_unittile_ddeg,         &
       CNVTOPO_oversampling_factor,   &
       CNVTOPO_smooth_hypdiff_niter,  &
       CNVTOPO_smooth_maxslope_ratio, &
       CNVTOPO_smooth_maxslope,       &
       CNVTOPO_smooth_local,          &
       CNVTOPO_smooth_itelim,         &
       CNVTOPO_smooth_type,           &
       CNVTOPO_copy_parent

    real(RP) :: minslope(IA,JA)
    real(RP) :: DXL(IA-1)
    real(RP) :: DYL(JA-1)
    real(RP) :: DZDX, DZDY

    real(RP) :: drad(IA,JA)
    real(RP) :: drad_min

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_setup",*) 'Setup'

    DXL(:) = FDX(:)
    DYL(:) = FDY(:)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVTOPO_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVTOPO_setup",*) 'Not appropriate names in namelist PARAM_CNVTOPO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVTOPO)

    select case(CNVTOPO_name)
    case('NONE')
       ! do nothing
    case('GTOPO30')
       CNVTOPO_UseGTOPO30   = .true.
       CNVTOPO_UseGMTED2010 = .false.
       CNVTOPO_UseDEM50M    = .false.
       CNVTOPO_UseUSERFILE  = .false.
    case('GMTED2010')
       CNVTOPO_UseGTOPO30   = .false.
       CNVTOPO_UseGMTED2010 = .true.
       CNVTOPO_UseDEM50M    = .false.
       CNVTOPO_UseUSERFILE  = .false.
    case('DEM50M')
       CNVTOPO_UseGTOPO30   = .false.
       CNVTOPO_UseGMTED2010 = .false.
       CNVTOPO_UseDEM50M    = .true.
       CNVTOPO_UseUSERFILE  = .false.
    case('COMBINE')
       CNVTOPO_UseGTOPO30   = .true.
       CNVTOPO_UseGMTED2010 = .true.
       CNVTOPO_UseDEM50M    = .true.
       CNVTOPO_UseUSERFILE  = .false.
    case('USERFILE')
       ! You can use GTOPO30, GMTED2010, DEM50M and combine User-defined file as you like
       CNVTOPO_UseUSERFILE  = .true.
    case default
       LOG_ERROR("CNVTOPO_setup",*) 'Unsupported TYPE: ', trim(CNVTOPO_name)
       call PRC_abort
    endselect

    CNVTOPO_DoNothing = .true.

    if ( CNVTOPO_UseGTOPO30 ) then
       CNVTOPO_DoNothing = .false.
       LOG_INFO("CNVTOPO_setup",*) 'Use GTOPO, global 30 arcsec. data'
       if ( CNVTOPO_UseGMTED2010 ) then
          LOG_INFO("CNVTOPO_setup",*) 'Use GMTED2010, new global 5 arcsec. data'
          LOG_INFO("CNVTOPO_setup",*) 'Overwrite Existing region'
       endif
       if ( CNVTOPO_UseDEM50M ) then
          LOG_INFO("CNVTOPO_setup",*) 'Use DEM 50m data for Japan region'
          LOG_INFO("CNVTOPO_setup",*) 'Overwrite Japan region'
       endif
    elseif ( CNVTOPO_UseGMTED2010 ) then
       CNVTOPO_DoNothing = .false.
       LOG_INFO("CNVTOPO_setup",*) 'Use GMTED2010, new global 5 arcsec. data'
       if ( CNVTOPO_UseDEM50M ) then
          LOG_INFO("CNVTOPO_setup",*) 'Use DEM 50m data for Japan region'
          LOG_INFO("CNVTOPO_setup",*) 'Overwrite Japan region'
       endif
    elseif ( CNVTOPO_UseDEM50M ) then
       CNVTOPO_DoNothing = .false.
       LOG_INFO("CNVTOPO_setup",*) 'Use DEM 50m data, Japan region only'
    elseif ( CNVTOPO_UseUSERFILE ) then
       CNVTOPO_DoNothing = .false.
       LOG_INFO("CNVTOPO_setup",*) 'Use user-defined file'
    endif

    if ( CNVTOPO_DoNothing ) then
       LOG_INFO("CNVTOPO_setup",*) 'Do nothing for topography data'
    else
       drad(:,:) = min( DLAT(:,:), DLON(:,:) )
       call COMM_horizontal_min( drad_min, drad(:,:) )

       if ( CNVTOPO_unittile_ddeg > 0.0_RP ) then
          CNVTOPO_oversampling_factor = ( drad_min / D2R ) / CNVTOPO_unittile_ddeg
       endif
       CNVTOPO_oversampling_factor = max( 1.0_RP, CNVTOPO_oversampling_factor )
       CNVTOPO_unittile_ddeg       = ( drad_min / D2R ) / CNVTOPO_oversampling_factor

       LOG_INFO("CNVTOPO_setup",*) 'The size of tile [deg] = ', CNVTOPO_unittile_ddeg
       LOG_INFO("CNVTOPO_setup",*) 'oversampling factor    = ', CNVTOPO_oversampling_factor
    endif

    if( CNVTOPO_smooth_maxslope > 0.0_RP ) then

      CNVTOPO_smooth_maxslope_limit = CNVTOPO_smooth_maxslope

    else
      minslope(:,:) = HUGE

      j = JS-1
      i = IS-1
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DYL(j) ) / D2R
         minslope(IS,JS) = min( minslope(IS,JS), DZDX, DZDY )
      enddo

      j = JS-1
      do i = IS, IE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DYL(j) ) / D2R
         minslope(i,JS) = min( minslope(i,JS), DZDX, DZDY )
      enddo
      enddo

      i = IS-1
      do j = JS, JE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DYL(j) ) / D2R
         minslope(IS,j) = min( minslope(IS,j), DZDX, DZDY )
      enddo
      enddo

      do j = JS, JE
      do i = IS, IE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DYL(j) ) / D2R
         minslope(i,j) = min( minslope(i,j), DZDX, DZDY )
      enddo
      enddo
      enddo

      call COMM_horizontal_min( CNVTOPO_smooth_maxslope_limit, minslope(:,:) )
    end if

    return
  end subroutine CNVTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVTOPO
    use scale_prc, only: &
       PRC_abort
    use scale_topography, only: &
       TOPO_fillhalo, &
       TOPO_Zsfc, &
       TOPO_write
    use mod_copytopo, only: &
       COPYTOPO
    implicit none
    !---------------------------------------------------------------------------

    if ( CNVTOPO_DoNothing ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'skip  convert topography data'
    else
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start convert topography data'

       if ( CNVTOPO_UseGTOPO30 ) then
          call CNVTOPO_GTOPO30
       endif

       if ( CNVTOPO_UseGMTED2010 ) then
          call CNVTOPO_GMTED2010
       endif

       if ( CNVTOPO_UseDEM50M ) then
          call CNVTOPO_DEM50M
       endif

       if ( CNVTOPO_UseUSERFILE ) then
          call CNVTOPO_USERFILE
       endif

       call CNVTOPO_smooth( TOPO_Zsfc(:,:) ) ! (inout)
       call TOPO_fillhalo( FILL_BND=.true. )

       if( CNVTOPO_copy_parent ) call COPYTOPO( TOPO_Zsfc )

       LOG_PROGRESS(*) 'end   convert topography data'

       ! output topography file
       call TOPO_write
    endif

    return
  end subroutine CNVTOPO

  !-----------------------------------------------------------------------------
  !> Convert from GTOPO30
  subroutine CNVTOPO_GTOPO30
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_atmos_grid_cartesC_real, only: &
       LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       LONUY => ATMOS_GRID_CARTESC_REAL_LONUY
    implicit none

    character(len=H_LONG) :: GTOPO30_IN_DIR       = '.' !< directory contains GTOPO30 files (GrADS format)
    character(len=H_LONG) :: GTOPO30_IN_CATALOGUE = ''  !< metadata files for GTOPO30

    namelist / PARAM_CNVTOPO_GTOPO30 / &
       GTOPO30_IN_DIR,       &
       GTOPO30_IN_CATALOGUE

    ! data catalogue list
    integer, parameter      :: TILE_nlim = 100
    integer                 :: TILE_nmax
    real(RP)                :: TILE_LATS (TILE_nlim)
    real(RP)                :: TILE_LATE (TILE_nlim)
    real(RP)                :: TILE_LONS (TILE_nlim)
    real(RP)                :: TILE_LONE (TILE_nlim)
    character(len=H_LONG)   :: TILE_fname(TILE_nlim)

    ! GTOPO30 data
    integer,  parameter     :: jsize_orig   = 6000
    integer,  parameter     :: isize_orig   = 4800
    real(RP), parameter     :: GTOPO30_DLAT = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.
    real(RP), parameter     :: GTOPO30_DLON = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.
    real(RP)                :: TILE_DLAT_orig, TILE_DLON_orig
    integer(2)              :: TILE_HEIGHT_orig(isize_orig,jsize_orig)

    ! GTOPO30 data (oversampling)
    integer                 :: jos, ios
    integer                 :: jsize, isize
    real(RP)                :: TILE_DLAT, TILE_DLON
    integer(2), allocatable :: TILE_HEIGHT(:,:)
    real(RP),   allocatable :: TILE_LATH  (:)
    real(RP),   allocatable :: TILE_LONH  (:)
    real(RP)                :: area, area_fraction

    integer  :: jloc, iloc
    real(RP) :: jfrac_b ! fraction for jloc
    real(RP) :: ifrac_l ! fraction for iloc

    real(RP) :: LONUY_mod(0:IA,JA)
    real(RP) :: DOMAIN_LATS, DOMAIN_LATE
    real(RP) :: DOMAIN_LONS, DOMAIN_LONE
    integer  :: DOMAIN_LONSLOC(2), DOMAIN_LONELOC(2)
    logical  :: check_IDL

    real(RP) :: topo_sum(IA,JA)
    real(RP) :: area_sum(IA,JA)
    real(RP) :: topo, mask

    character(len=H_LONG) :: fname

    real(RP) :: zerosw
    logical  :: hit_lat, hit_lon
    integer  :: index
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, iii, jjj, t
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_GTOPO30",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO_GTOPO30,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVTOPO_GTOPO30",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVTOPO_GTOPO30",*) 'Not appropriate names in namelist PARAM_CNVTOPO_GTOPO30. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVTOPO_GTOPO30)

    do j = 1, JA
    do i = 1, IA
       area_sum(i,j) = 0.0_RP
       topo_sum(i,j) = 0.0_RP
    enddo
    enddo

    LONUY_mod(:,:) = mod( LONUY(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI

    DOMAIN_LATS    = minval(LATXV    (:,:))
    DOMAIN_LATE    = maxval(LATXV    (:,:))
    DOMAIN_LONS    = minval(LONUY_mod(:,:))
    DOMAIN_LONE    = maxval(LONUY_mod(:,:))
    DOMAIN_LONSLOC = minloc(LONUY_mod(:,:))
    DOMAIN_LONELOC = maxloc(LONUY_mod(:,:))

    check_IDL = .false.
    if (      DOMAIN_LONS < LONUY_mod( 0,DOMAIN_LONSLOC(2)) &
         .OR. DOMAIN_LONE > LONUY_mod(IA,DOMAIN_LONELOC(2)) ) then
       check_IDL = .true.
       DOMAIN_LONS = minval(LONUY_mod(:,:),mask=(LONUY_mod>0.0_RP))
       DOMAIN_LONE = maxval(LONUY_mod(:,:),mask=(LONUY_mod<0.0_RP))
    endif

    jos   = nint( GTOPO30_DLAT / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    ios   = nint( GTOPO30_DLON / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    jsize = jsize_orig * jos
    isize = isize_orig * ios

    allocate( TILE_HEIGHT(isize,jsize) )
    allocate( TILE_LATH  (0:jsize)     )
    allocate( TILE_LONH  (0:isize)     )

    LOG_INFO_CONT(*) 'Oversampling (j) orig = ', jsize_orig, ', use = ', jsize
    LOG_INFO_CONT(*) 'Oversampling (i) orig = ', isize_orig, ', use = ', isize

    TILE_DLAT_orig = GTOPO30_DLAT * D2R
    TILE_DLON_orig = GTOPO30_DLON * D2R
    LOG_INFO_CONT(*) 'TILE_DLAT       :', TILE_DLAT_orig/D2R
    LOG_INFO_CONT(*) 'TILE_DLON       :', TILE_DLON_orig/D2R

    TILE_DLAT = TILE_DLAT_orig / jos
    TILE_DLON = TILE_DLON_orig / ios
    LOG_INFO_CONT(*) 'TILE_DLAT (OS)  :', TILE_DLAT/D2R
    LOG_INFO_CONT(*) 'TILE_DLON (OS)  :', TILE_DLON/D2R

    !---< READ from external files >---

    ! catalogue file
    fname = trim(GTOPO30_IN_DIR)//'/'//trim(GTOPO30_IN_CATALOGUE)

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_GTOPO30",*) 'Input catalogue file:', trim(fname)

    fid = IO_get_available_fid()
    open( fid,                  &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          LOG_ERROR("CNVTOPO_GTOPO30",*) 'catalogue file not found! ', trim(fname)
          call PRC_abort
       endif

       do t = 1, TILE_nlim
          read(fid,*,iostat=ierr) index, TILE_LATS(t), TILE_LATE(t), & ! South->North
                                         TILE_LONS(t), TILE_LONE(t), & ! WEST->EAST
                                         TILE_fname(t)
          if ( ierr /= 0 ) exit

          if ( TILE_LONS(t) >= 180.0_RP ) then
             TILE_LONS(t) = TILE_LONS(t) - 360.0_RP
             TILE_LONE(t) = TILE_LONE(t) - 360.0_RP
          endif
          if ( TILE_LONS(t) < -180.0_RP ) TILE_LONS(t) = TILE_LONS(t) + 360.0_RP
          if ( TILE_LONE(t) < -180.0_RP ) TILE_LONE(t) = TILE_LONE(t) + 360.0_RP

       enddo

       TILE_nmax = t - 1
    close(fid)

    ! data file
    do t = 1, TILE_nmax
       hit_lat = .false.
       hit_lon = .false.

       if (      ( TILE_LATS(t)*D2R >= DOMAIN_LATS  .AND. TILE_LATS(t)*D2R < DOMAIN_LATE  ) &
            .OR. ( TILE_LATE(t)*D2R >= DOMAIN_LATS  .AND. TILE_LATE(t)*D2R < DOMAIN_LATE  ) ) then
          hit_lat = .true.
       endif

       if (      ( DOMAIN_LATS  >= TILE_LATS(t)*D2R .AND. DOMAIN_LATS  < TILE_LATE(t)*D2R ) &
            .OR. ( DOMAIN_LATE  >= TILE_LATS(t)*D2R .AND. DOMAIN_LATE  < TILE_LATE(t)*D2R ) ) then
          hit_lat = .true.
       endif

       if ( check_IDL ) then
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < PI           ) &
               .OR. ( TILE_LONS(t)*D2R >= -PI          .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < PI           ) &
               .OR. ( TILE_LONE(t)*D2R >= -PI          .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
             hit_lon = .true.
          endif
       else
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
             hit_lon = .true.
          endif
       endif

       if (      ( DOMAIN_LONS  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS  < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE  < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(GTOPO30_IN_DIR)//'/'//trim(TILE_fname(t))

          LOG_NEWLINE
          LOG_INFO("CNVTOPO_GTOPO30",*) 'Input data file :', trim(fname)
          LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if ( check_IDL ) then
             LOG_INFO_CONT(*) '(Date line exists within the domain)'
          endif
          LOG_INFO_CONT(*) 'Tile   (LAT)    :', TILE_LATS(t), TILE_LATE(t)
          LOG_INFO_CONT(*) '       (LON)    :', TILE_LONS(t), TILE_LONE(t)

          fid = IO_get_available_fid()
          open( fid,                              &
                file   = trim(fname),             &
                form   = 'unformatted',           &
                access = 'direct',                &
                status = 'old',                   &
                recl   = isize_orig*jsize_orig*2, &
                iostat = ierr                     )

             if ( ierr /= 0 ) then
                LOG_ERROR("CNVTOPO_GTOPO30",*) 'data file not found!'
                call PRC_abort
             endif

             read(fid,rec=1) TILE_HEIGHT_orig(:,:)
          close(fid)

          ! oversampling
          do jj = 1, jsize_orig
          do ii = 1, isize_orig
             do j = 1, jos
             do i = 1, ios
                jjj = (jj-1) * jos + j
                iii = (ii-1) * ios + i

                TILE_HEIGHT(iii,jjj) = TILE_HEIGHT_orig(ii,jsize_orig-jj+1) ! reverse y-axis
             enddo
             enddo
          enddo
          enddo

          TILE_LATH(0) = TILE_LATS(t) * D2R
          do jj = 1, jsize
             TILE_LATH(jj) = TILE_LATH(jj-1) + TILE_DLAT
             !LOG_INFO("CNVTOPO_GTOPO30",*) jj, TILE_LATH(jj) / D2R
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, isize
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
             !LOG_INFO("CNVTOPO_GTOPO30",*) ii, TILE_LONH(ii) / D2R
          enddo

          ! match and calc fraction
          do jj = 1, jsize
          do ii = 1, isize

             iloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             ifrac_l = 1.0_RP

             jloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             jfrac_b = 1.0_RP

             if (      TILE_LATH(jj  ) < DOMAIN_LATS &
                  .OR. TILE_LATH(jj-1) > DOMAIN_LATE ) then
                cycle
             endif

             if ( check_IDL ) then
                if (       TILE_LONH(ii  ) < DOMAIN_LONS &
                     .AND. TILE_LONH(ii-1) > DOMAIN_LONE ) then
                   cycle
                endif
             else
                if (       TILE_LONH(ii  ) < DOMAIN_LONS &
                     .OR.  TILE_LONH(ii-1) > DOMAIN_LONE ) then
                   cycle
                endif
             endif

      jloop: do j = JS-1, JE+1
      iloop: do i = IS-1, IE+1
                if (       TILE_LONH(ii-1) >= LONUY_mod(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= LATXV    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  LATXV    (i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif

                if (       LONUY_mod(i-1,j) >= LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1)  >= LATXV    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1)  <  LATXV    (i  ,j  ) ) then ! across the IDL

                   if    (       TILE_LONH(ii-1) >= LONUY_mod(i-1,j) &
                           .AND. TILE_LONH(ii-1) <  PI               ) then

                      iloc    = i
                      ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1)+2.0_RP*PI, TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   elseif(       TILE_LONH(ii-1) >= -PI                  &
                           .AND. TILE_LONH(ii-1) <  LONUY_mod(i  ,j) ) then

                      iloc    = i
                      ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   endif

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             topo = real( TILE_HEIGHT(ii,jj), kind=RP )
             mask = 0.5_RP - sign( 0.5_RP, topo ) ! if Height is negative, mask = 1

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) ) * ( 1.0_RP - mask )

!             LOG_INFO("CNVTOPO_GTOPO30",*) ii, jj, area, iloc, jloc, ifrac_l, jfrac_b, TILE_HEIGHT(ii,jj)

             area_fraction = (       ifrac_l) * (       jfrac_b) * area
             area_sum(iloc  ,jloc  ) = area_sum(iloc  ,jloc  ) + area_fraction
             topo_sum(iloc  ,jloc  ) = topo_sum(iloc  ,jloc  ) + area_fraction * topo

             area_fraction = (1.0_RP-ifrac_l) * (       jfrac_b) * area
             area_sum(iloc+1,jloc  ) = area_sum(iloc+1,jloc  ) + area_fraction
             topo_sum(iloc+1,jloc  ) = topo_sum(iloc+1,jloc  ) + area_fraction * topo

             area_fraction = (       ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc  ,jloc+1) = area_sum(iloc  ,jloc+1) + area_fraction
             topo_sum(iloc  ,jloc+1) = topo_sum(iloc  ,jloc+1) + area_fraction * topo

             area_fraction = (1.0_RP-ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc+1,jloc+1) = area_sum(iloc+1,jloc+1) + area_fraction
             topo_sum(iloc+1,jloc+1) = topo_sum(iloc+1,jloc+1) + area_fraction * topo
          enddo
          enddo

       endif
    enddo ! tile loop

    do j = JS, JE
    do i = IS, IE
       zerosw = 0.5_RP - sign( 0.5_RP, area_sum(i,j)-EPS )
       TOPO_Zsfc(i,j) = topo_sum(i,j) * ( 1.0_RP-zerosw ) / ( area_sum(i,j)-zerosw )
    enddo
    enddo

    return
  end subroutine CNVTOPO_GTOPO30

  !-----------------------------------------------------------------------------
  !> Convert from GMTED 2010 7.5 arc-seconds mesh
  subroutine CNVTOPO_GMTED2010
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine CNVTOPO_GMTED2010

  !-----------------------------------------------------------------------------
  !> Convert from DEM 50m mesh
  subroutine CNVTOPO_DEM50M
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_atmos_grid_cartesC_real, only: &
       LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       LONUY => ATMOS_GRID_CARTESC_REAL_LONUY
    implicit none

    character(len=H_LONG) :: DEM50M_IN_DIR       = '.' !< directory contains DEM50M files (GrADS format)
    character(len=H_LONG) :: DEM50M_IN_CATALOGUE = ''  !< metadata files for DEM50M

    namelist / PARAM_CNVTOPO_DEM50M / &
       DEM50M_IN_DIR,      &
       DEM50M_IN_CATALOGUE

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 1000
    integer               :: TILE_nmax
    real(RP)              :: TILE_LATS (TILE_nlim)
    real(RP)              :: TILE_LATE (TILE_nlim)
    real(RP)              :: TILE_LONS (TILE_nlim)
    real(RP)              :: TILE_LONE (TILE_nlim)
    character(len=H_LONG) :: TILE_fname(TILE_nlim)

    ! DEM50M data
    integer,  parameter   :: jsize_orig = 1600
    integer,  parameter   :: isize_orig = 1600
    real(RP), parameter   :: DEM50M_DLAT = 5.0_RP / 60.0_RP / 200.0_RP ! 30 arc sec.
    real(RP), parameter   :: DEM50M_DLON = 7.5_RP / 60.0_RP / 200.0_RP ! 30 arc sec.
    real(SP)              :: TILE_HEIGHT_orig(isize_orig,jsize_orig)
    real(RP)              :: TILE_DLAT_orig, TILE_DLON_orig

    ! DEM50M data (oversampling)
    integer               :: jos, ios
    integer               :: jsize, isize
    real(RP)              :: TILE_DLAT, TILE_DLON
    real(SP), allocatable :: TILE_HEIGHT(:,:)
    real(RP), allocatable :: TILE_LATH  (:)
    real(RP), allocatable :: TILE_LONH  (:)
    real(RP)              :: area, area_fraction

    integer  :: jloc, iloc
    real(RP) :: jfrac_b ! fraction for jloc
    real(RP) :: ifrac_l ! fraction for iloc

    real(RP) :: LONUY_mod(0:IA,JA)
    real(RP) :: DOMAIN_LATS, DOMAIN_LATE
    real(RP) :: DOMAIN_LONS, DOMAIN_LONE
    integer  :: DOMAIN_LONSLOC(2), DOMAIN_LONELOC(2)
    logical  :: check_IDL

    real(RP) :: topo_sum(IA,JA)
    real(RP) :: area_sum(IA,JA)
    real(RP) :: topo, mask

    character(len=H_LONG) :: fname

    real(RP) :: zerosw
    logical  :: hit_lat, hit_lon
    integer  :: index
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, iii, jjj, t
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_DEM50M",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO_DEM50M,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVTOPO_DEM50M",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVTOPO_DEM50M",*) 'Not appropriate names in namelist PARAM_CNVTOPO_DEM50M. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVTOPO_DEM50M)

    do j = 1, JA
    do i = 1, IA
       area_sum(i,j) = 0.0_RP
       topo_sum(i,j) = 0.0_RP
    enddo
    enddo

    LONUY_mod(:,:) = mod( LONUY(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI

    DOMAIN_LATS    = minval(LATXV    (:,:))
    DOMAIN_LATE    = maxval(LATXV    (:,:))
    DOMAIN_LONS    = minval(LONUY_mod(:,:))
    DOMAIN_LONE    = maxval(LONUY_mod(:,:))
    DOMAIN_LONSLOC = minloc(LONUY_mod(:,:))
    DOMAIN_LONELOC = maxloc(LONUY_mod(:,:))

    check_IDL = .false.
    if (      DOMAIN_LONS < LONUY_mod( 0,DOMAIN_LONSLOC(2)) &
         .OR. DOMAIN_LONE > LONUY_mod(IA,DOMAIN_LONELOC(2)) ) then
       check_IDL = .true.
       DOMAIN_LONS = minval(LONUY_mod(:,:),mask=(LONUY_mod>0.0_RP))
       DOMAIN_LONE = maxval(LONUY_mod(:,:),mask=(LONUY_mod<0.0_RP))
    endif

    jos   = nint( DEM50M_DLAT / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    ios   = nint( DEM50M_DLON / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    jsize = jsize_orig * jos
    isize = isize_orig * ios

    allocate( TILE_HEIGHT(isize,jsize) )
    allocate( TILE_LATH  (0:jsize)     )
    allocate( TILE_LONH  (0:isize)     )

    LOG_INFO_CONT(*) 'Oversampling (j) orig = ', jsize_orig, ', use = ', jsize
    LOG_INFO_CONT(*) 'Oversampling (i) orig = ', isize_orig, ', use = ', isize

    TILE_DLAT_orig = DEM50M_DLAT * D2R
    TILE_DLON_orig = DEM50M_DLON * D2R
    LOG_INFO_CONT(*) 'TILE_DLAT       :', TILE_DLAT_orig/D2R
    LOG_INFO_CONT(*) 'TILE_DLON       :', TILE_DLON_orig/D2R

    TILE_DLAT = TILE_DLAT_orig / jos
    TILE_DLON = TILE_DLON_orig / ios
    LOG_INFO_CONT(*) 'TILE_DLAT (OS)  :', TILE_DLAT/D2R
    LOG_INFO_CONT(*) 'TILE_DLON (OS)  :', TILE_DLON/D2R

    !---< READ from external files >---

    ! catalogue file
    fname = trim(DEM50M_IN_DIR)//'/'//trim(DEM50M_IN_CATALOGUE)

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_DEM50M",*) 'Input catalogue file:', trim(fname)

    fid = IO_get_available_fid()
    open( fid,                  &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          LOG_ERROR("CNVTOPO_DEM50M",*) 'catalogue file not found! ', trim(fname)
          call PRC_abort
       endif

       do t = 1, TILE_nlim
          read(fid,*,iostat=ierr) index, TILE_LATS(t), TILE_LATE(t), & ! South->North
                                         TILE_LONS(t), TILE_LONE(t), & ! WEST->EAST
                                         TILE_fname(t)
          if ( ierr /= 0 ) exit

          if ( TILE_LONS(t) >= 180.0_RP ) then
             TILE_LONS(t) = TILE_LONS(t) - 360.0_RP
             TILE_LONE(t) = TILE_LONE(t) - 360.0_RP
          endif
          if ( TILE_LONS(t) < -180.0_RP ) TILE_LONS(t) = TILE_LONS(t) + 360.0_RP
          if ( TILE_LONE(t) < -180.0_RP ) TILE_LONE(t) = TILE_LONE(t) + 360.0_RP

       enddo

       TILE_nmax = t - 1
    close(fid)

    ! data file
    do t = 1, TILE_nmax
       hit_lat = .false.
       hit_lon = .false.

       if (      ( TILE_LATS(t)*D2R >= DOMAIN_LATS  .AND. TILE_LATS(t)*D2R < DOMAIN_LATE  ) &
            .OR. ( TILE_LATE(t)*D2R >= DOMAIN_LATS  .AND. TILE_LATE(t)*D2R < DOMAIN_LATE  ) ) then
          hit_lat = .true.
       endif

       if (      ( DOMAIN_LATS  >= TILE_LATS(t)*D2R .AND. DOMAIN_LATS  < TILE_LATE(t)*D2R ) &
            .OR. ( DOMAIN_LATE  >= TILE_LATS(t)*D2R .AND. DOMAIN_LATE  < TILE_LATE(t)*D2R ) ) then
          hit_lat = .true.
       endif

       if ( check_IDL ) then
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < PI           ) &
               .OR. ( TILE_LONS(t)*D2R >= -PI          .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < PI           ) &
               .OR. ( TILE_LONE(t)*D2R >= -PI          .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
             hit_lon = .true.
          endif
       else
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
             hit_lon = .true.
          endif
       endif

       if (      ( DOMAIN_LONS  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS  < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE  < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(DEM50M_IN_DIR)//'/'//trim(TILE_fname(t))

          LOG_NEWLINE
          LOG_INFO("CNVTOPO_DEM50M",*) 'Input data file :', trim(fname)
          LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if ( check_IDL ) then
             LOG_INFO_CONT(*) '(Date line exists within the domain)'
          endif
          LOG_INFO_CONT(*) 'Tile   (LAT)    :', TILE_LATS(t), TILE_LATE(t)
          LOG_INFO_CONT(*) '       (LON)    :', TILE_LONS(t), TILE_LONE(t)

          fid = IO_get_available_fid()
          open( fid,                              &
                file   = trim(fname),             &
                form   = 'unformatted',           &
                access = 'direct',                &
                status = 'old',                   &
                recl   = isize_orig*jsize_orig*4, &
                iostat = ierr                     )

             if ( ierr /= 0 ) then
                LOG_ERROR("CNVTOPO_DEM50M",*) 'data file not found!'
                call PRC_abort
             endif

             read(fid,rec=1) TILE_HEIGHT_orig(:,:)
          close(fid)

          ! oversampling
          do jj = 1, jsize_orig
          do ii = 1, isize_orig
             do j = 1, jos
             do i = 1, ios
                jjj = (jj-1) * jos + j
                iii = (ii-1) * ios + i

                TILE_HEIGHT(iii,jjj) = TILE_HEIGHT_orig(ii,jj)
             enddo
             enddo
          enddo
          enddo

          TILE_LATH(0) = TILE_LATS(t) * D2R
          do jj = 1, jsize
             TILE_LATH(jj) = TILE_LATH(jj-1) + TILE_DLAT
             !LOG_INFO("CNVTOPO_DEM50M",*) jj, TILE_LATH(jj) / D2R
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, isize
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
             !LOG_INFO("CNVTOPO_DEM50M",*) ii, TILE_LONH(ii) / D2R
          enddo

          ! match and calc fraction
          do jj = 1, jsize
          do ii = 1, isize

             iloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             ifrac_l = 1.0_RP

             jloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             jfrac_b = 1.0_RP

             if (      TILE_LATH(jj  ) < DOMAIN_LATS &
                  .OR. TILE_LATH(jj-1) > DOMAIN_LATE ) then
                cycle
             endif

             if ( check_IDL ) then
                if (       TILE_LONH(ii  ) < DOMAIN_LONS &
                     .AND. TILE_LONH(ii-1) > DOMAIN_LONE ) then
                   cycle
                endif
             else
                if (       TILE_LONH(ii  ) < DOMAIN_LONS &
                     .OR.  TILE_LONH(ii-1) > DOMAIN_LONE ) then
                   cycle
                endif
             endif

      jloop: do j = JS-1, JE+1
      iloop: do i = IS-1, IE+1
                if (       TILE_LONH(ii-1) >= LONUY_mod(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= LATXV    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  LATXV    (i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif

                if (       LONUY_mod(i-1,j) >= LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1)  >= LATXV    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1)  <  LATXV    (i  ,j  ) ) then ! across the IDL

                   if    (       TILE_LONH(ii-1) >= LONUY_mod(i-1,j) &
                           .AND. TILE_LONH(ii-1) <  PI               ) then

                      iloc    = i
                      ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1)+2.0_RP*PI, TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   elseif(       TILE_LONH(ii-1) >= -PI             &
                           .AND. TILE_LONH(ii-1) <  LONUY_mod(i  ,j) ) then

                      iloc    = i
                      ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   endif

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             topo = real( TILE_HEIGHT(ii,jj), kind=RP )
             mask = 0.5_RP - sign( 0.5_RP,topo ) ! if Height is negative, mask = 1

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) ) * ( 1.0_RP - mask )

!             LOG_INFO("CNVTOPO_DEM50M",*) ii, jj, area, iloc, jloc, ifrac_l, jfrac_b, TILE_HEIGHT(ii,jj)

             area_fraction = (       ifrac_l) * (       jfrac_b) * area
             area_sum(iloc  ,jloc  ) = area_sum(iloc  ,jloc  ) + area_fraction
             topo_sum(iloc  ,jloc  ) = topo_sum(iloc  ,jloc  ) + area_fraction * topo

             area_fraction = (1.0_RP-ifrac_l) * (       jfrac_b) * area
             area_sum(iloc+1,jloc  ) = area_sum(iloc+1,jloc  ) + area_fraction
             topo_sum(iloc+1,jloc  ) = topo_sum(iloc+1,jloc  ) + area_fraction * topo

             area_fraction = (       ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc  ,jloc+1) = area_sum(iloc  ,jloc+1) + area_fraction
             topo_sum(iloc  ,jloc+1) = topo_sum(iloc  ,jloc+1) + area_fraction * topo

             area_fraction = (1.0_RP-ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc+1,jloc+1) = area_sum(iloc+1,jloc+1) + area_fraction
             topo_sum(iloc+1,jloc+1) = topo_sum(iloc+1,jloc+1) + area_fraction * topo
          enddo
          enddo

       endif
    enddo ! tile loop

    do j = JS, JE
    do i = IS, IE
       mask   = 0.5_RP + sign( 0.5_RP, area_sum(i,j)-EPS ) ! if any data is found, overwrite
       zerosw = 0.5_RP - sign( 0.5_RP, area_sum(i,j)-EPS )
       topo   = topo_sum(i,j) * ( 1.0_RP-zerosw ) / ( area_sum(i,j)-zerosw )
       TOPO_Zsfc(i,j) = (        mask ) * topo &         ! overwrite
                      + ( 1.0_RP-mask ) * TOPO_Zsfc(i,j) ! keep existing value
    enddo
    enddo

    return
  end subroutine CNVTOPO_DEM50M

  !-----------------------------------------------------------------------------
  !> Convert from User-defined file
  subroutine CNVTOPO_USERFILE
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_atmos_grid_cartesC_real, only: &
       LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       LONUY => ATMOS_GRID_CARTESC_REAL_LONUY
    implicit none

    integer               :: USERFILE_NLAT         = -1       ! number of latitude  tile
    integer               :: USERFILE_NLON         = -1       ! number of longitude tile
    real(RP)              :: USERFILE_DLAT         = -1.0_RP  ! width  of latitude  tile [deg.]
    real(RP)              :: USERFILE_DLON         = -1.0_RP  ! width  of longitude tile [deg.]
    character(len=H_LONG) :: USERFILE_IN_DIR       = '.'      ! directory contains data files (GrADS format)
    character(len=H_LONG) :: USERFILE_IN_CATALOGUE = ''       ! catalogue file
    character(len=H_LONG) :: USERFILE_IN_FILENAME  = ''       ! single data file (GrADS format)
    character(len=H_LONG) :: USERFILE_IN_DATATYPE  = 'REAL4'  ! datatype (REAL4,REAL8,INT2)
    logical               :: USERFILE_LATORDER_N2S = .false.  ! data of the latitude direction is stored in ordar of North->South?
    real(RP)              :: USERFILE_LAT_START    = -90.0_RP ! (for single file) start latitude  of domain in input data
    real(RP)              :: USERFILE_LAT_END      =  90.0_RP ! (for single file) end   latitude  of domain in input data
    real(RP)              :: USERFILE_LON_START    =   0.0_RP ! (for single file) start longitude of domain in input data
    real(RP)              :: USERFILE_LON_END      = 360.0_RP ! (for single file) end   longitude of domain in input data

    namelist / PARAM_CNVTOPO_USERFILE / &
       USERFILE_NLAT,         &
       USERFILE_NLON,         &
       USERFILE_DLAT,         &
       USERFILE_DLON,         &
       USERFILE_IN_CATALOGUE, &
       USERFILE_IN_DIR,       &
       USERFILE_IN_FILENAME,  &
       USERFILE_IN_DATATYPE,  &
       USERFILE_LATORDER_N2S, &
       USERFILE_LAT_START,    &
       USERFILE_LAT_END,      &
       USERFILE_LON_START,    &
       USERFILE_LON_END

    ! data catalogue list
    integer, parameter      :: TILE_nlim = 1000
    integer                 :: TILE_nmax
    real(RP)                :: TILE_LATS (TILE_nlim)
    real(RP)                :: TILE_LATE (TILE_nlim)
    real(RP)                :: TILE_LONS (TILE_nlim)
    real(RP)                :: TILE_LONE (TILE_nlim)
    character(len=H_LONG)   :: TILE_fname(TILE_nlim)

    ! User-defined data
    real(RP),   allocatable :: TILE_HEIGHT_orig   (:,:)
    real(DP),   allocatable :: TILE_HEIGHT_orig_DP(:,:)
    real(SP),   allocatable :: TILE_HEIGHT_orig_SP(:,:)
    integer(2), allocatable :: TILE_HEIGHT_orig_I2(:,:)
    real(RP)                :: TILE_DLAT_orig, TILE_DLON_orig

    ! User-defined data (oversampling)
    integer                 :: jos, ios
    integer                 :: jsize, isize
    real(RP)                :: TILE_DLAT, TILE_DLON
    real(RP), allocatable   :: TILE_HEIGHT(:,:)
    real(RP), allocatable   :: TILE_LATH  (:)
    real(RP), allocatable   :: TILE_LONH  (:)
    real(RP)                :: area, area_fraction

    integer  :: jloc, iloc
    real(RP) :: jfrac_b ! fraction for jloc
    real(RP) :: ifrac_l ! fraction for iloc

    real(RP) :: LONUY_mod(0:IA,JA)
    real(RP) :: DOMAIN_LATS, DOMAIN_LATE
    real(RP) :: DOMAIN_LONS, DOMAIN_LONE
    integer  :: DOMAIN_LONSLOC(2), DOMAIN_LONELOC(2)
    logical  :: check_IDL

    real(RP) :: topo_sum(IA,JA)
    real(RP) :: area_sum(IA,JA)
    real(RP) :: topo, mask

    character(len=H_LONG) :: fname

    real(RP) :: zerosw
    logical  :: hit_lat, hit_lon
    integer  :: index
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, iii, jjj, t
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_USERFILE",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO_USERFILE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVTOPO_USERFILE",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVTOPO_USERFILE",*) 'Not appropriate names in namelist PARAM_CNVTOPO_USERFILE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVTOPO_USERFILE)

    if ( USERFILE_NLAT <= 0 ) then
       LOG_ERROR("CNVTOPO_USERFILE",*) 'USERFILE_NLAT (number of latitude tile)  should be positive. Check! ', USERFILE_NLAT
       call PRC_abort
    endif

    if ( USERFILE_NLON <= 0 ) then
       LOG_ERROR("CNVTOPO_USERFILE",*) 'USERFILE_NLON (number of longitude tile) should be positive. Check! ', USERFILE_NLON
       call PRC_abort
    endif

    if ( USERFILE_DLAT <= 0.0_RP ) then
       LOG_ERROR("CNVTOPO_USERFILE",*) 'USERFILE_DLAT (width (deg.) of latitude tile) should be positive. Check! ', USERFILE_DLAT
       call PRC_abort
    endif

    if ( USERFILE_DLON <= 0.0_RP ) then
       LOG_ERROR("CNVTOPO_USERFILE",*) 'USERFILE_DLON (width (deg.) of longitude tile) should be positive. Check! ', USERFILE_DLON
       call PRC_abort
    endif

    if (       USERFILE_IN_CATALOGUE == '' &
         .AND. USERFILE_IN_FILENAME  == '' ) then
       LOG_ERROR("CNVTOPO_USERFILE",*) 'Neither catalogue file nor single file do not specified. Check!'
       call PRC_abort
    endif

    if    ( USERFILE_IN_DATATYPE == 'REAL8' ) then

       LOG_INFO("CNVTOPO_USERFILE",*) 'type of input data : REAL8'
       allocate( TILE_HEIGHT_orig   (USERFILE_NLON,USERFILE_NLAT) )
       allocate( TILE_HEIGHT_orig_DP(USERFILE_NLON,USERFILE_NLAT) )

    elseif( USERFILE_IN_DATATYPE == 'REAL4' ) then

       LOG_INFO("CNVTOPO_USERFILE",*) 'type of input data : REAL4'
       allocate( TILE_HEIGHT_orig   (USERFILE_NLON,USERFILE_NLAT) )
       allocate( TILE_HEIGHT_orig_SP(USERFILE_NLON,USERFILE_NLAT) )

    elseif( USERFILE_IN_DATATYPE == 'INT2'  ) then

       LOG_INFO("CNVTOPO_USERFILE",*) 'type of input data : INT2'
       allocate( TILE_HEIGHT_orig   (USERFILE_NLON,USERFILE_NLAT) )
       allocate( TILE_HEIGHT_orig_I2(USERFILE_NLON,USERFILE_NLAT) )

    else
       LOG_ERROR("CNVTOPO_USERFILE",*) 'Not appropriate type for USERFILE_IN_DATATYPE. Check!'
       LOG_ERROR_CONT(*) 'REAL8, REAL4, INT2 are available. requested: ', trim(USERFILE_IN_DATATYPE)
       call PRC_abort
    endif

    if    ( USERFILE_LATORDER_N2S ) then
       LOG_INFO("CNVTOPO_USERFILE",*) 'data ordar of the latitude direction : North -> South'
    else
       LOG_INFO("CNVTOPO_USERFILE",*) 'data ordar of the latitude direction : South -> North'
    endif

    do j = 1, JA
    do i = 1, IA
       area_sum(i,j) = 0.0_RP
       topo_sum(i,j) = 0.0_RP
    enddo
    enddo

    LONUY_mod(:,:) = mod( LONUY(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI

    DOMAIN_LATS    = minval(LATXV    (:,:))
    DOMAIN_LATE    = maxval(LATXV    (:,:))
    DOMAIN_LONS    = minval(LONUY_mod(:,:))
    DOMAIN_LONE    = maxval(LONUY_mod(:,:))
    DOMAIN_LONSLOC = minloc(LONUY_mod(:,:))
    DOMAIN_LONELOC = maxloc(LONUY_mod(:,:))

    check_IDL = .false.
    if (      DOMAIN_LONS < LONUY_mod( 0,DOMAIN_LONSLOC(2)) &
         .OR. DOMAIN_LONE > LONUY_mod(IA,DOMAIN_LONELOC(2)) ) then
       check_IDL = .true.
       DOMAIN_LONS = minval(LONUY_mod(:,:),mask=(LONUY_mod>0.0_RP))
       DOMAIN_LONE = maxval(LONUY_mod(:,:),mask=(LONUY_mod<0.0_RP))
    endif

    jos   = nint( USERFILE_DLAT / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    ios   = nint( USERFILE_DLON / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    jsize = USERFILE_NLAT * jos
    isize = USERFILE_NLON * ios

    allocate( TILE_HEIGHT(isize,jsize) )
    allocate( TILE_LATH  (0:jsize)     )
    allocate( TILE_LONH  (0:isize)     )

    LOG_INFO_CONT(*) 'Oversampling (j) orig = ', USERFILE_NLAT, ', use = ', jsize
    LOG_INFO_CONT(*) 'Oversampling (i) orig = ', USERFILE_NLON, ', use = ', isize

    TILE_DLAT_orig = USERFILE_DLAT * D2R
    TILE_DLON_orig = USERFILE_DLON * D2R
    LOG_INFO_CONT(*) 'TILE_DLAT       :', TILE_DLAT_orig/D2R
    LOG_INFO_CONT(*) 'TILE_DLON       :', TILE_DLON_orig/D2R

    TILE_DLAT = TILE_DLAT_orig / jos
    TILE_DLON = TILE_DLON_orig / ios
    LOG_INFO_CONT(*) 'TILE_DLAT (OS)  :', TILE_DLAT/D2R
    LOG_INFO_CONT(*) 'TILE_DLON (OS)  :', TILE_DLON/D2R

    !---< READ from external files >---

    if ( USERFILE_IN_CATALOGUE /= '' ) then

       ! input from catalogue file
       fname = trim(USERFILE_IN_DIR)//'/'//trim(USERFILE_IN_CATALOGUE)

       LOG_NEWLINE
       LOG_INFO("CNVTOPO_USERFILE",*) 'Input catalogue file:', trim(fname)

       fid = IO_get_available_fid()
       open( fid,                  &
             file   = trim(fname), &
             form   = 'formatted', &
             status = 'old',       &
             iostat = ierr         )

          if ( ierr /= 0 ) then
             LOG_ERROR("CNVTOPO_USERFILE",*) 'catalogue file not found! ', trim(fname)
             call PRC_abort
          endif

          do t = 1, TILE_nlim
             read(fid,*,iostat=ierr) index, TILE_LATS(t), TILE_LATE(t), & ! South->North
                                            TILE_LONS(t), TILE_LONE(t), & ! WEST->EAST
                                            TILE_fname(t)
             if ( ierr /= 0 ) exit

             if ( TILE_LONS(t) >= 180.0_RP ) then
                TILE_LONS(t) = TILE_LONS(t) - 360.0_RP
                TILE_LONE(t) = TILE_LONE(t) - 360.0_RP
             endif
             if ( TILE_LONS(t) < -180.0_RP ) TILE_LONS(t) = TILE_LONS(t) + 360.0_RP
             if ( TILE_LONE(t) < -180.0_RP ) TILE_LONE(t) = TILE_LONE(t) + 360.0_RP

          enddo

          TILE_nmax = t - 1
       close(fid)

    elseif( USERFILE_IN_FILENAME /= '' ) then

       ! input from single file
       TILE_nmax     = 1
       TILE_fname(1) = USERFILE_IN_FILENAME
       TILE_LATS (1) = USERFILE_LAT_START
       TILE_LATE (1) = USERFILE_LAT_END
       TILE_LONS (1) = USERFILE_LON_START
       TILE_LONE (1) = USERFILE_LON_END

    endif

    ! data file
    do t = 1, TILE_nmax
       hit_lat = .false.
       hit_lon = .false.

       if (      ( TILE_LATS(t)*D2R >= DOMAIN_LATS  .AND. TILE_LATS(t)*D2R < DOMAIN_LATE  ) &
            .OR. ( TILE_LATE(t)*D2R >= DOMAIN_LATS  .AND. TILE_LATE(t)*D2R < DOMAIN_LATE  ) ) then
          hit_lat = .true.
       endif

       if (      ( DOMAIN_LATS  >= TILE_LATS(t)*D2R .AND. DOMAIN_LATS  < TILE_LATE(t)*D2R ) &
            .OR. ( DOMAIN_LATE  >= TILE_LATS(t)*D2R .AND. DOMAIN_LATE  < TILE_LATE(t)*D2R ) ) then
          hit_lat = .true.
       endif

       if ( check_IDL ) then
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < PI           ) &
               .OR. ( TILE_LONS(t)*D2R >= -PI          .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < PI           ) &
               .OR. ( TILE_LONE(t)*D2R >= -PI          .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
             hit_lon = .true.
          endif
       else
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
             hit_lon = .true.
          endif
       endif

       if (      ( DOMAIN_LONS  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS  < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE  < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(USERFILE_IN_DIR)//'/'//trim(TILE_fname(t))

          LOG_NEWLINE
          LOG_INFO("CNVTOPO_USERFILE",*) 'Input data file :', trim(fname)
          LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if ( check_IDL ) then
             LOG_INFO_CONT(*) '(Date line exists within the domain)'
          endif
          LOG_INFO_CONT(*) 'Tile   (LAT)    :', TILE_LATS(t), TILE_LATE(t)
          LOG_INFO_CONT(*) '       (LON)    :', TILE_LONS(t), TILE_LONE(t)

          if ( USERFILE_IN_DATATYPE == 'REAL8' ) then

             fid = IO_get_available_fid()
             open( fid,                                    &
                   file   = trim(fname),                   &
                   form   = 'unformatted',                 &
                   access = 'direct',                      &
                   status = 'old',                         &
                   recl   = USERFILE_NLON*USERFILE_NLAT*8, &
                   iostat = ierr                           )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNVTOPO_USERFILE",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=1) TILE_HEIGHT_orig_DP(:,:)
             close(fid)

             if ( USERFILE_LATORDER_N2S ) then
                do j = 1, USERFILE_NLAT
                do i = 1, USERFILE_NLON
                   TILE_HEIGHT_orig(i,j) = real( TILE_HEIGHT_orig_DP(i,USERFILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, USERFILE_NLAT
                do i = 1, USERFILE_NLON
                   TILE_HEIGHT_orig(i,j) = real( TILE_HEIGHT_orig_DP(i,j), kind=RP )
                enddo
                enddo
             endif

          elseif( USERFILE_IN_DATATYPE == 'REAL4' ) then

             fid = IO_get_available_fid()
             open( fid,                                    &
                   file   = trim(fname),                   &
                   form   = 'unformatted',                 &
                   access = 'direct',                      &
                   status = 'old',                         &
                   recl   = USERFILE_NLON*USERFILE_NLAT*4, &
                   iostat = ierr                           )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNVTOPO_USERFILE",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=1) TILE_HEIGHT_orig_SP(:,:)
             close(fid)

             if ( USERFILE_LATORDER_N2S ) then
                do j = 1, USERFILE_NLAT
                do i = 1, USERFILE_NLON
                   TILE_HEIGHT_orig(i,j) = real( TILE_HEIGHT_orig_SP(i,USERFILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, USERFILE_NLAT
                do i = 1, USERFILE_NLON
                   TILE_HEIGHT_orig(i,j) = real( TILE_HEIGHT_orig_SP(i,j), kind=RP )
                enddo
                enddo
             endif

          elseif( USERFILE_IN_DATATYPE == 'INT2' ) then

             fid = IO_get_available_fid()
             open( fid,                                    &
                   file   = trim(fname),                   &
                   form   = 'unformatted',                 &
                   access = 'direct',                      &
                   status = 'old',                         &
                   recl   = USERFILE_NLON*USERFILE_NLAT*2, &
                   iostat = ierr                           )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNVTOPO_USERFILE",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=1) TILE_HEIGHT_orig_I2(:,:)
             close(fid)

             if ( USERFILE_LATORDER_N2S ) then
                do j = 1, USERFILE_NLAT
                do i = 1, USERFILE_NLON
                   TILE_HEIGHT_orig(i,j) = real( TILE_HEIGHT_orig_I2(i,USERFILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, USERFILE_NLAT
                do i = 1, USERFILE_NLON
                   TILE_HEIGHT_orig(i,j) = real( TILE_HEIGHT_orig_I2(i,j), kind=RP )
                enddo
                enddo
             endif

          endif

          ! oversampling
          do jj = 1, USERFILE_NLAT
          do ii = 1, USERFILE_NLON
             do j = 1, jos
             do i = 1, ios
                jjj = (jj-1) * jos + j
                iii = (ii-1) * ios + i

                TILE_HEIGHT(iii,jjj) = TILE_HEIGHT_orig(ii,jj)
             enddo
             enddo
          enddo
          enddo

          TILE_LATH(0) = TILE_LATS(t) * D2R
          do jj = 1, jsize
             TILE_LATH(jj) = TILE_LATH(jj-1) + TILE_DLAT
             !LOG_INFO("CNVTOPO_USERFILE",*) jj, TILE_LATH(jj) / D2R
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, isize
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
             !LOG_INFO("CNVTOPO_USERFILE",*) ii, TILE_LONH(ii) / D2R
          enddo

          ! match and calc fraction
          do jj = 1, jsize
          do ii = 1, isize

             iloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             ifrac_l = 1.0_RP

             jloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             jfrac_b = 1.0_RP

             if (      TILE_LATH(jj  ) < DOMAIN_LATS &
                  .OR. TILE_LATH(jj-1) > DOMAIN_LATE ) then
                cycle
             endif

             if ( check_IDL ) then
                if (       TILE_LONH(ii  ) < DOMAIN_LONS &
                     .AND. TILE_LONH(ii-1) > DOMAIN_LONE ) then
                   cycle
                endif
             else
                if (       TILE_LONH(ii  ) < DOMAIN_LONS &
                     .OR.  TILE_LONH(ii-1) > DOMAIN_LONE ) then
                   cycle
                endif
             endif

      jloop: do j = JS-1, JE+1
      iloop: do i = IS-1, IE+1
                if (       TILE_LONH(ii-1) >= LONUY_mod(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= LATXV    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  LATXV    (i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif

                if (       LONUY_mod(i-1,j) >= LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1)  >= LATXV    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1)  <  LATXV    (i  ,j  ) ) then ! across the IDL

                   if    (       TILE_LONH(ii-1) >= LONUY_mod(i-1,j) &
                           .AND. TILE_LONH(ii-1) <  PI               ) then

                      iloc    = i
                      ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1)+2.0_RP*PI, TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   elseif(       TILE_LONH(ii-1) >= -PI              &
                           .AND. TILE_LONH(ii-1) <  LONUY_mod(i  ,j) ) then

                      iloc    = i
                      ifrac_l = min( LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   endif

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             topo = TILE_HEIGHT(ii,jj)
             mask = 0.5_RP - sign( 0.5_RP,topo ) ! if Height is negative, mask = 1

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) ) * ( 1.0_RP - mask )

!             LOG_INFO("CNVTOPO_USERFILE",*) ii, jj, area, iloc, jloc, ifrac_l, jfrac_b, TILE_HEIGHT(ii,jj)

             area_fraction = (       ifrac_l) * (       jfrac_b) * area
             area_sum(iloc  ,jloc  ) = area_sum(iloc  ,jloc  ) + area_fraction
             topo_sum(iloc  ,jloc  ) = topo_sum(iloc  ,jloc  ) + area_fraction * topo

             area_fraction = (1.0_RP-ifrac_l) * (       jfrac_b) * area
             area_sum(iloc+1,jloc  ) = area_sum(iloc+1,jloc  ) + area_fraction
             topo_sum(iloc+1,jloc  ) = topo_sum(iloc+1,jloc  ) + area_fraction * topo

             area_fraction = (       ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc  ,jloc+1) = area_sum(iloc  ,jloc+1) + area_fraction
             topo_sum(iloc  ,jloc+1) = topo_sum(iloc  ,jloc+1) + area_fraction * topo

             area_fraction = (1.0_RP-ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc+1,jloc+1) = area_sum(iloc+1,jloc+1) + area_fraction
             topo_sum(iloc+1,jloc+1) = topo_sum(iloc+1,jloc+1) + area_fraction * topo
          enddo
          enddo

       endif
    enddo ! tile loop

    do j = JS, JE
    do i = IS, IE
       mask   = 0.5_RP + sign( 0.5_RP, area_sum(i,j)-EPS ) ! if any data is found, overwrite
       zerosw = 0.5_RP - sign( 0.5_RP, area_sum(i,j)-EPS )
       topo   = topo_sum(i,j) * ( 1.0_RP-zerosw ) / ( area_sum(i,j)-zerosw )
       TOPO_Zsfc(i,j) = (        mask ) * topo &         ! overwrite
                      + ( 1.0_RP-mask ) * TOPO_Zsfc(i,j) ! keep existing value
    enddo
    enddo

    return
  end subroutine CNVTOPO_USERFILE

  !-----------------------------------------------------------------------------
  !> check slope
  subroutine CNVTOPO_smooth( &
       Zsfc )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       FDX => ATMOS_GRID_CARTESC_FDX, &
       FDY => ATMOS_GRID_CARTESC_FDY
    use scale_comm_cartesC, only: &
       COMM_horizontal_max
    use scale_statistics, only: &
       STATISTICS_detail
    use scale_topography, only: &
       TOPO_fillhalo
    implicit none

    real(RP), intent(inout) :: Zsfc(IA,JA)

    real(RP) :: DZsfc_DXY(IA,JA,2) ! d(Zsfc)/dx at u-position and d(Zsfc)/dy at v-position

    real(RP) :: DXL(IA-1)
    real(RP) :: DYL(JA-1)

    real(RP) :: FLX_X(IA,JA)
    real(RP) :: FLX_Y(IA,JA)

    real(RP) :: slope(IA,JA)
    real(RP) :: maxslope
    real(RP) :: flag

    character(len=8), parameter :: varname(2) = (/ "DZsfc_DX", "DZsfc_DY" /)


    integer :: ite
    integer :: i, j
    !---------------------------------------------------------------------------

    if ( CNVTOPO_smooth_type == 'OFF' ) then
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Do not apply smoothing.'

       return
    else
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Apply smoothing.'
       LOG_INFO_CONT(*) 'Slope limit       = ', CNVTOPO_smooth_maxslope_limit
       LOG_INFO_CONT(*) 'Smoothing type    = ', CNVTOPO_smooth_type
       LOG_INFO_CONT(*) 'Smoothing locally = ', CNVTOPO_smooth_local
       LOG_NEWLINE
    endif

    DXL(:) = FDX(:)
    DYL(:) = FDY(:)



    ! digital filter
    do ite = 1, CNVTOPO_smooth_itelim+1
       LOG_PROGRESS(*) 'smoothing itelation : ', ite

       call TOPO_fillhalo( Zsfc=Zsfc(:,:), FILL_BND=.true. )

       do j = 1, JA
       do i = 1, IA-1
          DZsfc_DXY(i,j,1) = atan2( ( Zsfc(i+1,j)-Zsfc(i,j) ), DXL(i) ) / D2R
       enddo
       enddo
       DZsfc_DXY(IA,:,1) = 0.0_RP
       do j = 1, JA-1
       do i = 1, IA
          DZsfc_DXY(i,j,2) = atan2( ( Zsfc(i,j+1)-Zsfc(i,j) ), DYL(j) ) / D2R
       enddo
       enddo
       DZsfc_DXY(:,JA,2) = 0.0_RP

       slope(:,:) = max( abs(DZsfc_DXY(:,:,1)), abs(DZsfc_DXY(:,:,2)) )
       call COMM_horizontal_max( maxslope, slope(:,:) )

       LOG_PROGRESS(*) 'maximum slope [deg] : ', maxslope

       if( maxslope < CNVTOPO_smooth_maxslope_limit ) exit

       call STATISTICS_detail( IA, IS, IE, JA, JS, JE, 2, &
                               varname(:), DZsfc_DXY(:,:,:) )

       select case( CNVTOPO_smooth_type )
       case( 'GAUSSIAN' )

          ! 3 by 3 gaussian filter
          do j = JS, JE
          do i = IS, IE
             Zsfc(i,j) = ( 0.2500_RP * Zsfc(i  ,j  ) &
                         + 0.1250_RP * Zsfc(i-1,j  ) &
                         + 0.1250_RP * Zsfc(i+1,j  ) &
                         + 0.1250_RP * Zsfc(i  ,j-1) &
                         + 0.1250_RP * Zsfc(i  ,j+1) &
                         + 0.0625_RP * Zsfc(i-1,j-1) &
                         + 0.0625_RP * Zsfc(i+1,j-1) &
                         + 0.0625_RP * Zsfc(i-1,j+1) &
                         + 0.0625_RP * Zsfc(i+1,j+1) )
          enddo
          enddo

       case( 'LAPLACIAN' )

          do j = JS  , JE
          do i = IS-1, IE
             FLX_X(i,j) = Zsfc(i+1,j) - Zsfc(i,j)
!             FLX_TMP(i,j) = Zsfc(i+1,j) - Zsfc(i,j)
          enddo
          enddo
!!$          call TOPO_fillhalo( FLX_TMP )
!!$          do j = JS  , JE
!!$          do i = IS-1, IE
!!$             FLX_X(i,j) = - ( FLX_TMP(i+1,j) - FLX_TMP(i,j) )
!!$          enddo
!!$          enddo

          do j = JS-1, JE
          do i = IS  , IE
             FLX_Y(i,j) = Zsfc(i,j+1) - Zsfc(i,j)
!             FLX_TMP(i,j) = Zsfc(i,j+1) - Zsfc(i,j)
          enddo
          enddo
!!$          call TOPO_fillhalo( FLX_TMP )
!!$          do j = JS-1, JE
!!$          do i = IS  , IE
!!$             FLX_Y(i,j) = - ( FLX_TMP(i,j+1) - FLX_TMP(i,j) )
!!$          enddo
!!$          enddo


          if ( CNVTOPO_smooth_local ) then
             do j = JS  , JE
             do i = IS-1, IE
                flag = 0.5_RP &
                     + sign(0.5_RP, max( abs(DZsfc_DXY(i+1,j  ,1)), &
                                         abs(DZsfc_DXY(i  ,j  ,1)), &
                                         abs(DZsfc_DXY(i-1,j  ,1)), &
                                         abs(DZsfc_DXY(i+1,j  ,2)), &
                                         abs(DZsfc_DXY(i+1,j-1,2)), &
                                         abs(DZsfc_DXY(i  ,j  ,2)), &
                                         abs(DZsfc_DXY(i  ,j-1,2))  &
                                       ) - CNVTOPO_smooth_maxslope_limit )
                FLX_X(i,j) = FLX_X(i,j) * flag
             enddo
             enddo
             do j = JS-1, JE
             do i = IS  , IE
                flag = 0.5_RP &
                     + sign(0.5_RP, max( abs(DZsfc_DXY(i  ,j+1,2)), &
                                         abs(DZsfc_DXY(i  ,j  ,2)), &
                                         abs(DZsfc_DXY(i  ,j-1,2)), &
                                         abs(DZsfc_DXY(i  ,j+1,1)), &
                                         abs(DZsfc_DXY(i-1,j+1,1)), &
                                         abs(DZsfc_DXY(i  ,j  ,1)), &
                                         abs(DZsfc_DXY(i-1,j  ,1))  &
                                       ) - CNVTOPO_smooth_maxslope_limit )
                FLX_Y(i,j) = FLX_Y(i,j) * flag
             enddo
             enddo
          endif

          do j = JS, JE
          do i = IS, IE
             Zsfc(i,j) = max( 0.0_RP, &
                              Zsfc(i,j) &
                              + 0.1_RP * ( ( FLX_X(i,j) - FLX_X(i-1,j) ) &
                                         + ( FLX_Y(i,j) - FLX_Y(i,j-1) ) ) )
          enddo
          enddo

       case default
          LOG_ERROR("CNVTOPO_smooth",*) 'Invalid smoothing type'
          call PRC_abort
       end select

    enddo

    if ( ite  > CNVTOPO_smooth_itelim ) then
       LOG_ERROR("CNVTOPO_smooth",*) 'not converged'
       call PRC_abort
    else
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'smoothing complete.'
    endif



    if ( CNVTOPO_smooth_hypdiff_niter > 0 ) then

       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Apply hyperdiffusion.'

       call CNVTOPO_hypdiff( Zsfc(:,:), CNVTOPO_smooth_hypdiff_niter )

       do j = 1, JA
       do i = 1, IA-1
          DZsfc_DXY(i,j,1) = atan2( ( Zsfc(i+1,j)-Zsfc(i,j) ), DXL(i) ) / D2R
       enddo
       enddo
       DZsfc_DXY(IA,:,1) = 0.0_RP
       do j = 1, JA-1
       do i = 1, IA
          DZsfc_DXY(i,j,2) = atan2( ( Zsfc(i,j+1)-Zsfc(i,j) ), DYL(j) ) / D2R
       enddo
       enddo
       DZsfc_DXY(:,JA,2) = 0.0_RP

       slope(:,:) = max( abs(DZsfc_DXY(:,:,1)), abs(DZsfc_DXY(:,:,2)) )
       call COMM_horizontal_max( maxslope, slope(:,:) )

       LOG_INFO("CNVTOPO_smooth",*) 'maximum slope [deg] : ', maxslope

    end if

    call TOPO_fillhalo( Zsfc=Zsfc(:,:), FILL_BND=.true. )

    call STATISTICS_detail( IA, IS, IE, JA, JS, JE, 2, &
                            varname(:), DZsfc_DXY(:,:,:) )

    LOG_NEWLINE

    return
  end subroutine CNVTOPO_smooth

  subroutine CNVTOPO_hypdiff( Zsfc, nite )
    use scale_topography, only: &
       TOPO_fillhalo
    real(RP), intent(inout) :: Zsfc(IA,JA)
    integer,  intent(in)    :: nite

    real(RP), pointer :: p1(:,:)
    real(RP), pointer :: p2(:,:)
    real(RP), target :: work1(IA,JA)
    real(RP), target :: work2(IA,JA)

    integer :: i, j
    integer :: ite, n

    ! reduce grid-scale variation
    do ite = 1, nite
       call TOPO_fillhalo( Zsfc=Zsfc(:,:), FILL_BND=.true. )
       work2(:,:) = Zsfc(:,:)
       p1 => work2
       p2 => work1
       do n = 1, 4 ! 8th derivative
!       do n = 1, 2 ! 4th derivative
          do j = JS, JE
          do i = IS, IE
             p2(i,j) = ( - p1(i+1,j) + p1(i,j)*2.0_RP - p1(i-1,j) &
                         - p1(i,j+1) + p1(i,j)*2.0_RP - p1(i,j-1) ) / 8.0_RP
          end do
          end do
          call TOPO_fillhalo( Zsfc=p2(:,:), FILL_BND=.true. )
          if ( mod(n,2) == 0 ) then
             p1 => work2
             p2 => work1
          else
             p1 => work1
             p2 => work2
          end if
       end do
       do j = JS, JE
       do i = IS, IE
          Zsfc(i,j) = max( Zsfc(i,j) - p1(i,j), 0.0_RP )
       end do
       end do
    end do

    return
  end subroutine CNVTOPO_hypdiff

end module mod_cnvtopo
