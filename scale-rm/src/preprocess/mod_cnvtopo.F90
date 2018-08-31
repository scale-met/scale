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
  integer,  private :: CNVTOPO_smooth_hypdiff_order  = 4
  integer,  private :: CNVTOPO_smooth_hypdiff_niter  = 20
  logical,  private :: CNVTOPO_smooth_local          = .true.
  integer,  private :: CNVTOPO_smooth_itelim         = 10000
  logical,  private :: CNVTOPO_smooth_trim_ocean     = .true.
  real(RP), private :: CNVTOPO_smooth_maxslope_ratio = 10.0_RP ! ratio of DZDX, DZDY
  real(RP), private :: CNVTOPO_smooth_maxslope       = -1.0_RP ! [deg]
  real(RP), private :: CNVTOPO_smooth_maxslope_limit

  logical,  private :: CNVTOPO_copy_parent           = .false.

  integer,  private :: CNVTOPO_interp_level          = 5

  real(RP), private :: DOMAIN_LATS, DOMAIN_LATE
  real(RP), private :: DOMAIN_LONS, DOMAIN_LONE

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
    use scale_statistics, only: &
       STATISTICS_horizontal_min
    use scale_atmos_grid_cartesC, only: &
       CDZ => ATMOS_GRID_CARTESC_CDZ, &
       FDX => ATMOS_GRID_CARTESC_FDX, &
       FDY => ATMOS_GRID_CARTESC_FDY
    implicit none

    character(len=H_SHORT) :: CNVTOPO_name = 'NONE' ! keep backward compatibility

    namelist / PARAM_CNVTOPO / &
       CNVTOPO_name,                  &
       CNVTOPO_UseGTOPO30,            &
!       CNVTOPO_UseGMTED2010,          &
       CNVTOPO_UseDEM50M,             &
       CNVTOPO_UseUSERFILE,           &
       CNVTOPO_smooth_trim_ocean,     &
       CNVTOPO_smooth_hypdiff_order,  &
       CNVTOPO_smooth_hypdiff_niter,  &
       CNVTOPO_smooth_maxslope_ratio, &
       CNVTOPO_smooth_maxslope,       &
       CNVTOPO_smooth_local,          &
       CNVTOPO_smooth_itelim,         &
       CNVTOPO_smooth_type,           &
       CNVTOPO_copy_parent,           &
       CNVTOPO_interp_level

    real(RP) :: minslope(IA,JA)
    real(RP) :: DXL(IA-1)
    real(RP) :: DYL(JA-1)
    real(RP) :: DZDX, DZDY

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
!!$    case('GMTED2010')
!!$       CNVTOPO_UseGTOPO30   = .false.
!!$       CNVTOPO_UseGMTED2010 = .true.
!!$       CNVTOPO_UseDEM50M    = .false.
!!$       CNVTOPO_UseUSERFILE  = .false.
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
      !$omp parallel do &
      !$omp private(DZDX,DZDY)
      do j = JS, JE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DYL(j) ) / D2R
         minslope(IS,j) = min( minslope(IS,j), DZDX, DZDY )
      enddo
      enddo

      !$omp parallel do &
      !$omp private(DZDX,DZDY)
      do j = JS, JE
      do i = IS, IE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * CDZ(k), DYL(j) ) / D2R
         minslope(i,j) = min( minslope(i,j), DZDX, DZDY )
      enddo
      enddo
      enddo

      call STATISTICS_horizontal_min( IA, IS, IE, JA, JS, JE, &
                                      minslope(:,:), CNVTOPO_smooth_maxslope_limit )
    end if

    return
  end subroutine CNVTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVTOPO
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R   => CONST_D2R
    use scale_prc, only: &
       PRC_abort
    use scale_topography, only: &
       TOPO_fillhalo, &
       TOPO_Zsfc, &
       TOPO_write
    use mod_copytopo, only: &
       COPYTOPO
    use scale_atmos_grid_cartesC_real, only: &
       LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       LONUY => ATMOS_GRID_CARTESC_REAL_LONUY
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( CNVTOPO_DoNothing ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'skip  convert topography data'
    else
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start convert topography data'

       DOMAIN_LATS = minval( LATXV(:,:) )
       DOMAIN_LATE = maxval( LATXV(:,:) )
       DOMAIN_LONS = minval( LONUY(:,:) )
       DOMAIN_LONE = maxval( LONUY(:,:) )

       LOG_INFO("CNVTOPO",*) 'Domain Information'
       LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
       LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R

       !$omp parallel do
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
          TOPO_Zsfc(i,j) = 0.0_RP
       end do
       end do

       if ( CNVTOPO_UseGTOPO30 ) then
          call CNVTOPO_GTOPO30( TOPO_Zsfc(:,:) ) ! [INOUT]
       endif

       if ( CNVTOPO_UseGMTED2010 ) then
          call CNVTOPO_GMTED2010( TOPO_Zsfc(:,:) ) ! [INOUT]
       endif

       if ( CNVTOPO_UseDEM50M ) then
          call CNVTOPO_DEM50M( TOPO_Zsfc(:,:) ) ! [INOUT]
       endif

       if ( CNVTOPO_UseUSERFILE ) then
          call CNVTOPO_USERFILE( TOPO_Zsfc(:,:) ) ! [INOUT]
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
  subroutine CNVTOPO_GTOPO30( TOPO_Zsfc )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF  => CONST_UNDEF, &
       D2R    => CONST_D2R
    use scale_atmos_grid_cartesC_real, only: &
       LAT   => ATMOS_GRID_CARTESC_REAL_LAT, &
       LON   => ATMOS_GRID_CARTESC_REAL_LON
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info, &
       FILE_TILEDATA_get_data
    use scale_interp, only: &
       INTERP_factor2d, &
       INTERP_interp2d
    implicit none
    real(RP), intent(inout) :: TOPO_Zsfc(IA,JA)

    character(len=H_LONG) :: GTOPO30_IN_DIR       = '.' !< directory contains GTOPO30 files (GrADS format)
    character(len=H_LONG) :: GTOPO30_IN_CATALOGUE = ''  !< metadata files for GTOPO30

    namelist / PARAM_CNVTOPO_GTOPO30 / &
       GTOPO30_IN_DIR,       &
       GTOPO30_IN_CATALOGUE

    ! GTOPO30 data
    real(RP), parameter   :: GTOPO30_DLAT = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.
    real(RP), parameter   :: GTOPO30_DLON = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.

    ! topo
    real(RP) :: Zsfc(IA,JA)

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 100
    integer               :: TILE_nmax
    character(len=H_LONG) :: TILE_fname(TILE_nlim)
    logical               :: TILE_hit  (TILE_nlim)
    integer               :: TILE_JS   (TILE_nlim)
    integer               :: TILE_JE   (TILE_nlim)
    integer               :: TILE_IS   (TILE_nlim)
    integer               :: TILE_IE   (TILE_nlim)
    real(RP)              :: TILE_DLAT, TILE_DLON

    real(RP), allocatable :: HEIGHT(:,:)
    real(RP), allocatable :: LATH  (:,:)
    real(RP), allocatable :: LONH  (:,:)
    integer               :: nLONH, nLATH

    integer  :: GLOBAL_IA

    ! interpolation
    integer,  allocatable :: idx_i(:,:,:)
    integer,  allocatable :: idx_j(:,:,:)
    real(RP), allocatable :: hfact(:,:,:)

    character(len=H_LONG) :: fname

    integer :: ish, ieh, jsh, jeh
    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

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


    TILE_DLAT = GTOPO30_DLAT * D2R
    TILE_DLON = GTOPO30_DLON * D2R

    ! catalogue file
    fname = trim(GTOPO30_IN_DIR)//'/'//trim(GTOPO30_IN_CATALOGUE)

    call FILE_TILEDATA_get_info( TILE_nlim,                                          & ! [IN]
                                 TILE_DLAT, TILE_DLON,                               & ! [IN]
                                 DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                 fname,                                              & ! [IN]
                                 GLOBAL_IA,                                          & ! [OUT]
                                 TILE_nmax,                                          & ! [OUT]
                                 TILE_fname(:), TILE_hit(:),                         & ! [OUT]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                 nLATH, nLONH, jsh, jeh, ish, ieh                    ) ! [OUT]

    allocate( HEIGHT(nLONH,nLATH) )
    allocate( LATH  (nLONH,nLATH) )
    allocate( LONH  (nLONH,nLATH) )

    call FILE_TILEDATA_get_data( nLATH, nLONH,                                   & ! [IN]
                                 GTOPO30_IN_DIR,                                 & ! [IN]
                                 GLOBAL_IA,                                      & ! [IN]
                                 TILE_nmax,                                      & ! [IN]
                                 TILE_DLAT, TILE_DLON,                           & ! [IN]
                                 TILE_fname(:), TILE_hit(:),                     & ! [IN]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:), & ! [IN]
                                 jsh, jeh, ish, ieh,                             & ! [IN]
                                 "INT2",                                         & ! [IN]
                                 HEIGHT(:,:), LATH(:,:), LONH(:,:),              & ! [OUT]
                                 min_value = -9000.0_RP, yrevers = .true.        ) ! [IN]

    ! interporation
    allocate( idx_i(IA,JA,CNVTOPO_interp_level) )
    allocate( idx_j(IA,JA,CNVTOPO_interp_level) )
    allocate( hfact(IA,JA,CNVTOPO_interp_level) )

    call INTERP_factor2d( CNVTOPO_interp_level,       & ! [IN]
                          nLONH, nLATH,               & ! [IN]
                          LONH(:,:), LATH(:,:),       & ! [IN]
                          IA, JA,                     & ! [IN]
                          LON(:,:), LAT(:,:),         & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [OUT]
                          hfact(:,:,:)                ) ! [OUT]

    call INTERP_interp2d( CNVTOPO_interp_level,       & ! [IN]
                          nLONH,nLATH,                & ! [IN]
                          IA, JA,                     & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [IN]
                          hfact(:,:,:),               & ! [IN]
                          HEIGHT(:,:),                & ! [IN]
                          Zsfc(:,:)                   ) ! [OUT]


    deallocate( HEIGHT, LATH, LONH )
    deallocate( idx_j, idx_i, hfact )

    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       if ( Zsfc(i,j) /= UNDEF ) TOPO_Zsfc(i,j) = Zsfc(i,j) ! replace data
    end do
    end do

    return
  end subroutine CNVTOPO_GTOPO30

  !-----------------------------------------------------------------------------
  !> Convert from GMTED 2010 7.5 arc-seconds mesh
  subroutine CNVTOPO_GMTED2010( TOPO_Zsfc )
    implicit none
    real(RP), intent(inout) :: TOPO_Zsfc(IA,JA)
    !---------------------------------------------------------------------------

    return
  end subroutine CNVTOPO_GMTED2010

  !-----------------------------------------------------------------------------
  !> Convert from DEM 50m mesh
  subroutine CNVTOPO_DEM50M( TOPO_Zsfc )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF  => CONST_UNDEF, &
       RADIUS => CONST_RADIUS, &
       D2R    => CONST_D2R
    use scale_atmos_grid_cartesC, only: &
       DX, &
       DY
    use scale_atmos_grid_cartesC_real, only: &
       LAT   => ATMOS_GRID_CARTESC_REAL_LAT, &
       LON   => ATMOS_GRID_CARTESC_REAL_LON
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info, &
       FILE_TILEDATA_get_data
    use scale_interp, only: &
       INTERP_factor2d, &
       INTERP_interp2d
    implicit none
    real(RP), intent(inout) :: TOPO_Zsfc(IA,JA)

    character(len=H_LONG) :: DEM50M_IN_DIR       = '.' !< directory contains DEM50M files (GrADS format)
    character(len=H_LONG) :: DEM50M_IN_CATALOGUE = ''  !< metadata files for DEM50M

    namelist / PARAM_CNVTOPO_DEM50M / &
       DEM50M_IN_DIR,      &
       DEM50M_IN_CATALOGUE

    real(RP), parameter   :: DEM50M_DLAT = 5.0_RP / 60.0_RP / 200.0_RP ! 30 arc sec.
    real(RP), parameter   :: DEM50M_DLON = 7.5_RP / 60.0_RP / 200.0_RP ! 30 arc sec.

    ! topo
    real(RP) :: Zsfc(IA,JA)

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 1000
    integer               :: TILE_nmax
    character(len=H_LONG) :: TILE_fname(TILE_nlim)
    logical               :: TILE_hit  (TILE_nlim)
    integer               :: TILE_JS   (TILE_nlim)
    integer               :: TILE_JE   (TILE_nlim)
    integer               :: TILE_IS   (TILE_nlim)
    integer               :: TILE_IE   (TILE_nlim)
    real(RP)              :: TILE_DLAT, TILE_DLON

    real(RP), allocatable :: HEIGHT(:,:)
    real(RP), allocatable :: LATH  (:,:)
    real(RP), allocatable :: LONH  (:,:)
    integer               :: nLONH, nLATH

    integer  :: GLOBAL_IA

    ! interpolation
    integer,  allocatable :: idx_i(:,:,:)
    integer,  allocatable :: idx_j(:,:,:)
    real(RP), allocatable :: hfact(:,:,:)
    real(RP) :: search_limit

    character(len=H_LONG) :: fname

    integer :: ish, ieh, jsh, jeh
    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    !--- read namelist
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

    TILE_DLAT = DEM50M_DLAT * D2R
    TILE_DLON = DEM50M_DLON * D2R

    ! catalogue file
    fname = trim(DEM50M_IN_DIR)//'/'//trim(DEM50M_IN_CATALOGUE)

    call FILE_TILEDATA_get_info( TILE_nlim,                                          & ! [IN]
                                 TILE_DLAT, TILE_DLON,                               & ! [IN]
                                 DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                 fname,                                              & ! [IN]
                                 GLOBAL_IA,                                          & ! [OUT]
                                 TILE_nmax,                                          & ! [OUT]
                                 TILE_fname(:), TILE_hit(:),                         & ! [OUT]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                 nLATH, nLONH, jsh, jeh, ish, ieh                    ) ! [OUT]

    if ( .not. any(TILE_hit(1:TILE_nmax) ) ) return

    allocate( HEIGHT(nLONH,nLATH) )
    allocate( LATH  (nLONH,nLATH) )
    allocate( LONH  (nLONH,nLATH) )

    call FILE_TILEDATA_get_data( nLATH, nLONH,                                   & ! [IN]
                                 DEM50M_IN_DIR,                                  & ! [IN]
                                 GLOBAL_IA,                                      & ! [IN]
                                 TILE_nmax,                                      & ! [IN]
                                 TILE_DLAT, TILE_DLON,                           & ! [IN]
                                 TILE_fname(:), TILE_hit(:),                     & ! [IN]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:), & ! [IN]
                                 jsh, jeh, ish, ieh,                             & ! [IN]
                                 "REAL4",                                        & ! [IN]
                                 HEIGHT(:,:), LATH(:,:), LONH(:,:),              & ! [OUT]
                                 min_value = -900.0_RP                           ) ! [IN]

    ! interporation
    allocate( idx_i(IA,JA,CNVTOPO_interp_level) )
    allocate( idx_j(IA,JA,CNVTOPO_interp_level) )
    allocate( hfact(IA,JA,CNVTOPO_interp_level) )

    search_limit = max( sqrt(TILE_DLON**2 + TILE_DLAT**2) * RADIUS * 1.5_RP, sqrt(DX**2 + DY**2) )
    call INTERP_factor2d( CNVTOPO_interp_level,       & ! [IN]
                          nLONH, nLATH,               & ! [IN]
                          LONH(:,:), LATH(:,:),       & ! [IN]
                          IA, JA,                     & ! [IN]
                          LON(:,:), LAT(:,:),         & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [OUT]
                          hfact(:,:,:),               & ! [OUT]
                          latlon_structure = .true.,  & ! [OUT]
                          search_limit = search_limit ) ! [IN]

    call INTERP_interp2d( CNVTOPO_interp_level,       & ! [IN]
                          nLONH,nLATH,                & ! [IN]
                          IA, JA,                     & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [IN]
                          hfact(:,:,:),               & ! [IN]
                          HEIGHT(:,:),                & ! [IN]
                          Zsfc(:,:)                   ) ! [OUT]


    deallocate( HEIGHT, LATH, LONH )
    deallocate( idx_j, idx_i, hfact )

    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       if ( Zsfc(i,j) /= UNDEF ) TOPO_Zsfc(i,j) = Zsfc(i,j) ! replace data
    end do
    end do

    return
  end subroutine CNVTOPO_DEM50M

  !-----------------------------------------------------------------------------
  !> Convert from User-defined file
  subroutine CNVTOPO_USERFILE( TOPO_Zsfc )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF  => CONST_UNDEF, &
       RADIUS => CONST_RADIUS, &
       D2R    => CONST_D2R
    use scale_atmos_grid_cartesC, only: &
       DX, &
       DY
    use scale_atmos_grid_cartesC_real, only: &
       LAT   => ATMOS_GRID_CARTESC_REAL_LAT, &
       LON   => ATMOS_GRID_CARTESC_REAL_LON
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info, &
       FILE_TILEDATA_get_data
    use scale_interp, only: &
       INTERP_factor2d, &
       INTERP_interp2d
    implicit none
    real(RP), intent(inout) :: TOPO_Zsfc(IA,JA)

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

    ! topo
    real(RP) :: Zsfc(IA,JA)

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 1000
    integer               :: TILE_nmax
    character(len=H_LONG) :: TILE_fname(TILE_nlim)
    logical               :: TILE_hit  (TILE_nlim)
    integer               :: TILE_JS   (TILE_nlim)
    integer               :: TILE_JE   (TILE_nlim)
    integer               :: TILE_IS   (TILE_nlim)
    integer               :: TILE_IE   (TILE_nlim)
    real(RP)              :: TILE_DLAT, TILE_DLON

    real(RP), allocatable :: HEIGHT(:,:)
    real(RP), allocatable :: LATH  (:,:)
    real(RP), allocatable :: LONH  (:,:)
    integer               :: nLONH, nLATH

    integer  :: GLOBAL_IA

    ! interpolation
    integer,  allocatable :: idx_i(:,:,:)
    integer,  allocatable :: idx_j(:,:,:)
    real(RP), allocatable :: hfact(:,:,:)
    real(RP) :: search_limit

    character(len=H_LONG) :: fname

    real(RP) :: LATS, LATE, LONS, LONE
    logical  :: yrevers

    integer :: ish, ieh, jsh, jeh
    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    !--- read namelist

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
       LOG_ERROR("CNVTOPO_USERFILE",*) 'Neither catalogue file nor single file specified. Check!'
       call PRC_abort
    endif

    if    ( USERFILE_LATORDER_N2S ) then
       LOG_INFO("CNVTOPO_USERFILE",*) 'data ordar of the latitude direction : North -> South'
       yrevers = .true.
    else
       LOG_INFO("CNVTOPO_USERFILE",*) 'data ordar of the latitude direction : South -> North'
       yrevers = .false.
    endif



    TILE_DLAT = USERFILE_DLAT * D2R
    TILE_DLON = USERFILE_DLON * D2R

    if ( USERFILE_IN_CATALOGUE /= "" ) then
       ! catalogue file
       fname = trim(USERFILE_IN_DIR)//'/'//trim(USERFILE_IN_CATALOGUE)

       call FILE_TILEDATA_get_info( TILE_nlim,                                          & ! [IN]
                                    TILE_DLAT, TILE_DLON,                               & ! [IN]
                                    DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                    fname,                                              & ! [IN]
                                    GLOBAL_IA,                                          & ! [OUT]
                                    TILE_nmax,                                          & ! [OUT]
                                    TILE_fname(:), TILE_hit(:),                         & ! [OUT]
                                    TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                    nLATH, nLONH, jsh, jeh, ish, ieh                    ) ! [OUT]
    else
       LATS = USERFILE_LAT_START * D2R
       LATE = USERFILE_LAT_END   * D2R
       LONS = USERFILE_LON_START * D2R
       LONE = USERFILE_LON_END   * D2R
       call FILE_TILEDATA_get_info( TILE_nlim,                                          & ! [IN]
                                    TILE_DLAT, TILE_DLON,                               & ! [IN]
                                    DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                    '',                                                 & ! [IN]
                                    GLOBAL_IA,                                          & ! [OUT]
                                    TILE_nmax,                                          & ! [OUT]
                                    TILE_fname(:), TILE_hit(:),                         & ! [OUT]
                                    TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                    nLATH, nLONH, jsh, jeh, ish, ieh,                   & ! [OUT]
                                    single_fname = USERFILE_IN_FILENAME,                & ! [IN]
                                    LATS = LATS, LATE = LATE, LONS = LONS, LONE = LONE  ) ! [IN]
    end if

    if ( .not. any(TILE_hit(1:TILE_nmax) ) ) return

    allocate( HEIGHT(nLONH,nLATH) )
    allocate( LATH  (nLONH,nLATH) )
    allocate( LONH  (nLONH,nLATH) )

    call FILE_TILEDATA_get_data( nLATH, nLONH,                                   & ! [IN]
                                 USERFILE_IN_DIR,                                & ! [IN]
                                 GLOBAL_IA,                                      & ! [IN]
                                 TILE_nmax,                                      & ! [IN]
                                 TILE_DLAT, TILE_DLON,                           & ! [IN]
                                 TILE_fname(:), TILE_hit(:),                     & ! [IN]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:), & ! [IN]
                                 jsh, jeh, ish, ieh,                             & ! [IN]
                                 USERFILE_IN_DATATYPE,                           & ! [IN]
                                 HEIGHT(:,:), LATH(:,:), LONH(:,:),              & ! [OUT]
                                 min_value = 0.0_RP, yrevers = yrevers           ) ! [IN]

    ! interporation
    allocate( idx_i(IA,JA,CNVTOPO_interp_level) )
    allocate( idx_j(IA,JA,CNVTOPO_interp_level) )
    allocate( hfact(IA,JA,CNVTOPO_interp_level) )

    search_limit = max( sqrt(TILE_DLON**2 + TILE_DLAT**2) * RADIUS * 1.5_RP, sqrt(DX**2 + DY**2) )
    call INTERP_factor2d( CNVTOPO_interp_level,       & ! [IN]
                          nLONH, nLATH,               & ! [IN]
                          LONH(:,:), LATH(:,:),       & ! [IN]
                          IA, JA,                     & ! [IN]
                          LON(:,:), LAT(:,:),         & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [OUT]
                          hfact(:,:,:),               & ! [OUT]
                          latlon_structure = .true.,  & ! [OUT]
                          search_limit = search_limit ) ! [IN]

    call INTERP_interp2d( CNVTOPO_interp_level,       & ! [IN]
                          nLONH,nLATH,                & ! [IN]
                          IA, JA,                     & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [IN]
                          hfact(:,:,:),               & ! [IN]
                          HEIGHT(:,:),                & ! [IN]
                          Zsfc(:,:)                   ) ! [OUT]

    deallocate( HEIGHT, LATH, LONH )
    deallocate( idx_j, idx_i, hfact )

    !$omp parallel do
    do j = 1, JA
    do i = 1, IA
       if ( Zsfc(i,j) /= UNDEF ) TOPO_Zsfc(i,j) = Zsfc(i,j) ! replace data
    end do
    end do

    return
  end subroutine CNVTOPO_USERFILE

  !-----------------------------------------------------------------------------
  !> check slope
  subroutine CNVTOPO_smooth( &
       Zsfc )
    use scale_const, only: &
       EPS => CONST_EPS, &
       D2R => CONST_D2R
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       FDX => ATMOS_GRID_CARTESC_FDX, &
       FDY => ATMOS_GRID_CARTESC_FDY
    use scale_statistics, only: &
       STATISTICS_detail, &
       STATISTICS_horizontal_max
    use scale_topography, only: &
       TOPO_fillhalo
    use scale_filter, only: &
       FILTER_hyperdiff
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    real(RP), intent(inout) :: Zsfc(IA,JA)

    real(RP) :: DZsfc_DXY(IA,JA,2) ! d(Zsfc)/dx at u-position and d(Zsfc)/dy at v-position

    real(RP) :: DXL(IA-1)
    real(RP) :: DYL(JA-1)

    real(RP) :: FLX_X(IA,JA)
    real(RP) :: FLX_Y(IA,JA)

    real(RP) :: slope(IA,JA)
    real(RP) :: maxslope
    real(RP), pointer :: TOPO_sign(:,:)
    real(RP) :: flag, ocean_flag

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

    if ( CNVTOPO_smooth_trim_ocean ) then
       allocate( TOPO_sign(IA,JA) )
       !$omp parallel do &
       !$omp private(ocean_flag)
       do j = 1, JA
       do i = 1, IA
          ocean_flag = ( 0.5_RP + sign( 0.5_RP, LANDUSE_fact_ocean(i,j) - 1.0_RP + EPS ) ) & ! fact_ocean==1
                     * ( 0.5_RP + sign( 0.5_RP, EPS - abs(Zsfc(i,j)) ) )                   ! |Zsfc| < EPS
          TOPO_sign(i,j) = sign( 1.0_RP, Zsfc(i,j) ) * ( 1.0_RP - ocean_flag )
       end do
       end do
    else
       TOPO_sign => NULL()
    end if

    ! digital filter
    do ite = 1, CNVTOPO_smooth_itelim+1
       LOG_PROGRESS(*) 'smoothing itelation : ', ite

       call TOPO_fillhalo( Zsfc=Zsfc(:,:), FILL_BND=.true. )

       !$omp parallel do
       do j = 1, JA
       do i = 1, IA-1
          DZsfc_DXY(i,j,1) = atan2( ( Zsfc(i+1,j)-Zsfc(i,j) ), DXL(i) ) / D2R
       enddo
       enddo
       DZsfc_DXY(IA,:,1) = 0.0_RP
       !$omp parallel do
       do j = 1, JA-1
       do i = 1, IA
          DZsfc_DXY(i,j,2) = atan2( ( Zsfc(i,j+1)-Zsfc(i,j) ), DYL(j) ) / D2R
       enddo
       enddo
       DZsfc_DXY(:,JA,2) = 0.0_RP

       slope(:,:) = max( abs(DZsfc_DXY(:,:,1)), abs(DZsfc_DXY(:,:,2)) )
       call STATISTICS_horizontal_max( IA, IS, IE, JA, JS, JE, &
                                       slope(:,:), maxslope )

       LOG_PROGRESS(*) 'maximum slope [deg] : ', maxslope

       if( maxslope < CNVTOPO_smooth_maxslope_limit ) exit

       call STATISTICS_detail( IA, IS, IE, JA, JS, JE, 2, &
                               varname(:), DZsfc_DXY(:,:,:) )

       select case( CNVTOPO_smooth_type )
       case( 'GAUSSIAN' )

          ! 3 by 3 gaussian filter
          !$omp parallel do
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

          !$omp parallel do
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

          !$omp parallel do
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
             !$omp parallel do &
             !$omp private(flag)
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
             !$omp parallel do &
             !$omp private(flag)
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

          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
             Zsfc(i,j) = Zsfc(i,j) &
                       + 0.1_RP * ( ( FLX_X(i,j) - FLX_X(i-1,j) ) &
                                  + ( FLX_Y(i,j) - FLX_Y(i,j-1) ) )
          enddo
          enddo

       case default
          LOG_ERROR("CNVTOPO_smooth",*) 'Invalid smoothing type'
          call PRC_abort
       end select

       if ( CNVTOPO_smooth_trim_ocean ) then
          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
             Zsfc(i,j) = sign( max( Zsfc(i,j) * TOPO_sign(i,j), 0.0_RP ), TOPO_sign(i,j) )
          end do
          end do
       end if

    enddo

    if ( ite  > CNVTOPO_smooth_itelim ) then
       LOG_ERROR("CNVTOPO_smooth",*) 'Smoothing did not converge until ', CNVTOPO_smooth_itelim,' times of iteration.'

       LOG_ERROR_CONT(*) 'Please try different parameters of PARAM_CNVTOPO.'
       LOG_ERROR_CONT(*) '- Number limit of iteration                (CNVTOPO_smooth_itelim)         = ', CNVTOPO_smooth_itelim
       LOG_ERROR_CONT(*) '- Maximum ratio of slope dZ/dX, dZ/dY      (CNVTOPO_smooth_maxslope_ratio) = ', CNVTOPO_smooth_maxslope_ratio
       LOG_ERROR_CONT(*) '  Or, Maximum of slope with degree         (CNVTOPO_smooth_maxslope)       = ', CNVTOPO_smooth_maxslope
       LOG_ERROR_CONT(*) '- Smoothing type LAPLACIAN/GAUSSIAN/OFF    (CNVTOPO_smooth_type)           = ', trim(CNVTOPO_smooth_type)
       LOG_ERROR_CONT(*) '- Number of using points for interpolation (CNVTOPO_interp_level)          = ', CNVTOPO_interp_level
       call PRC_abort
    else
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'smoothing complete.'
    endif



    if ( CNVTOPO_smooth_hypdiff_niter > 0 ) then

       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Apply hyperdiffusion.'

       call FILTER_hyperdiff( IA, IS, IE, JA, JS, JE, &
                              Zsfc(:,:), &
                              CNVTOPO_smooth_hypdiff_order, CNVTOPO_smooth_hypdiff_niter, &
                              limiter_sign = TOPO_sign(:,:)                               )

       !$omp parallel do
       do j = 1, JA
       do i = 1, IA-1
          DZsfc_DXY(i,j,1) = atan2( ( Zsfc(i+1,j)-Zsfc(i,j) ), DXL(i) ) / D2R
       enddo
       enddo
       DZsfc_DXY(IA,:,1) = 0.0_RP
       !$omp parallel do
       do j = 1, JA-1
       do i = 1, IA
          DZsfc_DXY(i,j,2) = atan2( ( Zsfc(i,j+1)-Zsfc(i,j) ), DYL(j) ) / D2R
       enddo
       enddo
       DZsfc_DXY(:,JA,2) = 0.0_RP

       slope(:,:) = max( abs(DZsfc_DXY(:,:,1)), abs(DZsfc_DXY(:,:,2)) )
       call STATISTICS_horizontal_max( IA, IS, IE, JA, JS, JE, &
                                       slope(:,:), maxslope )

       LOG_INFO("CNVTOPO_smooth",*) 'maximum slope [deg] : ', maxslope

    end if

    call TOPO_fillhalo( Zsfc=Zsfc(:,:), FILL_BND=.true. )

    call STATISTICS_detail( IA, IS, IE, JA, JS, JE, 2, &
                            varname(:), DZsfc_DXY(:,:,:) )

    LOG_NEWLINE

    return
  end subroutine CNVTOPO_smooth

end module mod_cnvtopo
