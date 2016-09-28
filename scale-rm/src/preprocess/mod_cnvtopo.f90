!-------------------------------------------------------------------------------
!> module Convert topography
!!
!! @par Description
!!          subroutines for preparing topography data (convert from external file)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_cnvtopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
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

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNVTOPO_GTOPO30
  private :: CNVTOPO_GMTED2010
  private :: CNVTOPO_DEM50M
  private :: CNVTOPO_smooth

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: CNVTOPO_smooth_type         = 'LAPLACIAN'
  logical,                private :: CNVTOPO_smooth_local        = .false.
  integer,                private :: CNVTOPO_smooth_itelim       = 10000

  logical,                private :: CNVTOPO_copy_parent         = .false.

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
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       D2R  => CONST_D2R, &
       HUGE => CONST_HUGE
    use scale_comm, only: &
       COMM_horizontal_min
    use scale_grid, only: &
       DX,       &
       DY,       &
       GRID_CDZ, &
       GRID_FDX, &
       GRID_FDY
    use scale_grid_real, only: &
       REAL_DLAT, &
       REAL_DLON
    implicit none

    character(len=H_SHORT) :: CNVTOPO_name = 'NONE' ! keep backward compatibility

    NAMELIST / PARAM_CNVTOPO / &
       CNVTOPO_name,                  &
       CNVTOPO_UseGTOPO30,            &
       CNVTOPO_UseGMTED2010,          &
       CNVTOPO_UseDEM50M,             &
       CNVTOPO_unittile_ddeg,         &
       CNVTOPO_oversampling_factor,   &
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[convert topo] / Categ[preprocess] / Origin[SCALE-RM]'

    if ( CNVTOPO_smooth_local ) then
       DXL(:) = DX
       DYL(:) = DY
    else
       DXL(:) = GRID_FDX(:)
       DYL(:) = GRID_FDY(:)
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CNVTOPO. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CNVTOPO)

    select case(CNVTOPO_name)
    case('NONE')
       ! do nothing
    case('GTOPO30')
       CNVTOPO_UseGTOPO30   = .true.
       CNVTOPO_UseGMTED2010 = .false.
       CNVTOPO_UseDEM50M    = .false.
    case('GMTED2010')
       CNVTOPO_UseGTOPO30   = .false.
       CNVTOPO_UseGMTED2010 = .true.
       CNVTOPO_UseDEM50M    = .false.
    case('DEM50M')
       CNVTOPO_UseGTOPO30   = .false.
       CNVTOPO_UseGMTED2010 = .false.
       CNVTOPO_UseDEM50M    = .true.
    case('COMBINE')
       CNVTOPO_UseGTOPO30   = .true.
       CNVTOPO_UseGMTED2010 = .true.
       CNVTOPO_UseDEM50M    = .true.
    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(CNVTOPO_name)
       call PRC_MPIstop
    endselect

    CNVTOPO_DoNothing = .true.

    if ( CNVTOPO_UseGTOPO30 ) then
       CNVTOPO_DoNothing = .false.
       if( IO_L ) write(IO_FID_LOG,*) '*** Use GTOPO, global 30 arcsec. data'
       if ( CNVTOPO_UseGMTED2010 ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Use GMTED2010, new global 5 arcsec. data'
          if( IO_L ) write(IO_FID_LOG,*) '*** Overwrite Existing region'
       endif
       if ( CNVTOPO_UseDEM50M ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Use DEM 50m data for Japan region'
          if( IO_L ) write(IO_FID_LOG,*) '*** Overwrite Japan region'
       endif
    elseif ( CNVTOPO_UseGMTED2010 ) then
       CNVTOPO_DoNothing = .false.
       if( IO_L ) write(IO_FID_LOG,*) '*** Use GMTED2010, new global 5 arcsec. data'
       if ( CNVTOPO_UseDEM50M ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Use DEM 50m data for Japan region'
          if( IO_L ) write(IO_FID_LOG,*) '*** Overwrite Japan region'
       endif
    elseif ( CNVTOPO_UseDEM50M ) then
       CNVTOPO_DoNothing = .false.
       if( IO_L ) write(IO_FID_LOG,*) '*** Use DEM 50m data, Japan region only'
    endif

    if ( CNVTOPO_DoNothing ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Do nothing for landuse index'
    else
       drad(:,:) = min( REAL_DLAT(:,:), REAL_DLON(:,:) )
       call COMM_horizontal_min( drad_min, drad(:,:) )

       if ( CNVTOPO_unittile_ddeg > 0.0_RP ) then
          CNVTOPO_oversampling_factor = ( drad_min / D2R ) / CNVTOPO_unittile_ddeg
       endif
       CNVTOPO_oversampling_factor = max( 1.0_RP, CNVTOPO_oversampling_factor )
       CNVTOPO_unittile_ddeg       = ( drad_min / D2R ) / CNVTOPO_oversampling_factor

       if( IO_L ) write(IO_FID_LOG,*) '*** The size of tile [deg] = ', CNVTOPO_unittile_ddeg
       if( IO_L ) write(IO_FID_LOG,*) '*** oversampling factor    = ', CNVTOPO_oversampling_factor
    endif

    if( CNVTOPO_smooth_maxslope > 0.0_RP ) then

      CNVTOPO_smooth_maxslope_limit = CNVTOPO_smooth_maxslope

    else
      minslope(:,:) = HUGE

      j = JS-1
      i = IS-1
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DYL(j) ) / D2R
         minslope(IS,JS) = min( minslope(IS,JS), DZDX, DZDY )
      enddo

      j = JS-1
      do i = IS, IE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DYL(j) ) / D2R
         minslope(i,JS) = min( minslope(i,JS), DZDX, DZDY )
      enddo
      enddo

      i = IS-1
      do j = JS, JE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DYL(j) ) / D2R
         minslope(IS,j) = min( minslope(IS,j), DZDX, DZDY )
      enddo
      enddo

      do j = JS, JE
      do i = IS, IE
      do k = KS, KE
         DZDX = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DXL(i) ) / D2R
         DZDY = atan2( CNVTOPO_smooth_maxslope_ratio * GRID_CDZ(k), DYL(j) ) / D2R
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
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI => CONST_PI
    use scale_grid, only: &
       CBFX => GRID_CBFX, &
       CBFY => GRID_CBFY
    use scale_topography, only: &
       TOPO_fillhalo, &
       TOPO_Zsfc, &
       TOPO_write
    use mod_copytopo, only: &
       COPYTOPO
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( CNVTOPO_DoNothing ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  CONVERT TOPOGRAPHY DATA ++++++'
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ START CONVERT TOPOGRAPHY DATA ++++++'

       if ( CNVTOPO_UseGTOPO30 ) then
          call CNVTOPO_GTOPO30
       endif

       if ( CNVTOPO_UseGMTED2010 ) then
          call CNVTOPO_GMTED2010
       endif

       if ( CNVTOPO_UseDEM50M ) then
          call CNVTOPO_DEM50M
       endif

       call CNVTOPO_smooth( TOPO_Zsfc(:,:) ) ! (inout)
       call TOPO_fillhalo

       if( CNVTOPO_copy_parent ) call COPYTOPO( TOPO_Zsfc )

       if( IO_L ) write(IO_FID_LOG,*) '++++++ END   CONVERT TOPOGRAPHY DATA ++++++'

       ! output topography file
       call TOPO_write
    endif

    return
  end subroutine CNVTOPO

  !-----------------------------------------------------------------------------
  !> Convert from GTOPO30
  subroutine CNVTOPO_GTOPO30
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_grid_real, only: &
       REAL_LATY, &
       REAL_LONX
    implicit none

    character(len=H_LONG) :: GTOPO30_IN_CATALOGUE = '' !< metadata files for GTOPO30
    character(len=H_LONG) :: GTOPO30_IN_DIR       = '' !< directory contains GTOPO30 files (GrADS format)

    NAMELIST / PARAM_CNVTOPO_GTOPO30 / &
       GTOPO30_IN_CATALOGUE, &
       GTOPO30_IN_DIR

    ! data catalogue list
    integer, parameter      :: TILE_nlim = 100
    integer                 :: TILE_nmax
    real(RP)                :: TILE_LATS (TILE_nlim)
    real(RP)                :: TILE_LATE (TILE_nlim)
    real(RP)                :: TILE_LONS (TILE_nlim)
    real(RP)                :: TILE_LONE (TILE_nlim)
    character(len=H_LONG)   :: TILE_fname(TILE_nlim)

    ! GTOPO30 data
    integer, parameter      :: isize_orig = 4800
    integer, parameter      :: jsize_orig = 6000
    integer(2)              :: TILE_HEIGHT_orig(isize_orig,jsize_orig)
    real(RP)                :: TILE_DLAT_orig, TILE_DLON_orig

    ! GTOPO30 data (oversampling)
    integer                 :: ios
    integer                 :: jos
    integer                 :: isize
    integer                 :: jsize
    integer(2), allocatable :: TILE_HEIGHT(:,:)
    real(RP),   allocatable :: TILE_LATH  (:)
    real(RP),   allocatable :: TILE_LONH  (:)
    real(RP)                :: TILE_DLAT, TILE_DLON
    real(RP)                :: area, area_fraction

    integer  :: iloc
    integer  :: jloc
    real(RP) :: ifrac_l ! fraction for iloc
    real(RP) :: jfrac_b ! fraction for jloc

    real(RP) :: REAL_LONX_mod(IA,JA)
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[convert GTOPO30] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO_GTOPO30,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CNVTOPO_GTOPO30. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CNVTOPO_GTOPO30)

    do j = 1, JA
    do i = 1, IA
       area_sum(i,j) = 0.0_RP
       topo_sum(i,j) = 0.0_RP
    enddo
    enddo

    REAL_LONX_mod(:,:) = mod( REAL_LONX(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI

    DOMAIN_LATS = minval(REAL_LATY    (:,:))
    DOMAIN_LATE = maxval(REAL_LATY    (:,:))
    DOMAIN_LONS = minval(REAL_LONX_mod(:,:))
    DOMAIN_LONE = maxval(REAL_LONX_mod(:,:))

    DOMAIN_LONSLOC = minloc(REAL_LONX_mod(:,:))
    DOMAIN_LONELOC = maxloc(REAL_LONX_mod(:,:))

    check_IDL = .false.
    if (      DOMAIN_LONS < REAL_LONX_mod(1 ,DOMAIN_LONSLOC(2)) &
         .OR. DOMAIN_LONE > REAL_LONX_mod(IA,DOMAIN_LONELOC(2)) ) then
       check_IDL = .true.
       DOMAIN_LONS = minval(REAL_LONX_mod(:,:),mask=(REAL_LONX_mod>0.0_RP))
       DOMAIN_LONE = maxval(REAL_LONX_mod(:,:),mask=(REAL_LONX_mod<0.0_RP))
    endif

    jos   = nint( 30.0_RP / 60.0_RP / 60.0_RP / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    ios   = nint( 30.0_RP / 60.0_RP / 60.0_RP / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    jsize = jsize_orig * jos
    isize = isize_orig * ios

    allocate( TILE_HEIGHT(isize,jsize) )
    allocate( TILE_LATH  (0:jsize)     )
    allocate( TILE_LONH  (0:isize)     )

    if( IO_L ) write(IO_FID_LOG,*) '*** Oversampling (j) orig = ', jsize_orig, ', use = ', jsize
    if( IO_L ) write(IO_FID_LOG,*) '*** Oversampling (i) orig = ', isize_orig, ', use = ', isize

    TILE_DLAT_orig = 30.0_RP / 60.0_RP / 60.0_RP * D2R
    TILE_DLON_orig = 30.0_RP / 60.0_RP / 60.0_RP * D2R
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLAT       :', TILE_DLAT_orig/D2R
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLON       :', TILE_DLON_orig/D2R

    TILE_DLAT = TILE_DLAT_orig / jos
    TILE_DLON = TILE_DLON_orig / ios
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLAT (OS)  :', TILE_DLAT/D2R
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLON (OS)  :', TILE_DLON/D2R

    !---< READ from external files >---

    ! catalogue file
    fname = trim(GTOPO30_IN_DIR)//'/'//trim(GTOPO30_IN_CATALOGUE)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input catalogue file:', trim(fname)

    fid = IO_get_available_fid()
    open( fid,                  &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          write(*,*) 'xxx catalogue file not found!', trim(fname)
          call PRC_MPIstop
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

          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '+++ Input data file :', trim(fname)
          if( IO_L ) write(IO_FID_LOG,*) '*** Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          if( IO_L ) write(IO_FID_LOG,*) '***        (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if ( check_IDL ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** (Date line exists within the domain)'
          endif
          if( IO_L ) write(IO_FID_LOG,*) '*** Tile   (LAT)    :', TILE_LATS(t), TILE_LATE(t)
          if( IO_L ) write(IO_FID_LOG,*) '***        (LON)    :', TILE_LONS(t), TILE_LONE(t)

          fid = IO_get_available_fid()
          open( fid,                              &
                file   = trim(fname),             &
                form   = 'unformatted',           &
                access = 'direct',                &
                status = 'old',                   &
                recl   = isize_orig*jsize_orig*2, &
                iostat = ierr                     )

             if ( ierr /= 0 ) then
                write(*,*) 'xxx data file not found!'
                call PRC_MPIstop
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
!             if( IO_L ) write(IO_FID_LOG,*) jj, TILE_LATH(jj)
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, isize
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
!             if( IO_L ) write(IO_FID_LOG,*) ii, TILE_LONH(ii)
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
                if (       TILE_LONH(ii-1) >= REAL_LONX_mod(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  REAL_LONX_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= REAL_LATY    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  REAL_LATY    (i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( REAL_LONX_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif

                if (       REAL_LONX_mod(i-1,j) >= REAL_LONX_mod(i  ,j  ) &
                     .AND. TILE_LATH    (jj-1)  >= REAL_LATY    (i  ,j-1) &
                     .AND. TILE_LATH    (jj-1)  <  REAL_LATY    (i  ,j  ) ) then ! across the IDL

                   if    (       TILE_LONH(ii-1) >= REAL_LONX_mod(i-1,j) &
                           .AND. TILE_LONH(ii-1) <  PI                   ) then

                      iloc    = i
                      ifrac_l = min( REAL_LONX_mod(i,j)-TILE_LONH(ii-1)+2.0_RP*PI, TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   elseif(       TILE_LONH(ii-1) >= -PI                  &
                           .AND. TILE_LONH(ii-1) <  REAL_LONX_mod(i  ,j) ) then

                      iloc    = i
                      ifrac_l = min( REAL_LONX_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   endif

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             topo = real( TILE_HEIGHT(ii,jj), kind=RP )
             mask = 0.5_RP - sign( 0.5_RP, topo ) ! if Height is negative, mask = 1

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) ) * ( 1.0_RP - mask )

!             if( IO_L ) write(IO_FID_LOG,*) ii, jj, area, iloc, jloc, ifrac_l, jfrac_b, TILE_HEIGHT(ii,jj)

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
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_grid_real, only: &
       REAL_LATY, &
       REAL_LONX
    implicit none

    character(len=H_LONG) :: DEM50M_IN_CATALOGUE = '' !< metadata files for DEM50M
    character(len=H_LONG) :: DEM50M_IN_DIR       = '' !< directory contains DEM50M files (GrADS format)

    NAMELIST / PARAM_CNVTOPO_DEM50M / &
       DEM50M_IN_CATALOGUE, &
       DEM50M_IN_DIR

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 1000
    integer               :: TILE_nmax
    real(RP)              :: TILE_LATS (TILE_nlim)
    real(RP)              :: TILE_LATE (TILE_nlim)
    real(RP)              :: TILE_LONS (TILE_nlim)
    real(RP)              :: TILE_LONE (TILE_nlim)
    character(len=H_LONG) :: TILE_fname(TILE_nlim)

    ! DEM50M data
    integer, parameter    :: isize_orig = 1600
    integer, parameter    :: jsize_orig = 1600
    real(SP)              :: TILE_HEIGHT_orig(isize_orig,jsize_orig)
    real(RP)              :: TILE_DLAT_orig, TILE_DLON_orig

    ! DEM50M data (oversampling)
    integer               :: ios
    integer               :: jos
    integer               :: isize
    integer               :: jsize
    real(SP), allocatable :: TILE_HEIGHT(:,:)
    real(RP), allocatable :: TILE_LATH  (:)
    real(RP), allocatable :: TILE_LONH  (:)
    real(RP)              :: TILE_DLAT, TILE_DLON
    real(RP)              :: area, area_fraction

    integer  :: iloc
    integer  :: jloc
    real(RP) :: ifrac_l ! fraction for iloc
    real(RP) :: jfrac_b ! fraction for jloc

    real(RP) :: REAL_LONX_mod(IA,JA)
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[convert DEM50M] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO_DEM50M,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CNVTOPO_DEM50M. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CNVTOPO_DEM50M)

    do j = 1, JA
    do i = 1, IA
       area_sum(i,j) = 0.0_RP
       topo_sum(i,j) = 0.0_RP
    enddo
    enddo

    REAL_LONX_mod(:,:) = mod( REAL_LONX(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI

    DOMAIN_LATS = minval(REAL_LATY    (:,:))
    DOMAIN_LATE = maxval(REAL_LATY    (:,:))
    DOMAIN_LONS = minval(REAL_LONX_mod(:,:))
    DOMAIN_LONE = maxval(REAL_LONX_mod(:,:))

    DOMAIN_LONSLOC = minloc(REAL_LONX_mod(:,:))
    DOMAIN_LONELOC = maxloc(REAL_LONX_mod(:,:))

    check_IDL = .false.
    if (      DOMAIN_LONS < REAL_LONX_mod(1 ,DOMAIN_LONSLOC(2)) &
         .OR. DOMAIN_LONE > REAL_LONX_mod(IA,DOMAIN_LONELOC(2)) ) then
       check_IDL = .true.
       DOMAIN_LONS = minval(REAL_LONX_mod(:,:),mask=(REAL_LONX_mod>0.0_RP))
       DOMAIN_LONE = maxval(REAL_LONX_mod(:,:),mask=(REAL_LONX_mod<0.0_RP))
    endif

    jos   = nint( 5.0_RP / 60.0_RP / 200.0_RP / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    ios   = nint( 7.5_RP / 60.0_RP / 200.0_RP / CNVTOPO_unittile_ddeg - 0.5_RP ) + 1
    jsize = jsize_orig * jos
    isize = isize_orig * ios

    allocate( TILE_HEIGHT(isize,jsize) )
    allocate( TILE_LATH  (0:jsize)     )
    allocate( TILE_LONH  (0:isize)     )

    if( IO_L ) write(IO_FID_LOG,*) '*** Oversampling (j) orig = ', jsize_orig, ', use = ', jsize
    if( IO_L ) write(IO_FID_LOG,*) '*** Oversampling (i) orig = ', isize_orig, ', use = ', isize

    TILE_DLAT_orig = 5.0_RP / 60.0_RP / 200.0_RP * D2R
    TILE_DLON_orig = 7.5_RP / 60.0_RP / 200.0_RP * D2R
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLAT       :', TILE_DLAT_orig/D2R
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLON       :', TILE_DLON_orig/D2R

    TILE_DLAT = TILE_DLAT_orig / jos
    TILE_DLON = TILE_DLON_orig / ios
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLAT (OS)  :', TILE_DLAT/D2R
    if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLON (OS)  :', TILE_DLON/D2R

    !---< READ from external files >---

    ! catalogue file
    fname = trim(DEM50M_IN_DIR)//'/'//trim(DEM50M_IN_CATALOGUE)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input catalogue file:', trim(fname)

    fid = IO_get_available_fid()
    open( fid,                  &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          write(*,*) 'xxx catalogue file not found!', trim(fname)
          call PRC_MPIstop
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

          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '+++ Input data file :', trim(fname)
          if( IO_L ) write(IO_FID_LOG,*) '*** Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          if( IO_L ) write(IO_FID_LOG,*) '***        (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if ( check_IDL ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** (Date line exists within the domain)'
          endif
          if( IO_L ) write(IO_FID_LOG,*) '*** Tile   (LAT)    :', TILE_LATS(t), TILE_LATE(t)
          if( IO_L ) write(IO_FID_LOG,*) '***        (LON)    :', TILE_LONS(t), TILE_LONE(t)

          fid = IO_get_available_fid()
          open( fid,                              &
                file   = trim(fname),             &
                form   = 'unformatted',           &
                access = 'direct',                &
                status = 'old',                   &
                recl   = isize_orig*jsize_orig*4, &
                iostat = ierr                     )

             if ( ierr /= 0 ) then
                write(*,*) 'xxx data file not found!'
                call PRC_MPIstop
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
!             if( IO_L ) write(IO_FID_LOG,*) jj, TILE_LATH(jj)
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, isize
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
!             if( IO_L ) write(IO_FID_LOG,*) ii, TILE_LONH(ii)
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
                if (       TILE_LONH(ii-1) >= REAL_LONX_mod(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  REAL_LONX_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= REAL_LATY    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  REAL_LATY    (i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( REAL_LONX_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif

                if (       REAL_LONX_mod(i-1,j) >= REAL_LONX_mod(i  ,j  ) &
                     .AND. TILE_LATH    (jj-1)  >= REAL_LATY    (i  ,j-1) &
                     .AND. TILE_LATH    (jj-1)  <  REAL_LATY    (i  ,j  ) ) then ! across the IDL

                   if    (       TILE_LONH(ii-1) >= REAL_LONX_mod(i-1,j) &
                           .AND. TILE_LONH(ii-1) <  PI                   ) then

                      iloc    = i
                      ifrac_l = min( REAL_LONX_mod(i,j)-TILE_LONH(ii-1)+2.0_RP*PI, TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   elseif(       TILE_LONH(ii-1) >= -PI                  &
                           .AND. TILE_LONH(ii-1) <  REAL_LONX_mod(i  ,j) ) then

                      iloc    = i
                      ifrac_l = min( REAL_LONX_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   endif

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             topo = real( TILE_HEIGHT(ii,jj), kind=RP )
             mask = 0.5_RP - sign( 0.5_RP,topo ) ! if Height is negative, mask = 1

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) ) * ( 1.0_RP - mask )

!             if( IO_L ) write(IO_FID_LOG,*) ii, jj, area, iloc, jloc, ifrac_l, jfrac_b, TILE_HEIGHT(ii,jj)

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
  !> check slope
  subroutine CNVTOPO_smooth( &
       Zsfc )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
       DX,       &
       DY,       &
       GRID_FDX, &
       GRID_FDY
    use scale_comm, only: &
       COMM_horizontal_max
    use scale_rm_statistics, only: &
       STAT_detail
    use scale_topography, only: &
       TOPO_fillhalo
    implicit none

    real(RP), intent(inout) :: Zsfc(IA,JA)

    real(RP) :: DZsfc_DX(1,IA,JA,1) ! d(Zsfc)/dx at u-position
    real(RP) :: DZsfc_DY(1,IA,JA,1) ! d(Zsfc)/dy at v-position

    real(RP) :: DXL(IA-1)
    real(RP) :: DYL(JA-1)

    real(RP) :: FLX_X(IA,JA)
    real(RP) :: FLX_Y(IA,JA)
    real(RP) :: FLX_TMP(IA,JA)

    real(RP) :: slope(IA,JA)
    real(RP) :: maxslope
    real(RP) :: flag

    character(len=H_SHORT) :: varname(1)

    integer :: ite
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Apply smoothing. Slope limit       = ', CNVTOPO_smooth_maxslope_limit
    if( IO_L ) write(IO_FID_LOG,*) '***                  Smoothing type    = ', CNVTOPO_smooth_type
    if( IO_L ) write(IO_FID_LOG,*) '***                  Smoothing locally = ', CNVTOPO_smooth_local

    if ( CNVTOPO_smooth_local ) then
       DXL(:) = DX
       DYL(:) = DY
    else
       DXL(:) = GRID_FDX(:)
       DYL(:) = GRID_FDY(:)
    endif

    ! digital filter
    do ite = 1, CNVTOPO_smooth_itelim+1
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Smoothing itelation : ', ite

       call TOPO_fillhalo( Zsfc )

       do j = 1, JA
       do i = 1, IA-1
          DZsfc_DX(1,i,j,1) = atan2( ( Zsfc(i+1,j)-Zsfc(i,j) ), DXL(i) ) / D2R
       enddo
       enddo
       DZsfc_DX(1,IA,:,1) = 0.0_RP
       do j = 1, JA-1
       do i = 1, IA
          DZsfc_DY(1,i,j,1) = atan2( ( Zsfc(i,j+1)-Zsfc(i,j) ), DYL(j) ) / D2R
       enddo
       enddo
       DZsfc_DY(1,:,JA,1) = 0.0_RP

       slope(:,:) = max( abs(DZsfc_DX(1,:,:,1)), abs(DZsfc_DY(1,:,:,1)) )
       call COMM_horizontal_max( maxslope, slope(:,:) )

       if( IO_L ) write(IO_FID_LOG,*) '*** maximum slope [deg] : ', maxslope

       if( maxslope < CNVTOPO_smooth_maxslope_limit ) exit

       varname(1) = "DZsfc_DX"
       call STAT_detail( DZsfc_DX(:,:,:,:), varname(:) )
       varname(1) = "DZsfc_DY"
       call STAT_detail( DZsfc_DY(:,:,:,:), varname(:) )

       select case ( CNVTOPO_smooth_type )
       case ( 'GAUSSIAN' )

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

       case ( 'LAPLACIAN' )

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
                     + sign(0.5_RP, max( abs(DZsfc_DX(1,i+1,j  ,1)), &
                                         abs(DZsfc_DX(1,i  ,j  ,1)), &
                                         abs(DZsfc_DX(1,i-1,j  ,1)), &
                                         abs(DZsfc_DY(1,i+1,j  ,1)), &
                                         abs(DZsfc_DY(1,i+1,j-1,1)), &
                                         abs(DZsfc_DY(1,i  ,j  ,1)), &
                                         abs(DZsfc_DY(1,i  ,j-1,1))  &
                                       ) - CNVTOPO_smooth_maxslope_limit )
                FLX_X(i,j) = FLX_X(i,j) &
                           * DXL(i) / GRID_FDX(i) &
                           * flag
             enddo
             enddo
             do j = JS-1, JE
             do i = IS  , IE
                flag = 0.5_RP &
                     + sign(0.5_RP, max( abs(DZsfc_DY(1,i  ,j+1,1)), &
                                         abs(DZsfc_DY(1,i  ,j  ,1)), &
                                         abs(DZsfc_DY(1,i  ,j-1,1)), &
                                         abs(DZsfc_DX(1,i  ,j+1,1)), &
                                         abs(DZsfc_DX(1,i-1,j+1,1)), &
                                         abs(DZsfc_DX(1,i  ,j  ,1)), &
                                         abs(DZsfc_DX(1,i-1,j  ,1))  &
                                       ) - CNVTOPO_smooth_maxslope_limit )
                FLX_Y(i,j) = FLX_Y(i,j) &
                           * DYL(j) / GRID_FDY(j) &
                           * flag
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
          write(*,*) 'xxx Invalid smoothing type'
          call PRC_MPIstop
       end select

    enddo

    if ( ite  > CNVTOPO_smooth_itelim ) then
       write(*,*) 'xxx not converged'
       call PRC_MPIstop
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** smoothing complete.'
    endif

    varname(1) = "DZsfc_DX"
    call STAT_detail( DZsfc_DX(:,:,:,:), varname(:) )
    varname(1) = "DZsfc_DY"
    call STAT_detail( DZsfc_DY(:,:,:,:), varname(:) )

    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine CNVTOPO_smooth

end module mod_cnvtopo
