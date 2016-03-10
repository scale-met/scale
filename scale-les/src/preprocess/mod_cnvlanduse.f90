!-------------------------------------------------------------------------------
!> module Convert LandUseIndex
!!
!! @par Description
!!          subroutines for preparing landuse index data (convert from external file)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_cnvlanduse
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
  public :: CNVLANDUSE_setup
  public :: CNVLANDUSE

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public            :: CNVLANDUSE_TYPE = -1

  integer, public, parameter :: I_IGNORE     =  0
  integer, public, parameter :: I_GLCCv2       =  1
  integer, public, parameter :: I_LU100M     =  2

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNVLANDUSE_GLCCv2
  private :: CNVLANDUSE_LU100M
  private :: CNVLANDUSE_cal_fact

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNVLANDUSE_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT) :: CNVLANDUSE_name = 'NONE'

    NAMELIST / PARAM_CNVLANDUSE / &
       CNVLANDUSE_name

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CNVLANDUSE]/Categ[CNVLANDUSE]'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVLANDUSE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CNVLANDUSE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CNVLANDUSE)

    select case(CNVLANDUSE_name)
    case('NONE')
       CNVLANDUSE_TYPE = I_IGNORE

    case('GLCCv2')
       CNVLANDUSE_TYPE = I_GLCCv2

    case('LU100M')
       CNVLANDUSE_TYPE = I_LU100M

    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(CNVLANDUSE_name)
       call PRC_MPIstop
    endselect

    return
  end subroutine CNVLANDUSE_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVLANDUSE
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_write
    implicit none
    !---------------------------------------------------------------------------

    if ( CNVLANDUSE_TYPE == I_IGNORE ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  CONVERT LANDUSE DATA ++++++'
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ START CONVERT LANDUSE DATA ++++++'

       select case(CNVLANDUSE_TYPE)
       case(I_GLCCv2)
           call CNVLANDUSE_GLCCv2

       case(I_LU100M)
           call CNVLANDUSE_LU100M

       case default
          write(*,*) ' xxx Unsupported TYPE:', CNVLANDUSE_TYPE
          call PRC_MPIstop
       endselect

       if( IO_L ) write(IO_FID_LOG,*) '++++++ END   CONVERT LANDUSE DATA ++++++'

       ! output landuse file
       call LANDUSE_write
    endif

    return
  end subroutine CNVLANDUSE

  !-----------------------------------------------------------------------------
  !> Convert from GLCCv2
  subroutine CNVLANDUSE_GLCCv2
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_landuse, only: &
       LANDUSE_frac_land,  &
       LANDUSE_frac_lake,  &
       LANDUSE_frac_urban, &
       LANDUSE_PFT_mosaic, &
       LANDUSE_PFT_nmax,   &
       LANDUSE_frac_PFT,   &
       LANDUSE_index_PFT
    use scale_grid_real, only: &
       REAL_LATY, &
       REAL_LONX
    implicit none

    character(len=H_LONG) :: GLCCv2_IN_CATALOGUE = '' !< metadata files for GLCCv2
    character(len=H_LONG) :: GLCCv2_IN_DIR       = '' !< directory contains GLCCv2 files (GrADS format)
    real(RP)              :: limit_urban_fraction = 1.0_RP !< fraction limiter for urban area

    NAMELIST / PARAM_CNVLANDUSE_GLCCv2 / &
       GLCCv2_IN_CATALOGUE, &
       GLCCv2_IN_DIR,       &
       limit_urban_fraction

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 100
    integer               :: TILE_nmax
    real(RP)              :: TILE_LATS (TILE_nlim)
    real(RP)              :: TILE_LATE (TILE_nlim)
    real(RP)              :: TILE_LONS (TILE_nlim)
    real(RP)              :: TILE_LONE (TILE_nlim)
    character(len=H_LONG) :: TILE_fname(TILE_nlim)

    ! GLCCv2 data
    integer(1) :: TILE_LANDUSE(3600,3600)
    real(RP) :: TILE_LATH   (0:3600)
    real(RP) :: TILE_LONH   (0:3600)
    real(RP) :: TILE_DLAT, TILE_DLON
    real(RP) :: area, area_fraction

    integer  :: iloc
    integer  :: jloc
    real(RP) :: ifrac_l ! fraction for iloc
    real(RP) :: jfrac_b ! fraction for jloc


    real(RP) :: DOMAIN_LATS, DOMAIN_LATE
    real(RP) :: DOMAIN_LONS, DOMAIN_LONE

    ! tentative: LANDUSE_PFT_nmax is assumed to be 4
    real(RP) :: categ_sum(IA,JA,-2:LANDUSE_PFT_nmax)

    integer  :: lookuptable(1:25)
    data lookuptable /  0, & !   1 Urban and Built-Up Land        ->  0 urban
                        7, & !   2 Dryland Cropland and Pasture   ->  7 Dryland Cropland and Pasture
                        8, & !   3 Irrigated Cropland and Pasture ->  8 Irrigated Cropland and Pasture
                        9, & !   4 Mixed Cropland and Pasture     ->  9 Mixed Cropland and Pasture
                        5, & !   5 Cropland/Grassland Mosaic      ->  5 Cropland/Grassland Mosaic
                        6, & !   6 Cropland/Woodland Mosaic       ->  6 Cropland/Woodland Mosaic
                        2, & !   7 Grassland                      ->  2 Grassland
                        3, & !   8 Shrubland                      ->  3 Shrubland
                        4, & !   9 Mixed Shrubland/Grassland      ->  4 Mixed Shrubland/Grassland
                        4, & !  10 Savanna                        ->  4 Mixed Shrubland/Grassland
                       11, & !  11 Deciduous Broadleaf Forest     -> 11 Deciduous Broadleaf Forest
                       12, & !  12 Deciduous Needleleaf Forest    -> 12 Deciduous Needleleaf Forest
                       13, & !  13 Evergreen Broadleaf Forest     -> 13 Deciduous Broadleaf Forest
                       14, & !  14 Evergreen Needleleaf Forest    -> 14 Deciduous Needleleaf Forest
                       15, & !  15 Mixed Forest                   -> 15 Mixed Forest
                       -2, & !  16 Water Bodies                   -> -2 Lake/River
                       10, & !  17 Herbaceous Wetland             -> 10 Paddy
                       10, & !  18 Wooded Wetland                 -> 10 Paddy
                        1, & !  19 Barren or Sparsely Vegetated   ->  1 Bare Ground
                        1, & !  20 Herbaceous Tundra              ->  1 Bare Ground
                        1, & !  21 Wooded Tundra                  ->  1 Bare Ground
                        1, & !  22 Mixed Tundra                   ->  1 Bare Ground
                        1, & !  23 Bare Ground Tundra             ->  1 Bare Ground
                        1, & !  24 Snow or Ice                    ->  1 Bare Ground
                       -1  / !  25 Missing Data                   -> -1 Sea Surface

    real(RP) :: categ_pftsum, allsum
    real(RP) :: PFT    (LANDUSE_PFT_nmax)
    integer  :: PFT_idx(LANDUSE_PFT_nmax)
    real(RP) :: temp
    integer  :: temp_idx

    character(len=H_LONG) :: fname

    real(RP) :: zerosw
    logical  :: hit_lat, hit_lon
    integer  :: index
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, t, p, pp
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[GLCCv2]/Categ[CNVLANDUSE]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVLANDUSE_GLCCv2,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CNVLANDUSE_GLCCv2. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CNVLANDUSE_GLCCv2)

    do p = -2, LANDUSE_PFT_nmax
    do j = 1, JA
    do i = 1, IA
       categ_sum(i,j,p) = 0.0_RP
    enddo
    enddo
    enddo

    DOMAIN_LATS = minval(REAL_LATY(:,:))
    DOMAIN_LATE = maxval(REAL_LATY(:,:))
    DOMAIN_LONS = minval(REAL_LONX(:,:))
    DOMAIN_LONE = maxval(REAL_LONX(:,:))

    !---< READ from external files >---

    fname = trim(GLCCv2_IN_DIR)//'/'//trim(GLCCv2_IN_CATALOGUE)

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
       enddo

       TILE_nmax = t - 1
    close(fid)

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

       if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
            .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
          hit_lon = .true.
       endif

       if (      ( DOMAIN_LONS  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS  < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE  < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(GLCCv2_IN_DIR)//'/'//trim(TILE_fname(t))

          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '+++ Input data      file :', trim(fname)
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** Domain(LAT):', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          if( IO_L ) write(IO_FID_LOG,*) '***       (LON):', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if( IO_L ) write(IO_FID_LOG,*) '*** Tile  (LAT):', TILE_LATS(t), TILE_LATE(t)
          if( IO_L ) write(IO_FID_LOG,*) '***       (LON):', TILE_LONS(t), TILE_LONE(t)

          TILE_DLAT = 30.0_RP / 60.0_RP / 60.0_RP * D2R
          TILE_DLON = 30.0_RP / 60.0_RP / 60.0_RP * D2R
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLAT  :', TILE_DLAT/D2R
          if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLON  :', TILE_DLON/D2R

          fid = IO_get_available_fid()
          open( fid,                    &
                file   = trim(fname),   &
                form   = 'unformatted', &
                access = 'direct',      &
                status = 'old',         &
                recl   = 3600*3600*1,   &
                iostat = ierr           )

             if ( ierr /= 0 ) then
                write(*,*) 'xxx data file not found!'
                call PRC_MPIstop
             endif

             read(fid,rec=1) TILE_LANDUSE(:,:)
          close(fid)

          TILE_LATH(0) = TILE_LATS(t) * D2R
          do jj = 1, 3600
             TILE_LATH(jj) = TILE_LATH(jj-1) + TILE_DLAT
!             if( IO_L ) write(IO_FID_LOG,*) jj, TILE_LATH(jj)
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, 3600
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
!             if( IO_L ) write(IO_FID_LOG,*) ii, TILE_LONH(ii)
          enddo

          ! match and calc fraction
          do jj = 1, 3600
          do ii = 1, 3600

             iloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             ifrac_l = 1.0_RP

             jloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             jfrac_b = 1.0_RP

             if (      TILE_LONH(ii-1) <  DOMAIN_LONS &
                  .OR. TILE_LONH(ii-1) >= DOMAIN_LONE &
                  .OR. TILE_LATH(jj-1) <  DOMAIN_LATS &
                  .OR. TILE_LATH(jj-1) >= DOMAIN_LATE ) then
                 cycle
             endif

      jloop: do j = JS-1, JE+1
      iloop: do i = IS-1, IE+1
                if (       TILE_LONH(ii-1) >= REAL_LONX(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  REAL_LONX(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= REAL_LATY(i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  REAL_LATY(i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( REAL_LONX(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) )

             pp = min( max( int(TILE_LANDUSE(ii,jj),kind=4), 0 ), 25 )
             p  = lookuptable(pp)

             !if( IO_L ) write(IO_FID_LOG,*) ii, jj, iloc, jloc, p, pp

             area_fraction = (       ifrac_l) * (       jfrac_b) * area
             categ_sum(iloc  ,jloc  ,p) = categ_sum(iloc  ,jloc  ,p) + area_fraction

             area_fraction = (1.0_RP-ifrac_l) * (       jfrac_b) * area
             categ_sum(iloc+1,jloc  ,p) = categ_sum(iloc+1,jloc  ,p) + area_fraction

             area_fraction = (       ifrac_l) * (1.0_RP-jfrac_b) * area
             categ_sum(iloc  ,jloc+1,p) = categ_sum(iloc  ,jloc+1,p) + area_fraction

             area_fraction = (1.0_RP-ifrac_l) * (1.0_RP-jfrac_b) * area
             categ_sum(iloc+1,jloc+1,p) = categ_sum(iloc+1,jloc+1,p) + area_fraction
          enddo
          enddo

       endif
    enddo ! tile loop

    do j = JS, JE
    do i = IS, IE
!       area = RADIUS * RADIUS * (REAL_LONX(i,j)-REAL_LONX(i-1,j)) * ( sin(REAL_LATY(i,j))-sin(REAL_LATY(i,j-1)) )
!       allsum = categ_sum(i,j,-2) + categ_sum(i,j,-1) + categ_sum(i,j,0) + categ_pftsum
!       if ( abs(allsum/area-1.0_RP) > EPS ) then
!          if( IO_L ) write(IO_FID_LOG,*) i,j,allsum/area
!       endif

       do p = 1, LANDUSE_PFT_nmax
          PFT    (p) = categ_sum(i,j,p)
          PFT_idx(p) = p
       enddo

       ! bubble sort
       do p  = 1,   LANDUSE_PFT_nmax-1
       do pp = p+1, LANDUSE_PFT_nmax
          if ( PFT(pp) > PFT(p) ) then
             temp    = PFT(p)
             PFT(p)  = PFT(pp)
             PFT(pp) = temp

             temp_idx    = PFT_idx(p)
             PFT_idx(p)  = PFT_idx(pp)
             PFT_idx(pp) = temp_idx
          endif
       enddo
       enddo

       categ_pftsum = sum( PFT(:) )

       ! land fraction : 1 - ocean / total
       allsum = categ_sum(i,j,-2) + categ_sum(i,j,-1) + categ_sum(i,j,0) + categ_pftsum
       zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
       LANDUSE_frac_land (i,j) = 1.0_RP-zerosw - categ_sum(i,j,-1) * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

       ! lake fraction : lake / ( total - ocean )
       allsum = categ_sum(i,j,-2) + categ_sum(i,j,0) + categ_pftsum
       zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
       LANDUSE_frac_lake (i,j) = categ_sum(i,j,-2) * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

       ! urban fraction : urban / ( total - ocean - lake )
       allsum = categ_sum(i,j,0) + categ_pftsum
       zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
       LANDUSE_frac_urban(i,j) = categ_sum(i,j,0) * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

       ! PFT fraction : PFT / sum( PFT(1:mosaic) )
       allsum = sum( PFT(1:LANDUSE_PFT_mosaic) )
       if ( abs(allsum) > EPS ) then
          do p = 1, LANDUSE_PFT_mosaic
             LANDUSE_frac_PFT (i,j,p) = PFT    (p) / allsum
             LANDUSE_index_PFT(i,j,p) = PFT_idx(p)
          enddo
          ! if no second PFT, set to same as PFT1
          if ( abs(LANDUSE_frac_PFT(i,j,1)-1.0_RP) <= EPS ) then
             LANDUSE_frac_PFT (i,j,:) = 0.0_RP
             LANDUSE_frac_PFT (i,j,1) = 1.0_RP
             LANDUSE_index_PFT(i,j,:) = PFT_idx(1)
          endif
       else ! if no PFT, set to bare ground
          LANDUSE_frac_PFT (i,j,:) = 0.0_RP
          LANDUSE_frac_PFT (i,j,1) = 1.0_RP
          LANDUSE_index_PFT(i,j,:) = 1
       endif

    enddo
    enddo

    if ( limit_urban_fraction < 1.0_RP ) then
       do j = JS, JE
       do i = IS, IE
          if ( LANDUSE_frac_urban(i,j) == 1.0_RP ) then ! if no PFT, set to grassland
             LANDUSE_frac_PFT (i,j,:) = 0.0_RP
             LANDUSE_frac_PFT (i,j,1) = 1.0_RP
             LANDUSE_index_PFT(i,j,:) = 2
          endif
          LANDUSE_frac_urban(i,j) = min( LANDUSE_frac_urban(i,j), limit_urban_fraction )
       enddo
       enddo
    endif

    ! calculate landuse factors
    call CNVLANDUSE_cal_fact

    return
  end subroutine CNVLANDUSE_GLCCv2

  !-----------------------------------------------------------------------------
  !> Convert from KSJ landuse 100m mesh
  subroutine CNVLANDUSE_LU100M
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_landuse, only: &
       LANDUSE_frac_land,  &
       LANDUSE_frac_lake,  &
       LANDUSE_frac_urban, &
       LANDUSE_PFT_mosaic, &
       LANDUSE_PFT_nmax,   &
       LANDUSE_frac_PFT,   &
       LANDUSE_index_PFT
    use scale_grid_real, only: &
       REAL_LATY, &
       REAL_LONX
    implicit none

    character(len=H_LONG) :: LU100M_IN_CATALOGUE = '' !< metadata files for LU100M
    character(len=H_LONG) :: LU100M_IN_DIR       = '' !< directory contains LU100M files (GrADS format)
    real(RP)              :: limit_urban_fraction = 1.0_RP !< fraction limiter for urban area

    NAMELIST / PARAM_CNVLANDUSE_LU100M / &
       LU100M_IN_CATALOGUE, &
       LU100M_IN_DIR,       &
       limit_urban_fraction

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 1000
    integer               :: TILE_nmax
    real(RP)              :: TILE_LATS (TILE_nlim)
    real(RP)              :: TILE_LATE (TILE_nlim)
    real(RP)              :: TILE_LONS (TILE_nlim)
    real(RP)              :: TILE_LONE (TILE_nlim)
    character(len=H_LONG) :: TILE_fname(TILE_nlim)

    ! LU100M data
    real(SP) :: TILE_LANDUSE(800,800)
    real(RP) :: TILE_LATH   (0:800)
    real(RP) :: TILE_LONH   (0:800)
    real(RP) :: TILE_DLAT, TILE_DLON
    real(RP) :: area, area_fraction

    integer  :: iloc
    integer  :: jloc
    real(RP) :: ifrac_l ! fraction for iloc
    real(RP) :: jfrac_b ! fraction for jloc


    real(RP) :: DOMAIN_LATS, DOMAIN_LATE
    real(RP) :: DOMAIN_LONS, DOMAIN_LONE

    ! tentative: LANDUSE_PFT_nmax is assumed to be 4
    real(RP) :: categ_sum(IA,JA,-2:LANDUSE_PFT_nmax)

    integer  :: lookuptable(0:16)
    data lookuptable / -1, & !  0 missing        -> -1 Sea Surface
                       10, & !  1 paddy          -> 10 Paddy
                        9, & !  2 cropland       ->  9 Mixed Cropland and Pasture
                        1, & !  3 UNDEF          ->  1 Bare Ground
                        1, & !  4 UNDEF          ->  1 Bare Ground
                       11, & !  5 forest         -> 11 Deciduous Broadleaf Forest
                        1, & !  6 bareground     ->  1 Bare Ground
                        0, & !  7 urban building ->  0 Urban and Built-up Land
                        1, & !  8 UNDEF          ->  1 Bare Ground
                        0, & !  9 motorway       ->  0 Urban and Built-up Land
                        0, & ! 10 urban ground   ->  0 Urban and Built-up Land
                       -2, & ! 11 lake,river     -> -2 Lake/River
                        1, & ! 12 UNDEF          ->  1 Bare Ground
                        1, & ! 13 UNDEF          ->  1 Bare Ground
                        1, & ! 14 seashore       ->  1 Bare Ground
                       -1, & ! 15 ocean          -> -1 Sea Surface
                        2  / ! 16 golf course    ->  2 Grassland

    real(RP) :: categ_pftsum, allsum
    real(RP) :: PFT    (LANDUSE_PFT_nmax)
    integer  :: PFT_idx(LANDUSE_PFT_nmax)
    real(RP) :: temp
    integer  :: temp_idx

    character(len=H_LONG) :: fname

    real(RP) :: zerosw
    logical  :: hit_lat, hit_lon
    integer  :: index
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, t, p, pp
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LU100M]/Categ[CNVLANDUSE]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVLANDUSE_LU100M,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CNVLANDUSE_LU100M. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CNVLANDUSE_LU100M)

    do p = -2, LANDUSE_PFT_nmax
    do j = 1, JA
    do i = 1, IA
       categ_sum(i,j,p) = 0.0_RP
    enddo
    enddo
    enddo

    DOMAIN_LATS = minval(REAL_LATY(:,:))
    DOMAIN_LATE = maxval(REAL_LATY(:,:))
    DOMAIN_LONS = minval(REAL_LONX(:,:))
    DOMAIN_LONE = maxval(REAL_LONX(:,:))

    !---< READ from external files >---

    fname = trim(LU100M_IN_DIR)//'/'//trim(LU100M_IN_CATALOGUE)

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
       enddo

       TILE_nmax = t - 1
    close(fid)

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

       if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONS(t)*D2R < DOMAIN_LONE  ) &
            .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS  .AND. TILE_LONE(t)*D2R < DOMAIN_LONE  ) ) then
          hit_lon = .true.
       endif

       if (      ( DOMAIN_LONS  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS  < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE  >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE  < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(LU100M_IN_DIR)//'/'//trim(TILE_fname(t))

          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '+++ Input data      file :', trim(fname)
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** Domain(LAT):', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          if( IO_L ) write(IO_FID_LOG,*) '***       (LON):', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if( IO_L ) write(IO_FID_LOG,*) '*** Tile  (LAT):', TILE_LATS(t), TILE_LATE(t)
          if( IO_L ) write(IO_FID_LOG,*) '***       (LON):', TILE_LONS(t), TILE_LONE(t)

          TILE_DLAT = 5.0_RP / 60.0_RP / 100.0_RP * D2R
          TILE_DLON = 7.5_RP / 60.0_RP / 100.0_RP * D2R
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLAT  :', TILE_DLAT/D2R
          if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLON  :', TILE_DLON/D2R

          fid = IO_get_available_fid()
          open( fid,                    &
                file   = trim(fname),   &
                form   = 'unformatted', &
                access = 'direct',      &
                status = 'old',         &
                recl   = 800*800*4,     &
                iostat = ierr           )

             if ( ierr /= 0 ) then
                write(*,*) 'xxx data file not found!'
                call PRC_MPIstop
             endif

             read(fid,rec=1) TILE_LANDUSE(:,:)
          close(fid)

          TILE_LATH(0) = TILE_LATS(t) * D2R
          do jj = 1, 800
             TILE_LATH(jj) = TILE_LATH(jj-1) + TILE_DLAT
!             if( IO_L ) write(IO_FID_LOG,*) jj, TILE_LATH(jj)
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, 800
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
!             if( IO_L ) write(IO_FID_LOG,*) ii, TILE_LONH(ii)
          enddo

          ! match and calc fraction
          do jj = 1, 800
          do ii = 1, 800

             iloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             ifrac_l = 1.0_RP

             jloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             jfrac_b = 1.0_RP

             if (      TILE_LONH(ii-1) <  DOMAIN_LONS &
                  .OR. TILE_LONH(ii-1) >= DOMAIN_LONE &
                  .OR. TILE_LATH(jj-1) <  DOMAIN_LATS &
                  .OR. TILE_LATH(jj-1) >= DOMAIN_LATE ) then
                 cycle
             endif

      jloop: do j = JS-1, JE+1
      iloop: do i = IS-1, IE+1
                if (       TILE_LONH(ii-1) >= REAL_LONX(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  REAL_LONX(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= REAL_LATY(i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  REAL_LATY(i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( REAL_LONX(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) )

             pp = int( max( TILE_LANDUSE(ii,jj), 0.0 ) )
             p  = lookuptable(pp)

             !if( IO_L ) write(IO_FID_LOG,*) ii, jj, iloc, jloc, p, pp

             area_fraction = (       ifrac_l) * (       jfrac_b) * area
             categ_sum(iloc  ,jloc  ,p) = categ_sum(iloc  ,jloc  ,p) + area_fraction

             area_fraction = (1.0_RP-ifrac_l) * (       jfrac_b) * area
             categ_sum(iloc+1,jloc  ,p) = categ_sum(iloc+1,jloc  ,p) + area_fraction

             area_fraction = (       ifrac_l) * (1.0_RP-jfrac_b) * area
             categ_sum(iloc  ,jloc+1,p) = categ_sum(iloc  ,jloc+1,p) + area_fraction

             area_fraction = (1.0_RP-ifrac_l) * (1.0_RP-jfrac_b) * area
             categ_sum(iloc+1,jloc+1,p) = categ_sum(iloc+1,jloc+1,p) + area_fraction
          enddo
          enddo

       endif
    enddo ! tile loop

    do j = JS, JE
    do i = IS, IE
!       area = RADIUS * RADIUS * (REAL_LONX(i,j)-REAL_LONX(i-1,j)) * ( sin(REAL_LATY(i,j))-sin(REAL_LATY(i,j-1)) )
!       allsum = categ_sum(i,j,-2) + categ_sum(i,j,-1) + categ_sum(i,j,0) + categ_pftsum
!       if ( abs(allsum/area-1.0_RP) > EPS ) then
!          if( IO_L ) write(IO_FID_LOG,*) i,j,allsum/area
!       endif

       do p = 1, LANDUSE_PFT_nmax
          PFT    (p) = categ_sum(i,j,p)
          PFT_idx(p) = p
       enddo

       ! bubble sort
       do p  = 1,   LANDUSE_PFT_nmax-1
       do pp = p+1, LANDUSE_PFT_nmax
          if ( PFT(pp) > PFT(p) ) then
             temp    = PFT(p)
             PFT(p)  = PFT(pp)
             PFT(pp) = temp

             temp_idx    = PFT_idx(p)
             PFT_idx(p)  = PFT_idx(pp)
             PFT_idx(pp) = temp_idx
          endif
       enddo
       enddo

       categ_pftsum = sum( PFT(:) )

       ! land fraction : 1 - ocean / total
       allsum = categ_sum(i,j,-2) + categ_sum(i,j,-1) + categ_sum(i,j,0) + categ_pftsum
       zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
       LANDUSE_frac_land (i,j) = 1.0_RP-zerosw - categ_sum(i,j,-1) * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

       ! lake fraction : lake / ( total - ocean )
       allsum = categ_sum(i,j,-2) + categ_sum(i,j,0) + categ_pftsum
       zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
       LANDUSE_frac_lake (i,j) = categ_sum(i,j,-2) * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

       ! urban fraction : urban / ( total - ocean - lake )
       allsum = categ_sum(i,j,0) + categ_pftsum
       zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
       LANDUSE_frac_urban(i,j) = categ_sum(i,j,0) * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

       ! PFT fraction : PFT / sum( PFT(1:mosaic) )
       allsum = sum( PFT(1:LANDUSE_PFT_mosaic) )
       if ( abs(allsum) > EPS ) then
          do p = 1, LANDUSE_PFT_mosaic
             LANDUSE_frac_PFT (i,j,p) = PFT    (p) / allsum
             LANDUSE_index_PFT(i,j,p) = PFT_idx(p)
          enddo
          ! if no second PFT, set to same as PFT1
          if ( abs(LANDUSE_frac_PFT(i,j,1)-1.0_RP) <= EPS ) then
             LANDUSE_frac_PFT (i,j,:) = 0.0_RP
             LANDUSE_frac_PFT (i,j,1) = 1.0_RP
             LANDUSE_index_PFT(i,j,:) = PFT_idx(1)
          endif
       else ! if no PFT, set to bare ground
          LANDUSE_frac_PFT (i,j,:) = 0.0_RP
          LANDUSE_frac_PFT (i,j,1) = 1.0_RP
          LANDUSE_index_PFT(i,j,:) = 1
       endif

    enddo
    enddo

    if ( limit_urban_fraction < 1.0_RP ) then
       do j = JS, JE
       do i = IS, IE
          if ( LANDUSE_frac_urban(i,j) == 1.0_RP ) then ! if no PFT, set to grassland
             LANDUSE_frac_PFT (i,j,:) = 0.0_RP
             LANDUSE_frac_PFT (i,j,1) = 1.0_RP
             LANDUSE_index_PFT(i,j,:) = 2
          endif
          LANDUSE_frac_urban(i,j) = min( LANDUSE_frac_urban(i,j), limit_urban_fraction )
       enddo
       enddo
    endif

    ! calculate landuse factors
    call CNVLANDUSE_cal_fact

    return
  end subroutine CNVLANDUSE_LU100M

  !-----------------------------------------------------------------------------
  !> calculate main factor
  subroutine CNVLANDUSE_cal_fact
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_frac_land,  &
       LANDUSE_frac_lake,  &
       LANDUSE_frac_urban, &
       LANDUSE_fact_ocean, &
       LANDUSE_fact_land,  &
       LANDUSE_fact_urban
    implicit none

    integer  :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ calculate landuse factor'

    ! tentative treatment for lake fraction
    LANDUSE_frac_land(:,:) = LANDUSE_frac_land(:,:) * ( 1.0_RP - LANDUSE_frac_lake(:,:) )

    ! make factors
    LANDUSE_fact_ocean(:,:) = ( 1.0_RP - LANDUSE_frac_land(:,:) )
    LANDUSE_fact_land (:,:) = (          LANDUSE_frac_land(:,:) ) * ( 1.0_RP - LANDUSE_frac_urban(:,:) )
    LANDUSE_fact_urban(:,:) = (          LANDUSE_frac_land(:,:) ) * (          LANDUSE_frac_urban(:,:) )

    return
  end subroutine CNVLANDUSE_cal_fact

end module mod_cnvlanduse
