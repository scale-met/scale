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
  integer, public            :: CNVTOPO_TYPE = -1

  integer, public, parameter :: I_IGNORE     =  0
  integer, public, parameter :: I_GTOPO30    =  1
  integer, public, parameter :: I_DEM50M     =  2
  integer, public, parameter :: I_GMTED2010  =  3

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNVTOPO_GTOPO30
  private :: CNVTOPO_DEM50M
  private :: CNVTOPO_GMTED2010
  private :: CNVTOPO_smooth

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: CNVTOPO_smooth_maxslope

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNVTOPO_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=H_SHORT) :: CNVTOPO_name = 'NONE'

    NAMELIST / PARAM_CNVTOPO / &
       CNVTOPO_name,            &
       CNVTOPO_smooth_maxslope

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CNVTOPO]/Categ[CNVTOPO]'

    CNVTOPO_smooth_maxslope = UNDEF

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
       CNVTOPO_TYPE = I_IGNORE

    case('GTOPO30')
       CNVTOPO_TYPE = I_GTOPO30

    case('DEM50M')
       CNVTOPO_TYPE = I_DEM50M

    case('GMTED2010')
       CNVTOPO_TYPE = I_GMTED2010

    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(CNVTOPO_name)
       call PRC_MPIstop
    endselect

    return
  end subroutine CNVTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVTOPO
    use scale_process, only: &
       PRC_MPIstop
    use scale_topography, only: &
       TOPO_write
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( CNVTOPO_TYPE == I_IGNORE ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  CONVERT TOPOGRAPHY DATA ++++++'
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ START CONVERT TOPOGRAPHY DATA ++++++'

       select case(CNVTOPO_TYPE)
       case(I_GTOPO30)
          call CNVTOPO_GTOPO30

       case(I_DEM50M)
          call CNVTOPO_DEM50M

       case(I_GMTED2010)
          call CNVTOPO_GMTED2010

       case default
          write(*,*) ' xxx Unsupported TYPE:', CNVTOPO_TYPE
          call PRC_MPIstop
       endselect

       call CNVTOPO_smooth

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
    use scale_topography, only: &
       TOPO_Zsfc
    implicit none

    character(len=H_LONG) :: TOPO_GTOPO30_IN_CATALOGUE = ''      !< metadata files for GTOPO30
    character(len=H_LONG) :: TOPO_GTOPO30_IN_DIR       = ''      !< directory contains GTOPO30 files (GrADS format)

    NAMELIST / PARAM_CNVTOPO_GTOPO30 / &
       TOPO_GTOPO30_IN_CATALOGUE, &
       TOPO_GTOPO30_IN_DIR

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[GTOPO30]/Categ[CNVTOPO]'
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

    return
  end subroutine CNVTOPO_GTOPO30

  !-----------------------------------------------------------------------------
  !> Convert from DEM 50m mesh
  subroutine CNVTOPO_DEM50M
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_grid_real, only: &
       REAL_LATY, &
       REAL_LONX
    implicit none

    character(len=H_LONG) :: DEM50M_IN_CATALOGUE = ''      !< metadata files for DEM50M
    character(len=H_LONG) :: DEM50M_IN_DIR       = ''      !< directory contains DEM50M files (GrADS format)

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
    real(SP) :: TILE_HEIGHT(1600,1600)
    real(RP) :: TILE_LATH  (0:1600)
    real(RP) :: TILE_LONH  (0:1600)
    real(RP) :: TILE_DLAT, TILE_DLON
    real(RP) :: area, area_fraction

    integer  :: iloc   (1600)
    integer  :: jloc   (1600)
    real(RP) :: ifrac_l(1600) ! fraction for iloc
    real(RP) :: jfrac_b(1600) ! fraction for jloc


    real(RP) :: DOMAIN_LATS, DOMAIN_LATE
    real(RP) :: DOMAIN_LONS, DOMAIN_LONE
    real(RP) :: area_sum(IA,JA)

    character(len=H_LONG) :: fname

    real(RP) :: zerosw
    logical  :: hit_lat, hit_lon
    integer  :: index
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, t
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DEM50M]/Categ[CNVTOPO]'
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

    do j = JS, JE
    do i = IS, IE
       area_sum (i,j) = 0.0_RP
       TOPO_Zsfc(i,j) = 0.0_RP
    enddo
    enddo

    DOMAIN_LATS = REAL_LATY(IS,JS-1) / D2R ! [rad->deg]
    DOMAIN_LATE = REAL_LATY(IE,JE  ) / D2R ! [rad->deg]
    DOMAIN_LONS = REAL_LONX(IS-1,JS) / D2R ! [rad->deg]
    DOMAIN_LONE = REAL_LONX(IE  ,JE) / D2R ! [rad->deg]

    !---< READ from external files >---

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
          if( IO_L ) write(*,*) 'xxx catalogue file not found!'
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

       if (      ( TILE_LATS(t) >= DOMAIN_LATS  .AND. TILE_LATS(t) < DOMAIN_LATE  ) &
            .OR. ( TILE_LATE(t) >= DOMAIN_LATS  .AND. TILE_LATE(t) < DOMAIN_LATE  ) ) then
          hit_lat = .true.
       endif

       if (      ( DOMAIN_LATS  >= TILE_LATS(t) .AND. DOMAIN_LATS  < TILE_LATE(t) ) &
            .OR. ( DOMAIN_LATE  >= TILE_LATS(t) .AND. DOMAIN_LATE  < TILE_LATE(t) ) ) then
          hit_lat = .true.
       endif

       if (      ( TILE_LONS(t) >= DOMAIN_LONS  .AND. TILE_LONS(t) < DOMAIN_LONE  ) &
            .OR. ( TILE_LONE(t) >= DOMAIN_LONS  .AND. TILE_LONE(t) < DOMAIN_LONE  ) ) then
          hit_lon = .true.
       endif

       if (      ( DOMAIN_LONS  >= TILE_LONS(t) .AND. DOMAIN_LONS  < TILE_LONE(t) ) &
            .OR. ( DOMAIN_LONE  >= TILE_LONS(t) .AND. DOMAIN_LONE  < TILE_LONE(t) ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(DEM50M_IN_DIR)//'/'//trim(TILE_fname(t))

          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '+++ Input data      file :', trim(fname)
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** Domain(LAT):', DOMAIN_LATS, DOMAIN_LATE
          if( IO_L ) write(IO_FID_LOG,*) '***       (LON):', DOMAIN_LONS, DOMAIN_LONE
          if( IO_L ) write(IO_FID_LOG,*) '*** Tile  (LAT):', TILE_LATS(t), TILE_LATE(t)
          if( IO_L ) write(IO_FID_LOG,*) '***       (LON):', TILE_LONS(t), TILE_LONE(t)

          TILE_DLAT = 5.0_RP / 60.0_RP / 200.0_RP * D2R
          TILE_DLON = 7.5_RP / 60.0_RP / 200.0_RP * D2R
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLAT  :', TILE_DLAT
          if( IO_L ) write(IO_FID_LOG,*) '*** TILE_DLON  :', TILE_DLON

          fid = IO_get_available_fid()
          open( fid,                    &
                file   = trim(fname),   &
                form   = 'unformatted', &
                access = 'direct',      &
                status = 'old',         &
                recl   = 1600*1600*4,   &
                iostat = ierr           )

             if ( ierr /= 0 ) then
                if( IO_L ) write(*,*) 'xxx data file not found!'
                call PRC_MPIstop
             endif

             read(fid,rec=1) TILE_HEIGHT(:,:)
          close(fid)

          TILE_LATH(0) = TILE_LATS(t) * D2R
          do jj = 1, 1600
             TILE_LATH(jj) = TILE_LATH(jj-1) + TILE_DLAT
!             if( IO_L ) write(IO_FID_LOG,*) jj, TILE_LATH(jj)
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, 1600
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
!             if( IO_L ) write(IO_FID_LOG,*) ii, TILE_LONH(ii)
          enddo

          ! match and calc fraction
          i = IS
          do jj = 1, 1600
             jloc   (jj) = 1 ! Z_sfc(1,1) is used for dummy grid
             jfrac_b(jj) = 1.0_RP

             do j = JS, JE
                if (       TILE_LATH(jj-1) >= REAL_LATY(i,j-1) &
                     .AND. TILE_LATH(jj-1) <  REAL_LATY(i,j  ) ) then

                   jloc   (jj) = j
                   jfrac_b(jj) = min( REAL_LATY(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT

                endif
             enddo
          enddo

          j = JS
          do ii = 1, 1600
             iloc   (ii) = 1 ! Z_sfc(1,1) is used for dummy grid
             ifrac_l(ii) = 1.0_RP

             do i = IS, IE
                if (       TILE_LONH(ii-1) >= REAL_LONX(i-1,j) &
                     .AND. TILE_LONH(ii-1) <  REAL_LONX(i  ,j) ) then

                   iloc   (ii) = i
                   ifrac_l(ii) = min( REAL_LONX(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                endif
             enddo
          enddo

          do jj = 1, 1600
          do ii = 1, 1600
             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) )
!             if( IO_L ) write(IO_FID_LOG,*) ii, jj, area, iloc(ii), jloc(jj), ifrac_l(ii), jfrac_b(jj), TILE_HEIGHT(ii,jj)

             TILE_HEIGHT(ii,jj) = max( TILE_HEIGHT(ii,jj), 0.0 )

             area_fraction = (       ifrac_l(ii)) * (       jfrac_b(jj)) * area
             area_sum (iloc(ii)  ,jloc(jj)  ) = area_sum (iloc(ii)  ,jloc(jj)  ) + area_fraction
             TOPO_Zsfc(iloc(ii)  ,jloc(jj)  ) = TOPO_Zsfc(iloc(ii)  ,jloc(jj)  ) + area_fraction * TILE_HEIGHT(ii,jj)

             area_fraction = (1.0_RP-ifrac_l(ii)) * (       jfrac_b(jj)) * area
             area_sum (iloc(ii)+1,jloc(jj)  ) = area_sum (iloc(ii)+1,jloc(jj)  ) + area_fraction
             TOPO_Zsfc(iloc(ii)+1,jloc(jj)  ) = TOPO_Zsfc(iloc(ii)+1,jloc(jj)  ) + area_fraction * TILE_HEIGHT(ii,jj)

             area_fraction = (       ifrac_l(ii)) * (1.0_RP-jfrac_b(jj)) * area
             area_sum (iloc(ii)  ,jloc(jj)+1) = area_sum (iloc(ii)  ,jloc(jj)+1) + area_fraction
             TOPO_Zsfc(iloc(ii)  ,jloc(jj)+1) = TOPO_Zsfc(iloc(ii)  ,jloc(jj)+1) + area_fraction * TILE_HEIGHT(ii,jj)

             area_fraction = (1.0_RP-ifrac_l(ii)) * (1.0_RP-jfrac_b(jj)) * area
             area_sum (iloc(ii)+1,jloc(jj)+1) = area_sum (iloc(ii)+1,jloc(jj)+1) + area_fraction
             TOPO_Zsfc(iloc(ii)+1,jloc(jj)+1) = TOPO_Zsfc(iloc(ii)+1,jloc(jj)+1) + area_fraction * TILE_HEIGHT(ii,jj)
          enddo
          enddo

       endif
    enddo ! tile loop

    do j = JS, JE
    do i = IS, IE
       zerosw = 0.5_RP - sign( 0.5_RP, area_sum(i,j)-EPS )
       TOPO_Zsfc(i,j) = TOPO_Zsfc(i,j) * ( 1.0_RP-zerosw ) / ( area_sum(i,j)-zerosw )
    enddo
    enddo

    return
  end subroutine CNVTOPO_DEM50M

  !-----------------------------------------------------------------------------
  !> Convert from DEM 50m mesh
  subroutine CNVTOPO_GMTED2010
    use scale_process, only: &
       PRC_MPIstop
    use scale_topography, only: &
       TOPO_Zsfc
    implicit none

    character(len=H_LONG) :: TOPO_GMTED2010_IN_CATALOGUE = ''      !< metadata files for GMTED2010
    character(len=H_LONG) :: TOPO_GMTED2010_IN_DIR       = ''      !< directory contains GMTED2010 files (GrADS format)

    NAMELIST / PARAM_CNVTOPO_GMTED2010 / &
       TOPO_GMTED2010_IN_CATALOGUE, &
       TOPO_GMTED2010_IN_DIR

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[GMTED2010]/Categ[CNVTOPO]'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVTOPO_GMTED2010,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CNVTOPO_GMTED2010. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CNVTOPO_GMTED2010)

    return
  end subroutine CNVTOPO_GMTED2010

  !-----------------------------------------------------------------------------
  !> check slope
  subroutine CNVTOPO_smooth
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_grid, only: &
       GRID_FDX, &
       GRID_FDY
    use scale_statistics, only: &
       STAT_detail
    use scale_topography, only: &
       TOPO_fillhalo, &
       TOPO_Zsfc
    implicit none

    real(RP) :: DZsfc_DX(1,IA,JA,1) ! d(Zsfc)/dx at u-position
    real(RP) :: DZsfc_DY(1,IA,JA,1) ! d(Zsfc)/dy at v-position

    character(len=H_SHORT) :: varname(1)

    integer :: i, j
    !---------------------------------------------------------------------------

    ! digital filter
!    if ( CNVTOPO_smooth_maxslope > UNDEF ) then
!    do j = JS, JE
!    do i = IS, IE
!       TOPO_Zsfc(i,j) = ( TOPO_Zsfc(i  ,j  ) &
!                        + TOPO_Zsfc(i-1,j  ) &
!                        + TOPO_Zsfc(i+1,j  ) &
!                        + TOPO_Zsfc(i  ,j-1) &
!                        + TOPO_Zsfc(i  ,j+1) ) / 5.0_RP
!    enddo
!    enddo
!    endif

    call TOPO_fillhalo

    do j = JS, JE
    do i = IS, IE
       DZsfc_DX(1,i,j,1) = atan2( ( TOPO_Zsfc(i+1,j)-TOPO_Zsfc(i,j) ), GRID_FDX(i) ) / D2R
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       DZsfc_DY(1,i,j,1) = atan2( ( TOPO_Zsfc(i,j+1)-TOPO_Zsfc(i,j) ), GRID_FDY(j) ) / D2R
    enddo
    enddo

    varname(1) = "DZsfc_DX"
    call STAT_detail( DZsfc_DX(:,:,:,:), varname(:) )
    varname(1) = "DZsfc_DY"
    call STAT_detail( DZsfc_DY(:,:,:,:), varname(:) )

    return
  end subroutine CNVTOPO_smooth

end module mod_cnvtopo
