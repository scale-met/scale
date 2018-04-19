!-------------------------------------------------------------------------------
!> module Convert 2D data
!!
!! @par Description
!!          subroutines for preparing 2D data (convert from external file)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cnv2d
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
  public :: CNV2D_setup
  public :: CNV2D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: CNV2D_DoNothing
  logical, public :: CNV2D_UseGrADS = .true.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNV2D_GrADS_bilinear
  private :: CNV2D_GrADS_areaweighted

  private :: CNV2D_write

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP),               private :: CNV2D_unittile_ddeg       =  0.0_RP                ! dx for unit tile [deg]
  real(RP),               private :: CNV2D_oversampling_factor =  2.0_RP                ! factor of min. dx against the unit tile

  character(len=H_SHORT), private :: CNV2D_interpolation_type  = 'bilinear'             ! interpolation type
                                                               ! 'bilinear'               bi-linear interpolation (coarse input->fine model grid)
                                                               ! 'areaweighted'           area-weighted average   (fine input->coarse model grid)
                                                               ! 'nearestneighbor'        nearest neighbor        (both)

  character(len=H_LONG),  private :: CNV2D_OUT_BASENAME        = ''                     ! basename of the output file
  character(len=H_MID),   private :: CNV2D_OUT_TITLE           = 'SCALE-RM 2D Boundary' ! title    of the output file
  character(len=H_SHORT), private :: CNV2D_OUT_VARNAME         = ''                     ! name  of the variable
  character(len=H_MID),   private :: CNV2D_OUT_VARDESC         = ''                     ! title of the variable
  character(len=H_SHORT), private :: CNV2D_OUT_VARUNIT         = ''                     ! units of the variable
  character(len=H_SHORT), private :: CNV2D_OUT_DTYPE           = 'DEFAULT'              ! REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNV2D_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       D2R  => CONST_D2R
    use scale_statistics, only: &
       STATISTICS_horizontal_min
    use scale_atmos_grid_cartesC_real, only: &
       REAL_DLAT => ATMOS_GRID_CARTESC_REAL_DLAT, &
       REAL_DLON => ATMOS_GRID_CARTESC_REAL_DLON
    implicit none

    namelist / PARAM_CNV2D / &
       CNV2D_UseGrADS,            &
       CNV2D_interpolation_type,  &
       CNV2D_unittile_ddeg,       &
       CNV2D_oversampling_factor, &
       CNV2D_OUT_BASENAME,        &
       CNV2D_OUT_TITLE,           &
       CNV2D_OUT_VARNAME,         &
       CNV2D_OUT_VARDESC,         &
       CNV2D_OUT_VARUNIT,         &
       CNV2D_OUT_DTYPE

    real(RP) :: drad(IA,JA)
    real(RP) :: drad_min

    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNV2D_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNV2D,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNV2D_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNV2D_setup",*) 'Not appropriate names in namelist PARAM_CNV2D. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNV2D)

    CNV2D_DoNothing = .true.

    if ( CNV2D_UseGrADS ) then
       CNV2D_DoNothing = .false.
       LOG_INFO("CNV2D_setup",*) 'Use GrADS format file for input'
    endif

    if    ( CNV2D_interpolation_type == 'bilinear'        ) then
       LOG_INFO("CNV2D_setup",*) 'Interpolation type : bi-linear interpolation (coarse input->fine model grid)'
    elseif( CNV2D_interpolation_type == 'areaweighted'    ) then
       LOG_INFO("CNV2D_setup",*) 'Interpolation type : area-weighted average   (fine input->coarse model grid)'
    elseif( CNV2D_interpolation_type == 'nearestneighbor' ) then
       LOG_INFO("CNV2D_setup",*) 'Interpolation type : nearest neighbor'
    else
       LOG_ERROR("CNV2D_setup",*) 'Not appropriate interpolation type. Check!', trim(CNV2D_interpolation_type)
       call PRC_abort
    endif

    if ( CNV2D_DoNothing ) then
       LOG_INFO("CNV2D_setup",*) 'Do nothing for 2D data conversion'
    else
       drad(:,:) = min( REAL_DLAT(:,:), REAL_DLON(:,:) )
       call STATISTICS_horizontal_min( IA, IS, IE, JA, JS, JE, &
                                       drad(:,:), drad_min )

       if ( CNV2D_unittile_ddeg > 0.0_RP ) then
          CNV2D_oversampling_factor = ( drad_min / D2R ) / CNV2D_unittile_ddeg
       endif
       CNV2D_oversampling_factor = max( 1.0_RP, CNV2D_oversampling_factor )
       CNV2D_unittile_ddeg       = ( drad_min / D2R ) / CNV2D_oversampling_factor

       LOG_INFO("CNV2D_setup",*) 'The size of tile [deg] = ', CNV2D_unittile_ddeg
       LOG_INFO("CNV2D_setup",*) 'oversampling factor    = ', CNV2D_oversampling_factor
    endif

    return
  end subroutine CNV2D_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNV2D
    implicit none
    !---------------------------------------------------------------------------

    if ( CNV2D_DoNothing ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'skip  convert topography data'
    else
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start convert topography data'

       if ( CNV2D_UseGrADS ) then
          call CNV2D_GrADS
       endif

       LOG_PROGRESS(*) 'end   convert topography data'
    endif

    return
  end subroutine CNV2D

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNV2D_GrADS
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       PI    => CONST_PI,    &
       D2R   => CONST_D2R
    use scale_calendar, only: &
       CALENDAR_unit2sec
    use scale_atmos_grid_cartesC_real, only: &
       REAL_LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       REAL_LONUY => ATMOS_GRID_CARTESC_REAL_LONUY
    implicit none

    integer                :: GrADS_NLAT         = -1        ! number of latitude  tile
    integer                :: GrADS_NLON         = -1        ! number of longitude tile
    real(RP)               :: GrADS_DLAT         = -1.0_RP   ! width  of latitude  tile [deg.]
    real(RP)               :: GrADS_DLON         = -1.0_RP   ! width  of longitude tile [deg.]
    character(len=H_LONG)  :: GrADS_IN_DIR       = '.'       ! directory contains data files (GrADS format)
    character(len=H_LONG)  :: GrADS_IN_CATALOGUE = ''        ! catalogue file
    character(len=H_LONG)  :: GrADS_IN_FILENAME  = ''        ! single data file (GrADS format)
    character(len=H_LONG)  :: GrADS_IN_DATATYPE  = 'REAL4'   ! datatype (REAL4,REAL8,INT2)
    logical                :: GrADS_LATORDER_N2S = .false.   ! data of the latitude direction is stored in ordar of North->South?
    real(RP)               :: GrADS_MISSINGVAL               ! missing value
    real(RP)               :: GrADS_LAT_START    =  -90.0_RP ! (for single file) start latitude  of domain in input data
    real(RP)               :: GrADS_LAT_END      =   90.0_RP ! (for single file) end   latitude  of domain in input data
    real(RP)               :: GrADS_LON_START    =    0.0_RP ! (for single file) start longitude of domain in input data
    real(RP)               :: GrADS_LON_END      =  360.0_RP ! (for single file) end   longitude of domain in input data
    integer                :: GrADS_NSTEP        = 1         ! number of steps
    real(DP)               :: GrADS_DT           =    0.0_DP ! time interval
    character(len=H_SHORT) :: GrADS_DT_UNIT      = "SEC"     ! time unit for GrADS_DT

    namelist / PARAM_CNV2D_GrADS / &
       GrADS_NLAT,         &
       GrADS_NLON,         &
       GrADS_DLAT,         &
       GrADS_DLON,         &
       GrADS_IN_CATALOGUE, &
       GrADS_IN_DIR,       &
       GrADS_IN_FILENAME,  &
       GrADS_IN_DATATYPE,  &
       GrADS_LATORDER_N2S, &
       GrADS_MISSINGVAL,   &
       GrADS_LAT_START,    &
       GrADS_LAT_END,      &
       GrADS_LON_START,    &
       GrADS_LON_END,      &
       GrADS_NSTEP,        &
       GrADS_DT,           &
       GrADS_DT_UNIT

    real(RP)              :: VAR2D(IA,JA) !< 2D data array
    real(DP)              :: VAR2D_DTSEC

    real(RP)              :: REAL_LONUY_mod(0:IA,JA)
    real(RP)              :: DOMAIN_LATS, DOMAIN_LATE
    real(RP)              :: DOMAIN_LONS, DOMAIN_LONE
    integer               :: DOMAIN_LONSLOC(2), DOMAIN_LONELOC(2)
    logical               :: check_IDL

    integer               :: TILE_NLAT, TILE_NLON
    real(RP)              :: TILE_DLAT, TILE_DLON

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 1000
    integer               :: TILE_nmax
    character(len=H_LONG) :: TILE_fname(TILE_nlim)
    real(RP)              :: TILE_LATS (TILE_nlim)
    real(RP)              :: TILE_LATE (TILE_nlim)
    real(RP)              :: TILE_LONS (TILE_nlim)
    real(RP)              :: TILE_LONE (TILE_nlim)

    character(len=H_LONG) :: fname
    integer               :: t, index
    integer               :: nowstep
    integer               :: fid, ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNV2D_GrADS",*) 'Setup'

    GrADS_MISSINGVAL = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNV2D_GrADS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNV2D_GrADS",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNV2D_GrADS",*) 'Not appropriate names in namelist PARAM_CNV2D_GrADS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNV2D_GrADS)

    if ( GrADS_NLAT <= 0 ) then
       LOG_ERROR("CNV2D_GrADS",*) 'GrADS_NLAT (number of latitude tile)  should be positive. Check!', GrADS_NLAT
       call PRC_abort
    endif

    if ( GrADS_NLON <= 0 ) then
       LOG_ERROR("CNV2D_GrADS",*) 'GrADS_NLON (number of longitude tile) should be positive. Check!', GrADS_NLON
       call PRC_abort
    endif

    if ( GrADS_DLAT <= 0.0_RP ) then
       LOG_ERROR("CNV2D_GrADS",*) 'GrADS_DLAT (width (deg.) of latitude tile) should be positive. Check!', GrADS_DLAT
       call PRC_abort
    endif

    if ( GrADS_DLON <= 0.0_RP ) then
       LOG_ERROR("CNV2D_GrADS",*) 'GrADS_DLON (width (deg.) of longitude tile) should be positive. Check!', GrADS_DLON
       call PRC_abort
    endif

    if (       GrADS_IN_CATALOGUE == '' &
         .AND. GrADS_IN_FILENAME  == '' ) then
       LOG_ERROR("CNV2D_GrADS",*) 'Neither catalogue file nor single file do not specified. Check!'
       call PRC_abort
    endif

    if    ( GrADS_IN_DATATYPE == 'REAL8' ) then
       LOG_INFO("CNV2D_GrADS",*) 'type of input data : REAL8'
    elseif( GrADS_IN_DATATYPE == 'REAL4' ) then
       LOG_INFO("CNV2D_GrADS",*) 'type of input data : REAL4'
    elseif( GrADS_IN_DATATYPE == 'INT2'  ) then
       LOG_INFO("CNV2D_GrADS",*) 'type of input data : INT2'
    else
       LOG_ERROR("CNV2D_GrADS",*) 'Not appropriate type for GrADS_IN_DATATYPE. Check!'
       LOG_ERROR_CONT(*) 'REAL8, REAL4, INT2 are available. requested:', trim(GrADS_IN_DATATYPE)
       call PRC_abort
    endif

    if    ( GrADS_LATORDER_N2S ) then
       LOG_INFO("CNV2D_GrADS",*) 'data ordar of the latitude direction : North -> South'
    else
       LOG_INFO("CNV2D_GrADS",*) 'data ordar of the latitude direction : South -> North'
    endif

    LOG_INFO("CNV2D_GrADS",*) 'Number of steps : ', GrADS_NSTEP
    if ( GrADS_NSTEP > 1 ) then
       call CALENDAR_unit2sec( VAR2D_DTSEC,  GrADS_DT,  GrADS_DT_UNIT  )
       LOG_INFO("CNV2D_GrADS",*) 'Time interval   : ', VAR2D_DTSEC
    else
       VAR2D_DTSEC = 0.0_DP
    endif



    REAL_LONUY_mod(:,:) = mod( REAL_LONUY(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI ! [-> 0..180,-180..0]

    DOMAIN_LATS    = minval(REAL_LATXV    (:,:))
    DOMAIN_LATE    = maxval(REAL_LATXV    (:,:))
    DOMAIN_LONS    = minval(REAL_LONUY_mod(:,:))
    DOMAIN_LONE    = maxval(REAL_LONUY_mod(:,:))
    DOMAIN_LONSLOC = minloc(REAL_LONUY_mod(:,:))
    DOMAIN_LONELOC = maxloc(REAL_LONUY_mod(:,:))

    check_IDL = .false.
    if (      DOMAIN_LONS < REAL_LONUY_mod( 0,DOMAIN_LONSLOC(2)) &
         .OR. DOMAIN_LONE > REAL_LONUY_mod(IA,DOMAIN_LONELOC(2)) ) then
       check_IDL = .true.
       DOMAIN_LONS = minval(REAL_LONUY_mod(:,:),mask=(REAL_LONUY_mod>0.0_RP))
       DOMAIN_LONE = maxval(REAL_LONUY_mod(:,:),mask=(REAL_LONUY_mod<0.0_RP))
    endif

    TILE_NLAT = GrADS_NLAT
    TILE_NLON = GrADS_NLON
    LOG_INFO("CNV2D_GrADS",*) 'Size of data in each tile (j) = ', TILE_NLAT
    LOG_INFO("CNV2D_GrADS",*) 'Size of data in each tile (i) = ', TILE_NLON

    TILE_DLAT = GrADS_DLAT * D2R
    TILE_DLON = GrADS_DLON * D2R
    LOG_INFO("CNV2D_GrADS",*) 'TILE_DLAT       :', TILE_DLAT/D2R
    LOG_INFO("CNV2D_GrADS",*) 'TILE_DLON       :', TILE_DLON/D2R

    !---< READ from external files >---

    if ( GrADS_IN_CATALOGUE /= '' ) then

       ! input from catalogue file
       fname = trim(GrADS_IN_DIR)//'/'//trim(GrADS_IN_CATALOGUE)

       LOG_NEWLINE
       LOG_INFO("CNV2D_GrADS",*) '+ Input catalogue file:', trim(fname)

       fid = IO_get_available_fid()
       open( fid,                  &
             file   = trim(fname), &
             form   = 'formatted', &
             status = 'old',       &
             iostat = ierr         )

          if ( ierr /= 0 ) then
             LOG_ERROR("CNV2D_GrADS",*) 'catalogue file not found!', trim(fname)
             call PRC_abort
          endif

          do t = 1, TILE_nlim
             read(fid,*,iostat=ierr) index, TILE_LATS(t), TILE_LATE(t), & ! South->North
                                            TILE_LONS(t), TILE_LONE(t), & ! WEST->EAST
                                            TILE_fname(t)
             if ( ierr /= 0 ) exit
          enddo

          TILE_nmax = t - 1
       close(fid)

    elseif( GrADS_IN_FILENAME /= '' ) then

       ! input from single file
       TILE_nmax     = 1
       TILE_fname(1) = GrADS_IN_FILENAME
       TILE_LATS (1) = GrADS_LAT_START
       TILE_LATE (1) = GrADS_LAT_END
       TILE_LONS (1) = GrADS_LON_START
       TILE_LONE (1) = GrADS_LON_END

    endif

    do t = 1, TILE_nmax
       if ( TILE_LONS(t) >= 180.0_RP ) then
          TILE_LONS(t) = TILE_LONS(t) - 360.0_RP
          TILE_LONE(t) = TILE_LONE(t) - 360.0_RP
       endif
       if ( TILE_LONS(t) < -180.0_RP ) TILE_LONS(t) = TILE_LONS(t) + 360.0_RP
       if ( TILE_LONE(t) < -180.0_RP ) TILE_LONE(t) = TILE_LONE(t) + 360.0_RP
    enddo

    do nowstep = 1, GrADS_NSTEP

       LOG_INFO("CNV2D_GrADS",*) 'step = ', nowstep

       VAR2D(:,:) = 0.0_RP

       if (      CNV2D_interpolation_type == 'bilinear'        &
            .OR. CNV2D_interpolation_type == 'nearestneighbor' ) then

          call CNV2D_GrADS_bilinear    ( DOMAIN_LATS,             & ! [IN]
                                         DOMAIN_LATE,             & ! [IN]
                                         DOMAIN_LONS,             & ! [IN]
                                         DOMAIN_LONE,             & ! [IN]
                                         TILE_nmax,               & ! [IN]
                                         GrADS_IN_DIR,            & ! [IN]
                                         TILE_fname(1:TILE_nmax), & ! [IN]
                                         TILE_LATS (1:TILE_nmax), & ! [IN]
                                         TILE_LATE (1:TILE_nmax), & ! [IN]
                                         TILE_LONS (1:TILE_nmax), & ! [IN]
                                         TILE_LONE (1:TILE_nmax), & ! [IN]
                                         TILE_NLAT,               & ! [IN]
                                         TILE_NLON,               & ! [IN]
                                         TILE_DLAT,               & ! [IN]
                                         TILE_DLON,               & ! [IN]
                                         GrADS_IN_DATATYPE,       & ! [IN]
                                         check_IDL,               & ! [IN]
                                         GrADS_LATORDER_N2S,      & ! [IN]
                                         GrADS_MISSINGVAL,        & ! [IN]
                                         nowstep,                 & ! [IN]
                                         VAR2D(:,:)               ) ! [OUT]

       elseif( CNV2D_interpolation_type == 'areaweighted' ) then

          call CNV2D_GrADS_areaweighted( DOMAIN_LATS,             & ! [IN]
                                         DOMAIN_LATE,             & ! [IN]
                                         DOMAIN_LONS,             & ! [IN]
                                         DOMAIN_LONE,             & ! [IN]
                                         TILE_nmax,               & ! [IN]
                                         GrADS_IN_DIR,            & ! [IN]
                                         TILE_fname(1:TILE_nmax), & ! [IN]
                                         TILE_LATS (1:TILE_nmax), & ! [IN]
                                         TILE_LATE (1:TILE_nmax), & ! [IN]
                                         TILE_LONS (1:TILE_nmax), & ! [IN]
                                         TILE_LONE (1:TILE_nmax), & ! [IN]
                                         TILE_NLAT,               & ! [IN]
                                         TILE_NLON,               & ! [IN]
                                         TILE_DLAT,               & ! [IN]
                                         TILE_DLON,               & ! [IN]
                                         GrADS_IN_DATATYPE,       & ! [IN]
                                         check_IDL,               & ! [IN]
                                         GrADS_LATORDER_N2S,      & ! [IN]
                                         GrADS_MISSINGVAL,        & ! [IN]
                                         nowstep,                 & ! [IN]
                                         VAR2D(:,:)               ) ! [OUT]

       endif

       ! output 2D file
       call CNV2D_write( VAR2D(:,:),  & ! [IN]
                         VAR2D_DTSEC, & ! [IN]
                         nowstep      ) ! [IN]
    enddo

    return
  end subroutine CNV2D_GrADS

  !-----------------------------------------------------------------------------
  !> Convert from GrADS format file with bi-linear or nearest neighbor interpolation
  subroutine CNV2D_GrADS_bilinear( &
       DOMAIN_LATS,       &
       DOMAIN_LATE,       &
       DOMAIN_LONS,       &
       DOMAIN_LONE,       &
       TILE_nmax,         &
       TILE_dir,          &
       TILE_fname,        &
       TILE_LATS,         &
       TILE_LATE,         &
       TILE_LONS,         &
       TILE_LONE,         &
       TILE_NLAT,         &
       TILE_NLON,         &
       TILE_DLAT,         &
       TILE_DLON,         &
       TILE_DATATYPE,     &
       TILE_check_IDL,    &
       TILE_LATORDER_N2S, &
       TILE_MISSINGVAL,   &
       nowstep,           &
       VAR2D              )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_atmos_grid_cartesC_real, only: &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT, &
       REAL_LON => ATMOS_GRID_CARTESC_REAL_LON
    implicit none

    real(RP),              intent(in)  :: DOMAIN_LATS
    real(RP),              intent(in)  :: DOMAIN_LATE
    real(RP),              intent(in)  :: DOMAIN_LONS
    real(RP),              intent(in)  :: DOMAIN_LONE
    integer,               intent(in)  :: TILE_nmax
    character(len=H_LONG), intent(in)  :: TILE_dir
    character(len=H_LONG), intent(in)  :: TILE_fname(TILE_nmax)
    real(RP),              intent(in)  :: TILE_LATS (TILE_nmax)
    real(RP),              intent(in)  :: TILE_LATE (TILE_nmax)
    real(RP),              intent(in)  :: TILE_LONS (TILE_nmax)
    real(RP),              intent(in)  :: TILE_LONE (TILE_nmax)
    integer,               intent(in)  :: TILE_NLAT
    integer,               intent(in)  :: TILE_NLON
    real(RP),              intent(in)  :: TILE_DLAT
    real(RP),              intent(in)  :: TILE_DLON
    character(len=H_LONG), intent(in)  :: TILE_DATATYPE
    logical,               intent(in)  :: TILE_check_IDL
    logical,               intent(in)  :: TILE_LATORDER_N2S
    real(RP),              intent(in)  :: TILE_MISSINGVAL
    integer,               intent(in)  :: nowstep
    real(RP),              intent(out) :: VAR2D(IA,JA)

    real(RP)   :: REAL_LON_mod(IA,JA)

    real(RP)   :: TILE_LAT     (0:TILE_NLAT+1)
    real(RP)   :: TILE_LON     (0:TILE_NLON+1)
    real(RP)   :: TILE_VALUE   (TILE_NLON,TILE_NLAT)
    real(DP)   :: TILE_VALUE_DP(TILE_NLON,TILE_NLAT)
    real(SP)   :: TILE_VALUE_SP(TILE_NLON,TILE_NLAT)
    integer(2) :: TILE_VALUE_I2(TILE_NLON,TILE_NLAT)

    integer    :: jloc_b
    integer    :: jloc_t
    integer    :: iloc_l
    integer    :: iloc_r
    real(RP)   :: jfrac_t ! fraction for jloc_t
    real(RP)   :: ifrac_r ! fraction for iloc_r

    character(len=H_LONG) :: fname

    logical  :: hit_lat, hit_lon
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, t
    !---------------------------------------------------------------------------

    REAL_LON_mod(:,:) = mod( REAL_LON(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI ! [-> 0..180,-180..0]

    ! data file
    do t = 1, TILE_nmax
       hit_lat = .false.
       hit_lon = .false.

       if (      ( TILE_LATS(t)*D2R >= DOMAIN_LATS .AND. TILE_LATS(t)*D2R < DOMAIN_LATE ) &
            .OR. ( TILE_LATE(t)*D2R >= DOMAIN_LATS .AND. TILE_LATE(t)*D2R < DOMAIN_LATE ) ) then
          hit_lat = .true.
       endif

       if (      ( DOMAIN_LATS >= TILE_LATS(t)*D2R .AND. DOMAIN_LATS < TILE_LATE(t)*D2R ) &
            .OR. ( DOMAIN_LATE >= TILE_LATS(t)*D2R .AND. DOMAIN_LATE < TILE_LATE(t)*D2R ) ) then
          hit_lat = .true.
       endif

       if ( TILE_check_IDL ) then
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS .AND. TILE_LONS(t)*D2R < PI          ) &
               .OR. ( TILE_LONS(t)*D2R >= -PI         .AND. TILE_LONS(t)*D2R < DOMAIN_LONE ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS .AND. TILE_LONE(t)*D2R < PI          ) &
               .OR. ( TILE_LONE(t)*D2R >= -PI         .AND. TILE_LONE(t)*D2R < DOMAIN_LONE ) ) then
             hit_lon = .true.
          endif
       else
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS .AND. TILE_LONS(t)*D2R < DOMAIN_LONE ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS .AND. TILE_LONE(t)*D2R < DOMAIN_LONE ) ) then
             hit_lon = .true.
          endif
       endif

       if (      ( DOMAIN_LONS >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if (      ( DOMAIN_LONS+2.0_RP*PI >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS+2.0_RP*PI < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE+2.0_RP*PI >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE+2.0_RP*PI < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(TILE_dir)//'/'//trim(TILE_fname(t))

          LOG_NEWLINE
          LOG_INFO("CNV2D_GrADS_bilinear",*) '+ Input data file :', trim(fname)
          LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if ( TILE_check_IDL ) then
             LOG_INFO_CONT(*) '(Date line exists within the domain)'
          endif
          LOG_INFO_CONT(*) 'Tile   (LAT)    :', TILE_LATS(t), TILE_LATE(t)
          LOG_INFO_CONT(*) '       (LON)    :', TILE_LONS(t), TILE_LONE(t)

          if ( TILE_DATATYPE == 'REAL8' ) then

             fid = IO_get_available_fid()
             open( fid,                            &
                   file   = trim(fname),           &
                   form   = 'unformatted',         &
                   access = 'direct',              &
                   status = 'old',                 &
                   recl   = TILE_NLON*TILE_NLAT*8, &
                   iostat = ierr                   )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNV2D_GrADS_bilinear",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=nowstep) TILE_VALUE_DP(:,:)
             close(fid)

             if ( TILE_LATORDER_N2S ) then
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_DP(i,TILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_DP(i,j), kind=RP )
                enddo
                enddo
             endif

          elseif( TILE_DATATYPE == 'REAL4' ) then

             fid = IO_get_available_fid()
             open( fid,                            &
                   file   = trim(fname),           &
                   form   = 'unformatted',         &
                   access = 'direct',              &
                   status = 'old',                 &
                   recl   = TILE_NLON*TILE_NLAT*4, &
                   iostat = ierr                   )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNV2D_GrADS_bilinear",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=nowstep) TILE_VALUE_SP(:,:)
             close(fid)

             if ( TILE_LATORDER_N2S ) then
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_SP(i,TILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_SP(i,j), kind=RP )
                enddo
                enddo
             endif

          elseif( TILE_DATATYPE == 'INT2' ) then

             fid = IO_get_available_fid()
             open( fid,                            &
                   file   = trim(fname),           &
                   form   = 'unformatted',         &
                   access = 'direct',              &
                   status = 'old',                 &
                   recl   = TILE_NLON*TILE_NLAT*2, &
                   iostat = ierr                   )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNV2D_GrADS_bilinear",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=nowstep) TILE_VALUE_I2(:,:)
             close(fid)

             if ( TILE_LATORDER_N2S ) then
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_I2(i,TILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_I2(i,j), kind=RP )
                enddo
                enddo
             endif

          endif

          call PROF_Valcheck('CNV2D','VAR',TILE_VALUE(:,:))

          TILE_LAT(0) = TILE_LATS(t) * D2R - TILE_DLAT
          do jj = 1, TILE_NLAT+1
             TILE_LAT(jj) = TILE_LAT(jj-1) + TILE_DLAT
             !LOG_INFO("CNV2D_GrADS_bilinear",*) jj, TILE_LAT(jj) / D2R
          enddo

          TILE_LON(0) = TILE_LONS(t) * D2R - TILE_DLON
          do ii = 1, TILE_NLON+1
             TILE_LON(ii) = TILE_LON(ii-1) + TILE_DLON
             if( TILE_LON(ii) > PI )  TILE_LON(ii) =  TILE_LON(ii) - 2.0_RP*PI
             !LOG_INFO("CNV2D_GrADS_bilinear",*) ii, TILE_LON(ii) / D2R
          enddo

          ! match and calc fraction
          do jj = 1, TILE_NLAT+1
          do ii = 1, TILE_NLON+1
             if (      TILE_LAT(jj  ) < DOMAIN_LATS &
                  .OR. TILE_LAT(jj-1) > DOMAIN_LATE ) then
                cycle
             endif

             if ( TILE_check_IDL .OR. TILE_LON(ii-1) >= TILE_LON(ii) ) then
                ! do nothing
             else
                if (       TILE_LON(ii  ) < DOMAIN_LONS &
                     .OR.  TILE_LON(ii-1) > DOMAIN_LONE ) then
                   cycle
                endif
             endif

             do j = 1, JA
             do i = 1, IA

                iloc_l = -1
                iloc_r = -1
                jloc_b = -1
                jloc_t = -1

                ! normal
                if (       REAL_LON_mod(i,j) >= TILE_LON(ii-1) &
                     .AND. REAL_LON_mod(i,j) <  TILE_LON(ii  ) &
                     .AND. REAL_LAT    (i,j) >= TILE_LAT(jj-1) &
                     .AND. REAL_LAT    (i,j) <  TILE_LAT(jj  ) ) then

                   iloc_l  = ii-1
                   iloc_r  = ii
                   ifrac_r = min( REAL_LON_mod(i,j)-TILE_LON(ii-1), TILE_DLON ) / TILE_DLON

                   jloc_b  = jj-1
                   jloc_t  = jj
                   jfrac_t = min( REAL_LAT(i,j)-TILE_LAT(jj-1), TILE_DLAT ) / TILE_DLAT

                endif

                ! across IDL
                if (       TILE_LON(ii-1) >= TILE_LON(ii)   &
                     .AND. REAL_LAT(i,j)  >= TILE_LAT(jj-1) &
                     .AND. REAL_LAT(i,j)  <  TILE_LAT(jj  ) ) then

                   if    (       REAL_LON_mod(i,j) >= TILE_LON(ii-1) &
                           .AND. REAL_LON_mod(i,j) <  PI             ) then

                      iloc_l  = ii-1
                      iloc_r  = ii
                      ifrac_r = min( REAL_LON_mod(i,j)-TILE_LON(ii-1), TILE_DLON ) / TILE_DLON

                      jloc_b  = jj-1
                      jloc_t  = jj
                      jfrac_t = min( REAL_LAT(i,j)-TILE_LAT(jj-1), TILE_DLAT ) / TILE_DLAT

                   elseif(       REAL_LON_mod(i,j) >= -PI            &
                           .AND. REAL_LON_mod(i,j) <  TILE_LON(ii  ) ) then

                      iloc_l  = ii-1
                      iloc_r  = ii
                      ifrac_r = min( REAL_LON_mod(i,j)+2.0_RP*PI-TILE_LON(ii-1), TILE_DLON ) / TILE_DLON

                      jloc_b  = jj-1
                      jloc_t  = jj
                      jfrac_t = min( REAL_LAT(i,j)-TILE_LAT(jj-1), TILE_DLAT ) / TILE_DLAT

                   endif

                endif

                if( iloc_r == -1 .AND. jloc_t == -1 ) cycle

                if( TILE_LON(0)           <=  0.0_RP    .AND. iloc_l == 0           ) iloc_l = TILE_NLON ! around prime meridian
                if( TILE_LON(TILE_NLON+1) >=  0.0_RP    .AND. iloc_r == TILE_NLON+1 ) iloc_r = 1         ! around prime meridian

                if( TILE_LAT(0)           <= -0.5_RP*PI .AND. jloc_b == 0           ) jloc_b = 1         ! around south pole
                if( TILE_LAT(TILE_NLAT+1) >=  0.5_RP*PI .AND. jloc_t == TILE_NLAT+1 ) jloc_t = TILE_NLAT ! around north pole

                if ( CNV2D_interpolation_type == 'nearestneighbor' ) then
                   if ( ifrac_r >= 0.5D0 ) then
                      iloc_l = iloc_r
                   else
                      iloc_r = iloc_l
                   endif

                   if ( jfrac_t >= 0.5D0 ) then
                      jloc_b = jloc_t
                   else
                      jloc_t = jloc_b
                   endif
                endif

                VAR2D(i,j) = (       ifrac_r) * (       jfrac_t) * TILE_VALUE(iloc_r,jloc_t) &
                           + (1.0_RP-ifrac_r) * (       jfrac_t) * TILE_VALUE(iloc_l,jloc_t) &
                           + (       ifrac_r) * (1.0_RP-jfrac_t) * TILE_VALUE(iloc_r,jloc_b) &
                           + (1.0_RP-ifrac_r) * (1.0_RP-jfrac_t) * TILE_VALUE(iloc_l,jloc_b) + VAR2D(i,j)

             enddo
             enddo

          enddo
          enddo

       endif
    enddo ! tile loop

    return
  end subroutine CNV2D_GrADS_bilinear

  !-----------------------------------------------------------------------------
  !> Convert from GrADS format file with area-weighted average
  subroutine CNV2D_GrADS_areaweighted( &
       DOMAIN_LATS,       &
       DOMAIN_LATE,       &
       DOMAIN_LONS,       &
       DOMAIN_LONE,       &
       TILE_nmax,         &
       TILE_dir,          &
       TILE_fname,        &
       TILE_LATS,         &
       TILE_LATE,         &
       TILE_LONS,         &
       TILE_LONE,         &
       TILE_NLAT,         &
       TILE_NLON,         &
       TILE_DLAT,         &
       TILE_DLON,         &
       TILE_DATATYPE,     &
       TILE_check_IDL,    &
       TILE_LATORDER_N2S, &
       TILE_MISSINGVAL,   &
       nowstep,           &
       VAR2D              )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       EPS    => CONST_EPS,    &
       D2R    => CONST_D2R
    use scale_atmos_grid_cartesC_real, only: &
       REAL_LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       REAL_LONUY => ATMOS_GRID_CARTESC_REAL_LONUY
    implicit none

    real(RP),              intent(in)  :: DOMAIN_LATS
    real(RP),              intent(in)  :: DOMAIN_LATE
    real(RP),              intent(in)  :: DOMAIN_LONS
    real(RP),              intent(in)  :: DOMAIN_LONE
    integer,               intent(in)  :: TILE_nmax
    character(len=H_LONG), intent(in)  :: TILE_dir
    character(len=H_LONG), intent(in)  :: TILE_fname(TILE_nmax)
    real(RP),              intent(in)  :: TILE_LATS (TILE_nmax)
    real(RP),              intent(in)  :: TILE_LATE (TILE_nmax)
    real(RP),              intent(in)  :: TILE_LONS (TILE_nmax)
    real(RP),              intent(in)  :: TILE_LONE (TILE_nmax)
    integer,               intent(in)  :: TILE_NLAT
    integer,               intent(in)  :: TILE_NLON
    real(RP),              intent(in)  :: TILE_DLAT
    real(RP),              intent(in)  :: TILE_DLON
    character(len=H_LONG), intent(in)  :: TILE_DATATYPE
    logical,               intent(in)  :: TILE_check_IDL
    logical,               intent(in)  :: TILE_LATORDER_N2S
    real(RP),              intent(in)  :: TILE_MISSINGVAL
    integer,               intent(in)  :: nowstep
    real(RP),              intent(out) :: VAR2D(IA,JA)

    real(RP)   :: REAL_LONUY_mod(0:IA,JA)

    real(RP)   :: TILE_LATH    (0:TILE_NLAT)
    real(RP)   :: TILE_LONH    (0:TILE_NLON)
    real(RP)   :: TILE_VALUE   (TILE_NLON,TILE_NLAT)
    real(DP)   :: TILE_VALUE_DP(TILE_NLON,TILE_NLAT)
    real(SP)   :: TILE_VALUE_SP(TILE_NLON,TILE_NLAT)
    integer(2) :: TILE_VALUE_I2(TILE_NLON,TILE_NLAT)

    integer    :: jloc, iloc
    real(RP)   :: jfrac_b ! fraction for jloc
    real(RP)   :: ifrac_l ! fraction for iloc

    real(RP)   :: val_sum (IA,JA)
    real(RP)   :: area_sum(IA,JA)
    real(RP)   :: val, mask
    real(RP)   :: area, area_fraction

    character(len=H_LONG) :: fname

    real(RP) :: zerosw
    logical  :: hit_lat, hit_lon
    integer  :: fid, ierr
    integer  :: i, j, ii, jj, t
    !---------------------------------------------------------------------------

    REAL_LONUY_mod(:,:) = mod( REAL_LONUY(:,:)+3.0_DP*PI, 2.0_DP*PI ) - PI ! [-> 0..180,-180..0]

    ! data file
    do t = 1, TILE_nmax
       hit_lat = .false.
       hit_lon = .false.

       if (      ( TILE_LATS(t)*D2R >= DOMAIN_LATS .AND. TILE_LATS(t)*D2R < DOMAIN_LATE ) &
            .OR. ( TILE_LATE(t)*D2R >= DOMAIN_LATS .AND. TILE_LATE(t)*D2R < DOMAIN_LATE ) ) then
          hit_lat = .true.
       endif

       if (      ( DOMAIN_LATS >= TILE_LATS(t)*D2R .AND. DOMAIN_LATS < TILE_LATE(t)*D2R ) &
            .OR. ( DOMAIN_LATE >= TILE_LATS(t)*D2R .AND. DOMAIN_LATE < TILE_LATE(t)*D2R ) ) then
          hit_lat = .true.
       endif

       if ( TILE_check_IDL ) then
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS .AND. TILE_LONS(t)*D2R < PI          ) &
               .OR. ( TILE_LONS(t)*D2R >= -PI         .AND. TILE_LONS(t)*D2R < DOMAIN_LONE ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS .AND. TILE_LONE(t)*D2R < PI          ) &
               .OR. ( TILE_LONE(t)*D2R >= -PI         .AND. TILE_LONE(t)*D2R < DOMAIN_LONE ) ) then
             hit_lon = .true.
          endif
       else
          if (      ( TILE_LONS(t)*D2R >= DOMAIN_LONS .AND. TILE_LONS(t)*D2R < DOMAIN_LONE ) &
               .OR. ( TILE_LONE(t)*D2R >= DOMAIN_LONS .AND. TILE_LONE(t)*D2R < DOMAIN_LONE ) ) then
             hit_lon = .true.
          endif
       endif

       if (      ( DOMAIN_LONS >= TILE_LONS(t)*D2R .AND. DOMAIN_LONS < TILE_LONE(t)*D2R ) &
            .OR. ( DOMAIN_LONE >= TILE_LONS(t)*D2R .AND. DOMAIN_LONE < TILE_LONE(t)*D2R ) ) then
          hit_lon = .true.
       endif

       if ( hit_lat .AND. hit_lon ) then
          fname = trim(TILE_dir)//'/'//trim(TILE_fname(t))

          LOG_NEWLINE
          LOG_INFO("CNV2D_GrADS_areaweighted",*) '+ Input data file :', trim(fname)
          LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
          LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R
          if ( TILE_check_IDL ) then
             LOG_INFO_CONT(*) '(Date line exists within the domain)'
          endif
          LOG_INFO_CONT(*) 'Tile   (LAT)    :', TILE_LATS(t), TILE_LATE(t)
          LOG_INFO_CONT(*) '       (LON)    :', TILE_LONS(t), TILE_LONE(t)

          if ( TILE_DATATYPE == 'REAL8' ) then

             fid = IO_get_available_fid()
             open( fid,                            &
                   file   = trim(fname),           &
                   form   = 'unformatted',         &
                   access = 'direct',              &
                   status = 'old',                 &
                   recl   = TILE_NLON*TILE_NLAT*8, &
                   iostat = ierr                   )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNV2D_GrADS_areaweighted",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=nowstep) TILE_VALUE_DP(:,:)
             close(fid)

             if ( TILE_LATORDER_N2S ) then
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_DP(i,TILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_DP(i,j), kind=RP )
                enddo
                enddo
             endif

          elseif( TILE_DATATYPE == 'REAL4' ) then

             fid = IO_get_available_fid()
             open( fid,                            &
                   file   = trim(fname),           &
                   form   = 'unformatted',         &
                   access = 'direct',              &
                   status = 'old',                 &
                   recl   = TILE_NLON*TILE_NLAT*4, &
                   iostat = ierr                   )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNV2D_GrADS_areaweighted",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=nowstep) TILE_VALUE_SP(:,:)
             close(fid)

             if ( TILE_LATORDER_N2S ) then
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_SP(i,TILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_SP(i,j), kind=RP )
                enddo
                enddo
             endif

          elseif( TILE_DATATYPE == 'INT2' ) then

             fid = IO_get_available_fid()
             open( fid,                            &
                   file   = trim(fname),           &
                   form   = 'unformatted',         &
                   access = 'direct',              &
                   status = 'old',                 &
                   recl   = TILE_NLON*TILE_NLAT*2, &
                   iostat = ierr                   )

                if ( ierr /= 0 ) then
                   LOG_ERROR("CNV2D_GrADS_areaweighted",*) 'data file not found!'
                   call PRC_abort
                endif

                read(fid,rec=nowstep) TILE_VALUE_I2(:,:)
             close(fid)

             if ( TILE_LATORDER_N2S ) then
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_I2(i,TILE_NLAT-j+1), kind=RP ) ! reverse order
                enddo
                enddo
             else
                do j = 1, TILE_NLAT
                do i = 1, TILE_NLON
                   TILE_VALUE(i,j) = real( TILE_VALUE_I2(i,j), kind=RP )
                enddo
                enddo
             endif

          endif

          call PROF_Valcheck('CNV2D','VAR',TILE_VALUE(:,:))

          TILE_LATH(0) = TILE_LATS(t) * D2R
          do jj = 1, TILE_NLAT
             TILE_LATH(jj) = TILE_LATH(jj-1) + TILE_DLAT
             !LOG_INFO("CNV2D_GrADS_areaweighted",*) jj, TILE_LATH(jj) / D2R
          enddo

          TILE_LONH(0) = TILE_LONS(t) * D2R
          do ii = 1, TILE_NLON
             TILE_LONH(ii) = TILE_LONH(ii-1) + TILE_DLON
             if( TILE_LONH(ii) > PI )  TILE_LONH(ii) =  TILE_LONH(ii) - 2.0_RP*PI
             !LOG_INFO("CNV2D_GrADS_areaweighted",*) ii, TILE_LONH(ii) / D2R
          enddo

          ! match and calc fraction
          do jj = 1, TILE_NLAT
          do ii = 1, TILE_NLON

             iloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             ifrac_l = 1.0_RP

             jloc    = 1 ! Z_sfc(1,1) is used for dummy grid
             jfrac_b = 1.0_RP

             if (      TILE_LATH(jj  ) < DOMAIN_LATS &
                  .OR. TILE_LATH(jj-1) > DOMAIN_LATE ) then
                cycle
             endif

             if ( TILE_check_IDL ) then
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
                if (       TILE_LONH(ii-1) >= REAL_LONUY_mod(i-1,j  ) &
                     .AND. TILE_LONH(ii-1) <  REAL_LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH(jj-1) >= REAL_LATXV    (i  ,j-1) &
                     .AND. TILE_LATH(jj-1) <  REAL_LATXV    (i  ,j  ) ) then

                   iloc    = i
                   ifrac_l = min( REAL_LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                   jloc    = j
                   jfrac_b = min( REAL_LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                   exit jloop

                endif

                if (       REAL_LONUY_mod(i-1,j) >= REAL_LONUY_mod(i  ,j  ) &
                     .AND. TILE_LATH     (jj-1)  >= REAL_LATXV    (i  ,j-1) &
                     .AND. TILE_LATH     (jj-1)  <  REAL_LATXV    (i  ,j  ) ) then ! across the IDL

                   if    (       TILE_LONH(ii-1) >= REAL_LONUY_mod(i-1,j) &
                           .AND. TILE_LONH(ii-1) <  PI                    ) then

                      iloc    = i
                      ifrac_l = min( REAL_LONUY_mod(i,j)-TILE_LONH(ii-1)+2.0_RP*PI, TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( REAL_LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   elseif(       TILE_LONH(ii-1) >= -PI                   &
                           .AND. TILE_LONH(ii-1) <  REAL_LONUY_mod(i  ,j) ) then

                      iloc    = i
                      ifrac_l = min( REAL_LONUY_mod(i,j)-TILE_LONH(ii-1), TILE_DLON ) / TILE_DLON

                      jloc    = j
                      jfrac_b = min( REAL_LATXV(i,j)-TILE_LATH(jj-1), TILE_DLAT ) / TILE_DLAT
                      exit jloop

                   endif

                endif
             enddo iloop
             enddo jloop

             if( iloc == 1 .AND. jloc == 1 ) cycle

             val  = TILE_VALUE(ii,jj)
             mask = 0.0_RP

             if ( abs(val-TILE_MISSINGVAL) < EPS ) mask = 1.0_RP! if missing value, mask = 1

             area = RADIUS * RADIUS * TILE_DLON * ( sin(TILE_LATH(jj))-sin(TILE_LATH(jj-1)) ) * ( 1.0_RP - mask )

!             LOG_INFO("CNV2D_GrADS_areaweighted",*) ii, jj, area, iloc, jloc, ifrac_l, jfrac_b, TILE_VALUE(ii,jj)

             area_fraction = (       ifrac_l) * (       jfrac_b) * area
             area_sum(iloc  ,jloc  ) = area_sum(iloc  ,jloc  ) + area_fraction
             val_sum (iloc  ,jloc  ) = val_sum (iloc  ,jloc  ) + area_fraction * val

             area_fraction = (1.0_RP-ifrac_l) * (       jfrac_b) * area
             area_sum(iloc+1,jloc  ) = area_sum(iloc+1,jloc  ) + area_fraction
             val_sum (iloc+1,jloc  ) = val_sum (iloc+1,jloc  ) + area_fraction * val

             area_fraction = (       ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc  ,jloc+1) = area_sum(iloc  ,jloc+1) + area_fraction
             val_sum (iloc  ,jloc+1) = val_sum (iloc  ,jloc+1) + area_fraction * val

             area_fraction = (1.0_RP-ifrac_l) * (1.0_RP-jfrac_b) * area
             area_sum(iloc+1,jloc+1) = area_sum(iloc+1,jloc+1) + area_fraction
             val_sum (iloc+1,jloc+1) = val_sum (iloc+1,jloc+1) + area_fraction * val
          enddo
          enddo

       endif
    enddo ! tile loop

    do j = JS, JE
    do i = IS, IE
       zerosw     = 0.5_RP - sign( 0.5_RP, area_sum(i,j)-EPS )
       VAR2D(i,j) = val_sum(i,j) * ( 1.0_RP-zerosw ) / ( area_sum(i,j)-zerosw )
    enddo
    enddo

    return
  end subroutine CNV2D_GrADS_areaweighted

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine CNV2D_fillhalo( VAR, FILL_BND )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(inout) :: VAR(IA,JA)

    logical,  intent(in), optional :: FILL_BND

    logical :: FILL_BND_
    !---------------------------------------------------------------------------

    FILL_BND_ = .false.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    call COMM_vars8( VAR(:,:), 1 )
    call COMM_wait ( VAR(:,:), 1, FILL_BND_ )

    return
  end subroutine CNV2D_fillhalo

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNV2D_write( &
       VAR2D,     &
       timeintv,  &
       istep      )
    use scale_time, only: &
       TIME_NOWDATE
    use scale_file_cartesC, only: &
       FILE_CARTESC_write
    implicit none

    real(RP), intent(in) :: VAR2D(IA,JA)
    real(DP), intent(in) :: timeintv
    integer,  intent(in) :: istep

    real(RP) :: work(IA,JA,1)
    real(DP) :: timeofs
    !---------------------------------------------------------------------------

    if ( CNV2D_OUT_BASENAME /= '' ) then
       LOG_NEWLINE
       LOG_INFO("CNV2D_write",*) 'Output converted 2D file '

       work(:,:,1) = VAR2D(:,:)

       call CNV2D_fillhalo( work(:,:,1), FILL_BND=.false. )

       timeofs = real(istep-1,kind=DP) * timeintv

       call FILE_CARTESC_write( work(:,:,:),        & ! [IN]
                          CNV2D_OUT_BASENAME, & ! [IN]
                          CNV2D_OUT_TITLE,    & ! [IN]
                          CNV2D_OUT_VARNAME,  & ! [IN]
                          CNV2D_OUT_VARDESC,  & ! [IN]
                          CNV2D_OUT_VARUNIT,  & ! [IN]
                          'XYT',              & ! [IN]
                          CNV2D_OUT_DTYPE,    & ! [IN]
                          timeintv,           & ! [IN]
                          TIME_NOWDATE,       & ! [IN]
                          timeofs=timeofs     ) ! [IN]
    endif

    return
  end subroutine CNV2D_write

end module mod_cnv2d
