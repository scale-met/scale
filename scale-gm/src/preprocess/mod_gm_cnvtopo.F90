!-------------------------------------------------------------------------------
!> module Convert topography
!!
!! @par Description
!!          subroutines for preparing topography data (convert from external file) (GM)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_gm_cnvtopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
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
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNVTOPO_GTOPO30
  private :: CNVTOPO_USERFILE
  private :: CNVTOPO_smooth

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: CNVTOPO_smooth_type = 'LAPLACIAN' ! smoothing type
                                                         ! 'OFF'         Do not apply smoothing
                                                         ! 'LAPLACIAN'   Laplacian filter

  logical,  private :: CNVTOPO_DoNothing
  logical,  private :: CNVTOPO_UseGTOPO30            = .false.
  logical,  private :: CNVTOPO_UseUSERFILE           = .false.

  integer,  private :: CNVTOPO_smooth_hypdiff_order  = 4
  integer,  private :: CNVTOPO_smooth_itelim         = 10000
  real(RP), private :: CNVTOPO_smooth_maxslope       = -1.0_RP ! [deg]
  real(RP), private :: CNVTOPO_smooth_maxslope_limit

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNVTOPO_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=H_SHORT) :: CNVTOPO_name = 'NONE' ! keep backward compatibility

    namelist / PARAM_CNVTOPO / &
       CNVTOPO_name,                  &
       CNVTOPO_UseGTOPO30,            &
       CNVTOPO_UseUSERFILE,           &
       CNVTOPO_smooth_hypdiff_order,  &
       CNVTOPO_smooth_maxslope,       &
       CNVTOPO_smooth_itelim,         &
       CNVTOPO_smooth_type

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_setup",*) 'Setup'

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
       CNVTOPO_UseUSERFILE  = .false.
    case('USERFILE')
       CNVTOPO_UseUSERFILE  = .true.
    case default
       LOG_ERROR("CNVTOPO_setup",*) 'Unsupported TYPE: ', trim(CNVTOPO_name)
       call PRC_abort
    endselect

    CNVTOPO_DoNothing = .true.

    if ( CNVTOPO_UseGTOPO30 ) then
       CNVTOPO_DoNothing = .false.
       LOG_INFO("CNVTOPO_setup",*) 'Use GTOPO, global 30 arcsec. data'
    elseif ( CNVTOPO_UseUSERFILE ) then
       CNVTOPO_DoNothing = .false.
       LOG_INFO("CNVTOPO_setup",*) 'Use user-defined file'
    endif

    if ( CNVTOPO_DoNothing ) then
       LOG_INFO("CNVTOPO_setup",*) 'Do nothing for topography data'
    endif

    if ( CNVTOPO_smooth_maxslope > 0.0_RP ) then

       CNVTOPO_smooth_maxslope_limit = CNVTOPO_smooth_maxslope

    else
       LOG_ERROR("CNVTOPO_setup",*) 'CNVTOPO_smooth_maxslope most be positive: ', CNVTOPO_smooth_maxslope
       call PRC_abort
    endif

    return
  end subroutine CNVTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVTOPO
    use scale_comm_icoA, only: &
       COMM_data_transfer
    use mod_gm_topography, only: &
       TOPOGRAPHY_fillhalo, &
       TOPOGRAPHY_Zsfc,     &
       TOPOGRAPHY_Zsfc_pl
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( CNVTOPO_DoNothing ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'skip  convert topography data'
    else
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start convert topography data'

       if ( CNVTOPO_UseGTOPO30 ) then
          call CNVTOPO_GTOPO30( TOPOGRAPHY_Zsfc   (:,ADM_KNONE,:,1), & ! [INOUT]
                                TOPOGRAPHY_Zsfc_pl(:,ADM_KNONE,:,1)  ) ! [INOUT]
       endif

       if ( CNVTOPO_UseUSERFILE ) then
          call CNVTOPO_USERFILE( TOPOGRAPHY_Zsfc   (:,ADM_KNONE,:,1), & ! [INOUT]
                                 TOPOGRAPHY_Zsfc_pl(:,ADM_KNONE,:,1)  ) ! [INOUT]
       endif

       call CNVTOPO_smooth( TOPOGRAPHY_Zsfc   (:,:,:,:), & ! [INOUT]
                            TOPOGRAPHY_Zsfc_pl(:,:,:,:)  ) ! [INOUT]

       call TOPOGRAPHY_fillhalo

       LOG_PROGRESS(*) 'end   convert topography data'
    endif

    return
  end subroutine CNVTOPO

  !-----------------------------------------------------------------------------
  !> Convert from GTOPO30
  subroutine CNVTOPO_GTOPO30( &
       topo,   &
       topo_pl )
    use scale_prc, only: &
       PRC_abort
    use mod_gm_cnv2d, only: &
       CNV2D_tile_init, &
       CNV2D_convert
    implicit none

    real(RP), intent(inout) :: topo   (ADM_gall   ,ADM_lall   )
    real(RP), intent(inout) :: topo_pl(ADM_gall_pl,ADM_lall_pl)

    character(len=H_LONG) :: GTOPO30_IN_DIR       = '.' !< directory contains GTOPO30 files (GrADS format)
    character(len=H_LONG) :: GTOPO30_IN_CATALOGUE = ''  !< metadata files for GTOPO30

    namelist / PARAM_CNVTOPO_GTOPO30 / &
       GTOPO30_IN_DIR,       &
       GTOPO30_IN_CATALOGUE

    ! GTOPO30 data
    real(RP), parameter   :: GTOPO30_DLAT = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.
    real(RP), parameter   :: GTOPO30_DLON = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.

    integer :: ierr
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

    call CNV2D_tile_init( 'INT2',               & ! [IN]
                          GTOPO30_DLAT,         & ! [IN]
                          GTOPO30_DLON,         & ! [IN]
                          GTOPO30_IN_DIR,       & ! [IN]
                          GTOPO30_IN_CATALOGUE, & ! [IN]
                          'LINEAR'              ) ! [IN]

    call CNV2D_convert( topo   (:,:),           & ! [OUT]
                        topo_pl(:,:),           & ! [OUT]
                        min_value = -9000.0_RP, & ! [IN]
                        yrevers   = .true.      ) ! [IN]

    return
  end subroutine CNVTOPO_GTOPO30

  !-----------------------------------------------------------------------------
  !> Convert from User-defined file
  subroutine CNVTOPO_USERFILE( &
       topo,   &
       topo_pl )
    use scale_prc, only: &
       PRC_abort
    use mod_gm_cnv2d, only: &
       CNV2D_tile_init,  &
       CNV2D_grads_init, &
       CNV2D_convert
    implicit none

    real(RP), intent(inout) :: topo   (ADM_gall   ,ADM_lall   )
    real(RP), intent(inout) :: topo_pl(ADM_gall_pl,ADM_lall_pl)

    character(len=H_SHORT) :: USERFILE_TYPE           = ''       ! "TILE" or "GrADS"
    ! TILE data
    character(len=H_SHORT) :: USERFILE_DTYPE          = 'REAL4'  ! datatype (REAL4,REAL8,INT2,INT4)
    real(RP)               :: USERFILE_DLAT           = -1.0_RP  ! width  of latitude  tile [deg.]
    real(RP)               :: USERFILE_DLON           = -1.0_RP  ! width  of longitude tile [deg.]
    character(len=H_LONG)  :: USERFILE_DIR            = '.'      ! directory contains data files (GrADS format)
    character(len=H_LONG)  :: USERFILE_CATALOGUE      = ''       ! catalogue file
    logical                :: USERFILE_yrevers        = .false.  ! data of the latitude direction is stored in ordar of North->South?
    real(RP)               :: USERFILE_MINVAL         = 0.0_RP
    ! GrADS data
    character(len=H_LONG)  :: USERFILE_GrADS_FILENAME = ''       ! single data file (GrADS format)
    character(len=H_SHORT) :: USERFILE_GrADS_VARNAME  = 'topo'
    character(len=H_SHORT) :: USERFILE_GrADS_LATNAME  = 'lat'
    character(len=H_SHORT) :: USERFILE_GrADS_LONNAME  = 'lon'
    character(len=H_SHORT) :: USERFILE_INTERP_TYPE    = 'LINEAR'
    integer                :: USERFILE_INTERP_level   = 5

    namelist / PARAM_CNVTOPO_USERFILE / &
       USERFILE_TYPE,           &
       USERFILE_DTYPE,          &
       USERFILE_DLAT,           &
       USERFILE_DLON,           &
       USERFILE_CATALOGUE,      &
       USERFILE_DIR,            &
       USERFILE_yrevers,        &
       USERFILE_MINVAL,         &
       USERFILE_GrADS_FILENAME, &
       USERFILE_GrADS_VARNAME,  &
       USERFILE_GrADS_LATNAME,  &
       USERFILE_GrADS_LONNAME,  &
       USERFILE_INTERP_TYPE,    &
       USERFILE_INTERP_LEVEL

    integer :: ierr
    !---------------------------------------------------------------------------

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

    select case(USERFILE_TYPE)
    case('TILE')

       if ( USERFILE_DLAT <= 0.0_RP ) then
          LOG_ERROR("CNVTOPO_USERFILE",*) 'USERFILE_DLAT (width (deg.) of latitude tile) should be positive. Check! ', USERFILE_DLAT
          call PRC_abort
       endif
       if ( USERFILE_DLON <= 0.0_RP ) then
          LOG_ERROR("CNVTOPO_USERFILE",*) 'USERFILE_DLON (width (deg.) of longitude tile) should be positive. Check! ', USERFILE_DLON
          call PRC_abort
       endif
       if ( USERFILE_CATALOGUE == '' ) then
          LOG_ERROR("CNVTOPO_USERFILE",*) 'Catalogue file does not specified. Check!'
          call PRC_abort
       endif

       call CNV2D_tile_init( USERFILE_DTYPE,     & ! [IN]
                             USERFILE_DLAT,      & ! [IN]
                             USERFILE_DLON,      & ! [IN]
                             USERFILE_DIR,       & ! [IN]
                             USERFILE_CATALOGUE, & ! [IN]
                             'LINEAR'            ) ! [IN]

    case('GrADS')

       if ( USERFILE_GrADS_FILENAME == '' ) then
          LOG_ERROR("CNVTOPO_USERFILE",*) 'GrADS file name does not specified. Check!'
          call PRC_abort
       endif

       call CNV2D_grads_init( USERFILE_GrADS_FILENAME,             &
                              USERFILE_GrADS_VARNAME,              &
                              USERFILE_GrADS_LATNAME,              &
                              USERFILE_GrADS_LONNAME,              &
                              USERFILE_INTERP_TYPE,                &
                              interp_level = USERFILE_INTERP_LEVEL )

    case default
       LOG_ERROR("CNVTOPO_USERFILE",*) 'USERFILE_TYPE is invalid: ',trim(USERFILE_TYPE)
       LOG_ERROR_CONT(*) 'It must be "TILE" or "GrADS"'
       call PRC_abort
    end select

    call CNV2D_convert( topo   (:,:),                & ! [OUT]
                        topo_pl(:,:),                & ! [OUT]
                        min_value = USERFILE_MINVAL, & ! [IN]
                        yrevers   = USERFILE_yrevers ) ! [IN]

    return
  end subroutine CNVTOPO_USERFILE

  !-----------------------------------------------------------------------------
  !> check slope
  subroutine CNVTOPO_smooth( &
       topo,   &
       topo_pl )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(inout) :: topo   (ADM_gall   ,ADM_KNONE,ADM_lall   ,1)
    real(RP), intent(inout) :: topo_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,1)

    !---------------------------------------------------------------------------

    if ( CNVTOPO_smooth_type == 'OFF' ) then
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Do not apply smoothing.'

       return
    else
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Apply smoothing.'
       LOG_INFO_CONT(*) 'Slope limit    = ', CNVTOPO_smooth_maxslope_limit
       LOG_INFO_CONT(*) 'Smoothing type = ', CNVTOPO_smooth_type
       LOG_NEWLINE
    endif

    return
  end subroutine CNVTOPO_smooth

end module mod_gm_cnvtopo
