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

  integer,  private :: CNVTOPO_smooth_power_factor   = 4
  integer,  private :: CNVTOPO_smooth_itelim         = 10000
  real(RP), private :: CNVTOPO_smooth_maxslope       = -1.0_RP
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
       CNVTOPO_smooth_power_factor,  &
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
    use mod_gm_statistics, only: &
       GTL_max, &
       GTL_min
    implicit none

    real(RP) :: mintopo, maxtopo
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

       mintopo = GTL_min(TOPOGRAPHY_Zsfc(:,:,:,1),TOPOGRAPHY_Zsfc_pl(:,:,:,1),1,1,1)
       maxtopo = GTL_max(TOPOGRAPHY_Zsfc(:,:,:,1),TOPOGRAPHY_Zsfc_pl(:,:,:,1),1,1,1)

       LOG_NEWLINE
       LOG_INFO_CONT('(A,ES24.16)') 'Minimum topography height = ', mintopo
       LOG_INFO_CONT('(A,ES24.16)') 'Maximum topography height = ', maxtopo
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

    character(len=H_LONG)  :: GTOPO30_IN_DIR       = '.' !< directory contains GTOPO30 files (GrADS format)
    character(len=H_LONG)  :: GTOPO30_IN_CATALOGUE = ''  !< metadata files for GTOPO30
    character(len=H_SHORT) :: GTOPO30_INTERP_TYPE  = 'LINEAR'
    integer                :: GTOPO30_INTERP_LEVEL = 5

    namelist / PARAM_CNVTOPO_GTOPO30 / &
       GTOPO30_IN_DIR,       &
       GTOPO30_IN_CATALOGUE, &
       GTOPO30_INTERP_TYPE,  &
       GTOPO30_INTERP_LEVEL

    ! GTOPO30 data
    real(RP), parameter   :: GTOPO30_DLAT = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.
    real(RP), parameter   :: GTOPO30_DLON = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_GTOPO30",*) 'Setup input data'

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
                          GTOPO30_INTERP_TYPE,  & ! [IN]
                          GTOPO30_INTERP_LEVEL  ) ! [IN]

    LOG_NEWLINE
    LOG_INFO("CNVTOPO_GTOPO30",*) 'Convert from GTOPO30'

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
    integer                :: USERFILE_INTERP_LEVEL   = 5

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

       call CNV2D_tile_init( USERFILE_DTYPE,       & ! [IN]
                             USERFILE_DLAT,        & ! [IN]
                             USERFILE_DLON,        & ! [IN]
                             USERFILE_DIR,         & ! [IN]
                             USERFILE_CATALOGUE,   & ! [IN]
                             USERFILE_INTERP_TYPE, & ! [IN]
                             USERFILE_INTERP_LEVEL ) ! [IN]

    case('GrADS')

       if ( USERFILE_GrADS_FILENAME == '' ) then
          LOG_ERROR("CNVTOPO_USERFILE",*) 'GrADS file name does not specified. Check!'
          call PRC_abort
       endif

       call CNV2D_grads_init( USERFILE_GrADS_FILENAME, & ! [IN]
                              USERFILE_GrADS_VARNAME,  & ! [IN]
                              USERFILE_GrADS_LATNAME,  & ! [IN]
                              USERFILE_GrADS_LONNAME,  & ! [IN]
                              USERFILE_INTERP_TYPE,    & ! [IN]
                              USERFILE_INTERP_LEVEL    ) ! [IN]

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
    use scale_const, only: &
       EPS => CONST_EPS, &
       PI  => CONST_PI
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_have_pl
    use scale_comm_icoA, only: &
       COMM_data_transfer
    use mod_gmtr, only: &
       GMTR_area,   &
       GMTR_area_pl
    use mod_oprt, only: &
       OPRT_gradient,     &
       OPRT_coef_grad,    &
       OPRT_coef_grad_pl, &
       OPRT_laplacian,    &
       OPRT_coef_lap,     &
       OPRT_coef_lap_pl,  &
       OPRT_diffusion,    &
       OPRT_coef_intp,    &
       OPRT_coef_intp_pl, &
       OPRT_coef_diff,    &
       OPRT_coef_diff_pl
    use mod_gm_statistics, only: &
       GTL_max, &
       GTL_min
    implicit none

    real(RP), intent(inout) :: topo   (ADM_gall   ,ADM_KNONE,ADM_lall   ,1)
    real(RP), intent(inout) :: topo_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,1)

    integer,  parameter :: nfact   = 10
    real(RP), parameter :: gamma_h = 1.0_RP / 16.0_RP / 10.0_RP

    real(RP) :: coef        (ADM_gall   ,ADM_KNONE,ADM_lall   )
    real(RP) :: coef_pl     (ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP) :: zgrad_abs_ij(ADM_gall   ,ADM_KNONE,ADM_lall   ,1)
    real(RP) :: dv_ij       (ADM_gall   ,ADM_KNONE,ADM_lall   ,1)

    real(RP) :: zsfc        (ADM_iall,ADM_jall,ADM_KNONE,ADM_lall   )
    real(RP) :: zsfc_pl     (ADM_gall_pl      ,ADM_KNONE,ADM_lall_pl)
    real(RP) :: zgrad       (ADM_iall,ADM_jall,ADM_KNONE,ADM_lall   ,ADM_nxyz)
    real(RP) :: zgrad_pl    (ADM_gall_pl      ,ADM_KNONE,ADM_lall_pl,ADM_nxyz)
    real(RP) :: zgrad_abs   (ADM_iall,ADM_jall,ADM_KNONE,ADM_lall   ,1)
    real(RP) :: zgrad_abs_pl(ADM_gall_pl      ,ADM_KNONE,ADM_lall_pl,1)
    real(RP) :: kh          (ADM_iall,ADM_jall,ADM_KNONE,ADM_lall   )
    real(RP) :: kh_pl       (ADM_gall_pl      ,ADM_KNONE,ADM_lall_pl)
    real(RP) :: dv          (ADM_iall,ADM_jall,ADM_KNONE,ADM_lall   ,1)
    real(RP) :: dv_pl       (ADM_gall_pl      ,ADM_KNONE,ADM_lall_pl,1)
    real(RP) :: maxcoef, mincoef
    real(RP) :: maxgrad

    integer  :: n, ite
    !---------------------------------------------------------------------------

    if ( CNVTOPO_smooth_type == 'OFF' ) then
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Do not apply smoothing.'

       return
    else
       LOG_NEWLINE
       LOG_INFO("CNVTOPO_smooth",*) 'Apply LAPLACIAN smoothing.'
       LOG_INFO_CONT(*) 'Slope limit = ', CNVTOPO_smooth_maxslope_limit

       coef(:,ADM_KNONE,:) = gamma_h * ( 2.0_RP * sqrt(GMTR_area(:,:)/PI) )**2

       if ( PRC_have_pl ) then
          coef_pl(:,ADM_KNONE,:) = gamma_h * ( 2.0_RP * sqrt(GMTR_area_pl(:,:)/PI) )**2
       endif

       mincoef = GTL_min(coef(:,:,:),coef_pl(:,:,:),1,1,1)
       maxcoef = GTL_max(coef(:,:,:),coef_pl(:,:,:),1,1,1)
       LOG_INFO_CONT(*) '### coef (min/max) = ', mincoef, maxcoef

       do ite = 1, CNVTOPO_smooth_itelim

          zsfc = reshape(topo,shape(zsfc))

          if ( PRC_have_pl ) then
             zsfc_pl(:,:,:) = topo_pl(:,:,:,1)
          endif

          call OPRT_gradient( zgrad         (:,:,:,:,:), zgrad_pl         (:,:,:,:), & ! [OUT]
                              zsfc          (:,:,:,:),   zsfc_pl          (:,:,:),   & ! [IN]
                              OPRT_coef_grad(:,:,:,:,:), OPRT_coef_grad_pl(:,:,:)    ) ! [IN]

          zgrad_abs(:,:,:,:,1) = sqrt( zgrad(:,:,:,:,1)**2 &
                                     + zgrad(:,:,:,:,2)**2 &
                                     + zgrad(:,:,:,:,3)**2 )

          if ( PRC_have_pl ) then
             zgrad_abs_pl(:,:,:,1) = sqrt( zgrad_pl(:,:,:,1)**2 &
                                         + zgrad_pl(:,:,:,2)**2 &
                                         + zgrad_pl(:,:,:,3)**2 )
          endif

          zgrad_abs_ij = reshape(zgrad_abs,shape(zgrad_abs_ij))

          call COMM_data_transfer( zgrad_abs_ij, zgrad_abs_pl )

          maxgrad = GTL_max(zgrad_abs_ij(:,:,:,1),zgrad_abs_pl(:,:,:,1),1,1,1)

          zgrad_abs = reshape(zgrad_abs_ij,shape(zgrad_abs))

          kh(:,:,:,:) = ( abs(zgrad_abs(:,:,:,:,1)/maxgrad)+EPS )**CNVTOPO_smooth_power_factor

          if ( PRC_have_pl ) then
             kh_pl(:,:,:) = ( abs(zgrad_abs_pl(:,:,:,1)/maxgrad)+EPS )**CNVTOPO_smooth_power_factor
          endif

          do n = 1, nfact

             zsfc = reshape(topo,shape(zsfc))

             if ( PRC_have_pl ) then
                zsfc_pl(:,:,:) = topo_pl(:,:,:,1)
             endif

             if ( CNVTOPO_smooth_power_factor /= 0.0_RP ) then
                call OPRT_diffusion( dv            (:,:,:,:,1),   dv_pl            (:,:,:,1), & ! [OUT]
                                     zsfc          (:,:,:,:),     zsfc_pl          (:,:,:),   & ! [IN]
                                     kh            (:,:,:,:),     kh_pl            (:,:,:),   & ! [IN]
                                     OPRT_coef_intp(:,:,:,:,:,:), OPRT_coef_intp_pl(:,:,:,:), & ! [IN]
                                     OPRT_coef_diff(:,:,:,:,:),   OPRT_coef_diff_pl(:,:,:)    ) ! [IN]
             else
                call OPRT_laplacian( dv           (:,:,:,:,1), dv_pl           (:,:,:,1), & ! [OUT]
                                     zsfc         (:,:,:,:),   zsfc_pl         (:,:,:),   & ! [IN]
                                     OPRT_coef_lap(:,:,:,:),   OPRT_coef_lap_pl(:,:)      ) ! [IN]
             endif

             dv_ij = reshape(dv,shape(dv_ij))

             call COMM_data_transfer( dv_ij, dv_pl )

             topo(:,:,:,1) = max( topo(:,:,:,1) + coef(:,:,:) * dv_ij(:,:,:,1), 0.0_RP )

             if ( PRC_have_pl ) then
                topo_pl(:,:,:,1) = max( topo_pl(:,:,:,1) + coef_pl(:,:,:) * dv_pl(:,:,:,1), 0.0_RP )
             endif

             call COMM_data_transfer( topo, topo_pl )

!            minzsfc = GTL_min(zsfc(:,:,:,1),zsfc_pl(:,:,:,1),1,1,1)
!            maxzsfc = GTL_max(zsfc(:,:,:,1),zsfc_pl(:,:,:,1),1,1,1)
!            write(IO_FID_LOG,*) '### zsfc (min/max) = ', minzsfc, maxzsfc
          enddo

          LOG_INFO_CONT('(A,I6,A,ES24.16)') 'ite = ', ite, ', maxgrad = ', maxgrad

          if( maxgrad <= CNVTOPO_smooth_maxslope_limit ) exit
          if( maxgrad >= 1.E+5     ) exit
       enddo

       LOG_INFO("CNVTOPO_smooth",*) 'Finish smoothing.'

    endif

    return
  end subroutine CNVTOPO_smooth

end module mod_gm_cnvtopo
